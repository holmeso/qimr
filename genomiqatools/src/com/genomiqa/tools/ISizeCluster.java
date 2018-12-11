package com.genomiqa.tools;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicIntegerArray;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

//import org.apache.commons.math3.fitting.PolynomialCurveFitter;
//import org.apache.commons.math3.fitting.WeightedObservedPoints;
//import org.apache.commons.math3.ml.clustering.Cluster;
//import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
//import org.apache.commons.math3.ml.clustering.DoublePoint;
import org.qcmg.common.model.ChrPosition;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

public class ISizeCluster {
	
	public static final double MIN_PERCENTAGE = 99.5d;
	
	String [] qpXmlFiles;
	String geneModelFile;
	long genomeLength;
	List<ChrPosition> regions = new ArrayList<>();
	Map<String, ChrPosition> geneCPMap;
	
	private void engage() throws IOException {
		try {
			for (String s : qpXmlFiles) {
				loadISizes(s);
			}
		} catch (ParserConfigurationException | SAXException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void loadISizes(String qpXmlFile) throws IOException, ParserConfigurationException, SAXException {
		
		File file = new File(qpXmlFile);
        DocumentBuilderFactory documentBuilderFactory = DocumentBuilderFactory.newInstance();
         
        DocumentBuilder documentBuilder = documentBuilderFactory.newDocumentBuilder();
        Document document = documentBuilder.parse(file);
         
        document.getDocumentElement().normalize();
         
//        System.out.println("Root element :" + document.getDocumentElement().getNodeName());
        NodeList nodeList = document.getElementsByTagName("ISIZE");
//        System.out.println("nodeList length: " + nodeList.getLength());
        Node isizeNode = nodeList.item(0);
        NodeList iSizeReadGroups = isizeNode.getChildNodes();
//        System.out.println("iSizeReadGroups length: " + iSizeReadGroups.getLength());
        for (int i = 0 ; i < iSizeReadGroups.getLength() ; i++) {
        		Node rgNode = iSizeReadGroups.item(i);
//        		System.out.println("rgNode :" + rgNode.getNodeName());
        		if ("ReadGroup".equals(rgNode.getNodeName()) || "RG".equals(rgNode.getNodeName())) {
//        			System.out.println(rgNode.getAttributes().getNamedItem("id"));
        			String readGroupName = getReadGroupName(rgNode);
        			if ( ! "overall".equals(readGroupName)) {
        				System.out.println("rgName: " + readGroupName);
        			
	        			AtomicIntegerArray aia =  getAIAFromNode(rgNode);
	        			getLowerAndUpperBounds(aia, readGroupName);
        			}
        			
        		}
        }
	}
	
	public static int[] getLowerAndUpperBounds(AtomicIntegerArray aia, String rg) {
		
		if ("d1194f81-e28c-41d0-9d12-92c7b53002a8".equals(rg)) {
			StringBuilder sb = new StringBuilder();
			int modal = getModalValue(aia);
			for (int j = 0 ; j < aia.length() ; j++) {
				if (j > modal && aia.get(j) == 0) {
					break;
				}
				sb.append(aia.get(j)).append(",");
			}
			System.out.println(sb.toString());
		}
				
		int [] singlePointValues = singlePointDifferenceMethod(aia, rg);
		int [] percentageMethodValues = percentageMethod(aia, rg, 99.9);
		int [] multiplePointValues = multiplePointDifferenceMethod(aia, rg, 10, 9);
		System.out.println("singlePointValues: " + singlePointValues[0] + "," + singlePointValues[1]);
		System.out.println("percentageMethodValues: " + percentageMethodValues[0] + "," + percentageMethodValues[1]);
		System.out.println("multiplePointValues: " + multiplePointValues[0] + "," + multiplePointValues[1]);
		/*
		 * get the maximum lowerBound and the minimum upperBound
		 */
//		int lb = Math.max(multiplePointValues[0], Math.max(singlePointValues[0], percentageMethodValues[0]));
		int lb = singlePointValues[0];
		if (lb != 18) {
			/*
			 * check to see if we have 18 from the percentageValues method
			 */
			if (18 == percentageMethodValues[0]) {
				lb = 18;
			}
			
		}
		int drawTheLineUpperBound = drawTheLineMethod(aia, lb);
		System.out.println("drawTheLineMethod: " + lb + "," + drawTheLineUpperBound);
		
		
//		int ub = Math.min(drawTheLineUpperBound, Math.min(multiplePointValues[1], Math.min(singlePointValues[1], percentageMethodValues[1])));
		/*
		 *if the drawTheLine method is larger than the multiplePoint method, then use that, otherwise use min of the other methods
		 */
//		int ub =  Math.min(multiplePointValues[1], Math.min(singlePointValues[1], percentageMethodValues[1]));
//		if (drawTheLineUpperBound > ub) {
//			ub = drawTheLineUpperBound;
//		}
		
		int ub = drawTheLineUpperBound;
		
		
		double percent = getAreaUndeCurve(aia, lb, ub);
		System.out.println("LB: " + lb + ", UB: " + ub + ", percent: " + percent);
		System.out.println("RG=" + rg + " LB=" + lb + " UB=" + ub);
		return new int[] {lb, ub};
	}
	
	public static int drawTheLineMethod(AtomicIntegerArray aia, int lowerBound) {
		
		int countAtLower = aia.get(lowerBound);
		int modal = getModalValue(aia);
		for (int i = modal + 1, len = aia.length() ; i < len ; i++) {
			int count = aia.get(i);
			if (count <= countAtLower) {
				/*
				 * check that we have met the percentage requirement
				 */
				if (getAreaUndeCurve(aia, lowerBound, i) >= MIN_PERCENTAGE) {
					return i;
				}
			}
		}
		
		return -1;
	}
	
	public static String getReadGroupName(Node n) {
		if ("ReadGroup".equals(n.getNodeName())) {
			return n.getAttributes().getNamedItem("id").getNodeValue();
		} else if ("RG".equals(n.getNodeName())) {
			return n.getAttributes().getNamedItem("value").getNodeValue();
		}
		return null;
	}
	
	public static AtomicIntegerArray getAIAFromNode(Node n) {
		AtomicIntegerArray aia = new AtomicIntegerArray(5000);
		NodeList rtItem = n.getChildNodes();
		
		if ("RG".equals(n.getNodeName())) {
//			System.out.println("rtItem length: " + rtItem.getLength());
			for (int j = 0 ; j < rtItem.getLength() ; j++) {
				Node rtiNode = rtItem.item(j);
				
				NodeList nl = rtiNode.getChildNodes();
				for (int k = 0 ; k < nl.getLength() && k < 5000; k++) {
				
					Node node = nl.item(k);
					
					if (null != node && node.getNodeType() == Node.ELEMENT_NODE) {
						Element eElement = (Element) node;
						int start = Integer.parseInt( eElement.getAttribute("start"));
						int end = Integer.parseInt( eElement.getAttribute("end"));
//						System.out.println("start: " + start + ", end: " + end);
						if (start == end) {
							int count = Integer.parseInt( eElement.getAttribute("count"));
							aia.getAndSet(start, count);
						}
					}
				}
			}
			
		} else {
		
			for (int j = 0 ; j < rtItem.getLength() ; j++) {
				Node rtiNode = rtItem.item(j);
				
				if (rtiNode.getNodeType() == Node.ELEMENT_NODE) {
					Element eElement = (Element) rtiNode;
//					System.out.println(eElement.toString());
					int start = Integer.parseInt( eElement.getAttribute("start"));
					int end = Integer.parseInt( eElement.getAttribute("end"));
	//				System.out.println("start: " + start + ", end: " + end);
					if (start == end) {
						int count = Integer.parseInt( eElement.getAttribute("count"));
						aia.getAndSet(start, count);
					}
				}
			}
		}
		return aia;
	}
	
	public static int getModalValue(AtomicIntegerArray aia) {
		/*
		 * find modal
		 */
		int modal = 0;
		int modalCount = 0;
		for (int i = 1, len = aia.length() ; i < len ; i++ ) {
			int count = aia.get(i);
			if (count > modalCount) {
				modal = i;
				modalCount = count;
			}
		}
		return modal;
	}
	
	public static int[] percentageMethod(AtomicIntegerArray aia, String rgName, double percentage) {
		int modal = getModalValue(aia);
		int prevCount = aia.get(modal);
		int lowerBound = 0;
		for (int iSize = modal - 1 ; iSize >= 0 ; iSize -- ) {
			 
			 int thisCount = aia.get(iSize);
			 /*
			  * get difference between thisCount and prevCount as a percentage
			  */
			 double diff = (double) thisCount / prevCount;
			 if (diff <= 0.70) {
				 lowerBound = iSize;
//				 System.out.println("we've got one!!!! isize: " + iSize + ", diff: " + diff + ", thisCount: " + thisCount + ", prevCount: " + prevCount);
				 break;
			 }
			 
			 /*
			  * set prevCount to thisCount
			  */
			 prevCount = thisCount;
		}
		
		/*
		 * once lower bound has been calculated, we go forward until a certain percentage has been reached
		 */
		int upperBound = getUpperBoundBasedOnPercentage(aia, lowerBound, percentage);
//		System.out.println("percentage method: RG=" + rgName + " lowerBound=" + lowerBound + " upperBound=" + upperBound);
//		System.out.println("RG=" + rgName + " lowerBound=" + lowerBound + " upperBound=" + upperBound);
		return new int[] {lowerBound, upperBound};
	}
	
	public static double[] getStdDevAndMean(double [] array) {
		double sum = 0.0, stdDev = 0.0;
        int length = array.length;

        for (double num : array) {
            sum += num;
        }

        double mean = sum / length;

        for (double num : array) {
            stdDev += Math.pow(num - mean, 2);
        }

        return new double[] {mean, Math.sqrt(stdDev / length)};
	}
	
	public static int[] multiplePointDifferenceMethod(AtomicIntegerArray aia, String rgName, int windowSize, int stdDevMultiplier) {
		/*
		 * find modal
		 */
		int modal = getModalValue(aia);
		/*
		 * setup array to hold the binned values across the window size for each iSize
		 */
		int [] iSizeWindowCountsLB = new int[ modal];
		double [] percentageDiffsLB = new double[modal];
		for (int iSize = modal -1 ; (iSize - windowSize) > 0 ; iSize --) {
			/*
			 * sum the next windowSize counts and put into array
			 */
			int tally = aia.get(iSize);
			for (int i = 1 ; i < windowSize && (iSize - i ) >= 0 ; i++) {
				tally += aia.get(iSize - i);
			}
			iSizeWindowCountsLB[iSize] = tally;
		}
		
		int previousTally = iSizeWindowCountsLB[iSizeWindowCountsLB.length -1];
		int lowerBound = 0;
		for (int i = iSizeWindowCountsLB.length - 1 ; i >= 0 ; i--) {
			double diff = (100d * iSizeWindowCountsLB[i]) / previousTally;
			if (i <  (iSizeWindowCountsLB.length - 10)) {
				double [] meanAndStdDev = getStdDevAndMean(Arrays.copyOfRange(percentageDiffsLB, i, modal -1));
				double mean = meanAndStdDev[0];
				double stdDev = meanAndStdDev[1];
				double lowerLimit = mean - (stdDevMultiplier * stdDev);
				if (diff < lowerLimit) {
//					if (diff < lowerLimit || diff > upperLimit) {
					lowerBound = i - (windowSize / 2);
					
					System.out.println("diff at position " + (i + modal + (windowSize / 2)) + ": " + diff + " is more than " + stdDevMultiplier + "*stdDevs from the mean! mean: " + mean + ", stdDev: " + stdDev);
					break;
						
				}
			}
			percentageDiffsLB[i] = diff;
			previousTally = 	iSizeWindowCountsLB[i];
		}
		
		
		
		/*
		 * setup array to hold the binned values across the window size for each iSize
		 */
		int [] iSizeWindowCountsUB = new int[aia.length() - modal];
		double [] percentageDiffsUB = new double[aia.length() - modal];
		for (int iSize = modal, len = aia.length() - windowSize ; iSize < len ; iSize ++) {
			/*
			 * sum the next windowSize counts and put into array
			 */
			int tally = aia.get(iSize);
			for (int i = 1 ; i < windowSize ; i++) {
				tally += aia.get(iSize + i);
			}
			iSizeWindowCountsUB[iSize - modal] = tally;
		}
		
		/*
		 * examine array, and find point at which difference between current and previous is significant.....
		 */
		previousTally = iSizeWindowCountsUB[0];
		int upperBound = Integer.MAX_VALUE;
		for (int i = 1 ; i < iSizeWindowCountsUB.length ; i++) {
			double diff = (100d * iSizeWindowCountsUB[i]) / previousTally;
			if (i > 10) {
				double [] meanAndStdDev = getStdDevAndMean(Arrays.copyOfRange(percentageDiffsUB, 1, i));
				double mean = meanAndStdDev[0];
				double stdDev = meanAndStdDev[1];
				double lowerLimit = mean - (stdDevMultiplier * stdDev);
				if (diff < lowerLimit) {
//					if (diff < lowerLimit || diff > upperLimit) {
					upperBound = i + modal + (windowSize / 2);
					
					/*
					 * lets see if this gives us a percentage value greater that the min
					 */
					if (getAreaUndeCurve(aia, lowerBound, upperBound) >= MIN_PERCENTAGE) {
						System.out.println("diff at position " + (i + modal + (windowSize / 2)) + ": " + diff + " is more than " + stdDevMultiplier + "*stdDevs from the mean! mean: " + mean + ", stdDev: " + stdDev);
						break;
						
					} else {
						// keep looking
					}
					
				}
			}
			percentageDiffsUB[i] = diff;
			previousTally = 	iSizeWindowCountsUB[i];
		}
		return new int[] {lowerBound, upperBound};
	}
	
	public static int[] singlePointDifferenceMethod(AtomicIntegerArray aia, String rgName) {
		/*
		 * find modal
		 */
		int modal = getModalValue(aia);
		
		/*
		 * get starting point
		 */
		int prevCount = aia.get(modal);
		
		 /*
		  * look for difference between previous position - flag if count drops by more than a certain percentage (maybe 30%)
		  * Start at modal and work down both sides of mountain
		  */
		 prevCount = aia.get(modal);
		 int lowerBound = 0;
		 int upperBound = Integer.MAX_VALUE;
		 double [] diffValues = new double[modal];
		 for (int iSize = modal - 1 ; iSize >= 0 ; iSize -- ) {
			 
			 int thisCount = aia.get(iSize);
			 /*
			  * get difference between thisCount and prevCount as a percentage
			  */
			 double diff = (double) thisCount / prevCount;
			 diffValues[iSize] = diff;
//			 if (diff <= 0.70) {
//				 lowerBound = iSize;
//				 System.out.println("we've got one!!!! isize: " + iSize + ", diff: " + diff + ", thisCount: " + thisCount + ", prevCount: " + prevCount);
//				 break;
//			 }
			 
			 /*
			  * set prevCount to thisCount
			  */
			 prevCount = thisCount;
		}
		 
		 /*
		  * loop through diffValues, and the one with the lowest value wins!
		  */
		 double minValue = Double.MAX_VALUE;
		 int minValueISize = modal;
		 for (int i = diffValues.length -1 ; i > 0  ; i--) {
			 double d = diffValues[i];
			 if (d > 0.0 && d < minValue) {
				 minValue = d;
				 minValueISize = i;
			 }
		 }
		 lowerBound = minValueISize;
		 
		 
		 /*
		  * now for the other slope
		  */
//		 diffValues = new double[aia.length() - modal];
		 prevCount = aia.get(modal);
		 boolean ignoreNextISize = false;
		 for (int iSize = modal + 1 ; iSize < aia.length() ; iSize ++ ) {
			 int thisCount = aia.get(iSize);
			 
			 double diff = (double) thisCount / prevCount;
			 /*
			  * if we have had an upward spike, set ignoreNextISize to true as we don't want the presence of a spike to skew this very scientific method
			  */
			 if (ignoreNextISize) {
				 ignoreNextISize = false;
			 }
			 else if (diff > 1.1) {
//				 System.out.println("got a spike! at iSize: " + iSize);
				 ignoreNextISize = true;
			 }
			 /*
			  * get difference between thisCount and prevCount as a percentage
			  */
			 else if (diff <= 0.85) {
				 
				 /*
				  * is the count of this iSize less that that of the lower bound?
				  */
				 if (thisCount <= aia.get(lowerBound)) {
					 
					 /*
					  * check to see if percentage is good
					  */
					 double perc = getAreaUndeCurve(aia, lowerBound, iSize); 
					 if (perc >= 99.75) {
					 
						 upperBound = iSize;
						 System.out.println("we've got one!!!! isize: " + iSize + ", diff: " + diff + ", thisCount: " + thisCount + ", prevCount: " + prevCount);
						 break;
					 } else {
						 System.out.println("Min coverage level not yet reached: " + perc);
					 }
				 } else {
					 System.out.println("lower bound count not yet reached: " + thisCount + ", lowerBound count: " +  aia.get(lowerBound));
				 }
			 }
			 
//			 diffValues[iSize - modal] = diff;
			 
			 /*
			  * set prevCount to thisCount
			  */
			 prevCount = thisCount;
		 }
		 
//		 /*
//		  * loop through diffValues, and the one with the lowest value wins!
//		  */
//		 minValue = Double.MAX_VALUE;
//		 minValueISize = modal;
//		 for (int i = 0 ; i <  diffValues.length ; i++) {
//			 double d = diffValues[i];
//			 if (d > 0.0 && d < minValue) {
//				 minValue = d;
//				 minValueISize = i;
//			 }
//		 }
//		 upperBound = minValueISize;
		 
		 
		 
		 /*
		  * Perform some checks
		  *  - area under the curve must be > 99.5%
		  *  - lower bound must be > 1
		  *  - upper bound count must be less than lower bound count
		  */
//		 double percentageCutoff = 99.5;
//		 double percentage = getAreaUndeCurve(aia, lowerBound, upperBound);
//		 if (percentage >= 99.5) {
//			 System.out.println("PASSED percentage check (" + percentage +" is >= " + percentageCutoff + ")");
//		 } else {
//			 System.out.println("FAILED percentage check (" + percentage +" is NOT >= " + percentageCutoff + ")");
//			 
//		 }
//		 if (lowerBound > 1) {
//			 System.out.println("PASSED lower bound check (" + lowerBound +" is > 1)");
// 		 } else {
// 			 System.out.println("FAILED lower bound check (" + lowerBound +" is NOT > 1)");
// 		 }
//		 int lowerBC = aia.get(lowerBound);
//		 int upperBC = aia.get(upperBound);
//		 if (lowerBC > upperBC) {
//			 System.out.println("PASSED lower bound count must be greater than upper bound count check (" + lowerBC +" is > " + upperBC + ")");
//		 } else {
//			 System.out.println("FAILED lower bound count must be greater than upper bound count check (" + lowerBC +" is NOT > " + upperBC + ")");
//		 }
		 
//		 System.out.println("single difference method:" + rgName + " lowerBound=" + lowerBound + " upperBound=" + upperBound);
//		 System.out.println("RG=" + rgName + " lowerBound=" + lowerBound + " upperBound=" + upperBound);
		 return new int[]{lowerBound, upperBound};
		 
	}
	
	public static double getAreaUndeCurve(AtomicIntegerArray aia, int start, int end) {
		int withinISizeTally = 0;
		int totalTally = 0;
		 for (int i = 1 ; i < aia.length() ; i++) {
			 if (i >= start && i <= end) {
				 withinISizeTally += aia.get(i);
			 }
			 totalTally += aia.get(i);
		 }
		 return ((100d*withinISizeTally) / totalTally);
	}
	public static int getUpperBoundBasedOnPercentage(AtomicIntegerArray aia, int start, double percentage) {
		int withinISizeTally = 0;
		int totalTally = 0;
		for (int i = 1 ; i < aia.length() ; i++) {
			totalTally += aia.get(i);
		}
		for (int i = start ; i < aia.length() ; i++) {
			withinISizeTally += aia.get(i);
			if ((100d*withinISizeTally) / totalTally >= percentage) {
				/*
				 * we have found our upper bound
				 */
				return i;
			}
		}
		return Integer.MAX_VALUE;
	}
	
	public static void main(String[] args) {
		ISizeCluster crs = new ISizeCluster();
		crs.qpXmlFiles = args;
//		crs.geneModelFile = args[1];
		try {
			crs.engage();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}