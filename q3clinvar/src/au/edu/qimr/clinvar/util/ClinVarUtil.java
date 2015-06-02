package au.edu.qimr.clinvar.util;

import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.procedure.TIntIntProcedure;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.util.ChrPositionUtils;
import org.qcmg.common.util.Constants;
import org.qcmg.common.util.Pair;
import org.qcmg.common.vcf.VcfRecord;

import au.edu.qimr.clinvar.model.Bin;
import au.edu.qimr.clinvar.model.Probe;

public class ClinVarUtil {
	
	private static final QLogger logger = QLoggerFactory.getLogger(ClinVarUtil.class);
	
	public static int [] getDoubleEditDistance(String read1, String read2, String primer1, String primer2, int editDistanceCutoff) {
		
		int editDistance = getEditDistance(read1, primer1, editDistanceCutoff + 1);
		int editDistance2 = Integer.MAX_VALUE;
		
		if (editDistance <= editDistanceCutoff) {
			// get read2 edit distance
			editDistance2 = getEditDistance(read2, primer2);
		}
		
		return new int [] {editDistance, editDistance2};
	}
	
	public static int  getEditDistance(String read, String primer) {
		if (org.qcmg.common.string.StringUtils.isNullOrEmpty(read)
				|| org.qcmg.common.string.StringUtils.isNullOrEmpty(primer)) {
			throw new IllegalArgumentException("read or primer (or both) supplied to ClinVarUtil.getEditDistance were null. read: " + read + ", primer: " + primer);
		}
		
//		//TODO see if doing a basic edit distance here and only running levenshtein of bed > cutoff would save time
//		int bed = getBasicEditDistance(primer, read.substring(0, primer.length()));
//		if (bed <= 1) {
//			return bed;
//		}
		
		return StringUtils.getLevenshteinDistance(primer, read.substring(0, primer.length()));
	}
	
	public static int  getEditDistance(String read, String primer, int editDistanceCutoff) {
//		if (org.qcmg.common.string.StringUtils.isNullOrEmpty(read)
//				|| org.qcmg.common.string.StringUtils.isNullOrEmpty(primer)) {
//			throw new IllegalArgumentException("read or primer (or both) supplied to ClinVarUtil.getEditDistance were null. read: " + read + ", primer: " + primer);
//		}
		int led = StringUtils.getLevenshteinDistance(primer, read.substring(0, primer.length()), editDistanceCutoff);
		
		return led >= 0 ? led : Integer.MAX_VALUE;
	}
	
	public static String breakdownEditDistanceDistribution(List<Integer> editDistances) {
		TIntIntHashMap dist = new TIntIntHashMap();
		final StringBuilder sb = new StringBuilder();
		for (Integer ed : editDistances) {
			int existingValue = dist.get(ed);
			dist.put(ed, existingValue + 1);
		}
		dist.forEachEntry(new TIntIntProcedure() {

			@Override
			public boolean execute(int arg0, int arg1) {
				sb.append(arg0).append(":").append(arg1).append(",");
				return true;
			}

		});
		return sb.toString();
	}
	
	
	public static int getBasicEditDistance(CharSequence s, CharSequence t) {
		if (StringUtils.isEmpty(s) || StringUtils.isEmpty(t)) {
			throw new IllegalArgumentException("Null string passed to ClinVarUtil.getBasicEditDistance(). s: " + s + ", t: " + t);
		}
		
		// s and t need to be the same length
		if (s.length() != t.length()) {
			throw new IllegalArgumentException("Strings passed to ClinVarUtil.getBasicEditDistance() are not the same length s: " + s + ", t: " + t);
		}
		
		int ed = 0;
		for (int i = 0 , len = s.length() ; i < len ; i++) {
			if (s.charAt(i) != t.charAt(i)) {
				ed ++;
			}
		}
		
		return ed;
	}
	
	public static int getOutwardSlideCount(final String overlap1, final String overlap2) {
		// it is assumed that the 2 char sequences do not match as they are
		if (StringUtils.isEmpty(overlap1) || StringUtils.isEmpty(overlap2)) {
			throw new IllegalArgumentException("Null string passed to ClinVarUtil.getOutwardSlideCount(). s: " + overlap1 + ", t: " + overlap2);
		}
		
		// s and t need to be the same length
		if (overlap1.length() != overlap2.length()) {
			throw new IllegalArgumentException("Strings passed to ClinVarUtil.getOutwardSlideCount() are not the same length s: " + overlap1 + ", t: " + overlap2);
		}
		String s = overlap1;
		String t = overlap2;
		
		int initialLength = overlap1.length();
		int noOfSlides = 0;
		
		while (noOfSlides < initialLength &&  ! t.equals(s)) {
			noOfSlides++;
			s = s.substring(1);
			t = t.substring(0, t.length() -1);
		}
		// want a -ve number for outward slide
		return 0 - noOfSlides;
	}
	
//	public static String getFragmentUsingSlide(final String fullSeq1, final String fullSeq2, int expectedOverlap) {
//		// it is assumed that the 2 char sequences do not match as they are
//		if (StringUtils.isEmpty(fullSeq1) || StringUtils.isEmpty(fullSeq2)) {
//			throw new IllegalArgumentException("Null string passed to ClinVarUtil.getInwardSlideCount(). s: " + fullSeq1 + ", t: " + fullSeq2);
//		}
//		
//		if (fullSeq1.length() <= expectedOverlap || fullSeq2.length() <= expectedOverlap) {
//			throw new IllegalArgumentException("Strings passed to ClinVarUtil.getInwardSlideCount(). s: " + fullSeq1 + ", t: " + fullSeq2 +", are not longer than the expectedOverlap: " + expectedOverlap);
//		}
//		
////		String s = fullSeq1.substring(fullSeq1.length() - expectedOverlap);
//		String t = fullSeq2.substring(0, expectedOverlap);
//		
//		// get index of t occurring in fullSeq1
//		int index = fullSeq1.lastIndexOf(t);
//		int noOfShifts = 0;
//		while (noOfShifts < expectedOverlap &&  index == -1) {
//			noOfShifts++;
//			t = t.substring(0, t.length() -1);
//			index = fullSeq1.lastIndexOf(t);
//		}
//		if (index == -1) {
//			// can't create fragment
//			logger.info("unable to create fragment for fullSeq1: " + fullSeq1 + ", fullSeq2: " + fullSeq2 + ", expectedOverlap: " + expectedOverlap);
//			return null;
//		}
//		
//		index = fullSeq1.length() - expectedOverlap - index;
//		
//		return index;
		
//		int initialLength = fullSeq1.length();
//		int noOfSlides = 0;
//		
//		while (noOfSlides < initialLength &&  ! t.equals(s)) {
//			noOfSlides++;
//			s = s.substring(1);
//			t = t.substring(0, t.length() -1);
//		}
//		// want a -ve number for outward slide
//		return 0 - noOfSlides;
//	}
	
	
	/**
	 * 
	 * @param smithWatermanDiffs
	 * @return
	 */
	public static List<Pair<Integer, String>> getPositionRefAndAltFromSW(String [] smithWatermanDiffs) {
		List<Pair<Integer, String>> mutations = new ArrayList<>();
		
		String refSeq = smithWatermanDiffs[0];
		String diffs = smithWatermanDiffs[1];
		String binSeq = smithWatermanDiffs[2];
		
		
		if (null != diffs && ! diffs.isEmpty()) {
		
			if (diffs.charAt(0) == ' ') {
				logger.warn("First char in diffs string is empty string!!!");
			}
			int position = 0;
			int span = 0;
			int indelStartPosDiff = 0;
			int indelStartPosRef = 0;
			int insertionBaseCount = 0;
			for (char c : diffs.toCharArray()) {
				if (c != ' ') {
					if (span >0) {
						// create indel
						
						int start = Math.max(0, indelStartPosDiff - 1);
						String ref = refSeq.substring(start, indelStartPosDiff + span);
						String alt = binSeq.substring(start, indelStartPosDiff + span);
						mutations.add(new Pair<Integer, String>(indelStartPosRef - 1, ref.replaceAll("-","") +"/" +  alt.replaceAll("-","")));
						// reset span
						span = 0;
					}
					if (c == '.') {
						// snp
						char ref = refSeq.charAt(position);
						char alt = binSeq.charAt(position);
						mutations.add(new Pair<Integer, String>(position - insertionBaseCount, ref + "/" + alt));
					}
					
				} else {
					if (span == 0) {
						indelStartPosRef = position - insertionBaseCount;
						indelStartPosDiff = position;
					}
					span++;
					// indel
					// if this is an insertion, update insertionBaseCount
					if (refSeq.charAt(position) == '-') {
						insertionBaseCount++;
					}
				}
				position++;
			}
			if (span >0) {
				// create indel
				
				int start = Math.max(0, indelStartPosDiff - 1);
				String ref = refSeq.substring(start, indelStartPosDiff + span);
				String alt = binSeq.substring(start, indelStartPosDiff + span);
				mutations.add(new Pair<Integer, String>(indelStartPosRef - 1, ref.replaceAll("-","") +"/" +  alt.replaceAll("-","")));
			}
		}
		return mutations;
	}
	
	public static List<Probe> getAmpliconsOverlappingPosition(ChrPosition cp, Set<Probe> probes) {
		List<Probe> overlappingProbes = new ArrayList<>();
		for (Probe p : probes) {
			if (ChrPositionUtils.isChrPositionContained(p.getCp(), cp)) {
				overlappingProbes.add(p);
			}
		}
		return overlappingProbes;
	}
	
	
	public static String getSortedBBString(String bb, String ref) {
		String [] params = bb.split(";");
		int length = params.length;
		if (length == 1) {
			return bb;
		}
		
		List<String> sizeSortedList = new ArrayList<>();
		for (String s : params) {
			String sRef = s.substring(0, s.indexOf(","));
			if ( ! sRef.equals(ref)) {
				sizeSortedList.add(s);
			}
		}
		
		Collections.sort(sizeSortedList, new Comparator<String>(){

			@Override
			public int compare(String arg0, String arg1) {
				int arg0Tally = Integer.valueOf(arg0.split(",")[1]);
				int arg1Tally = Integer.valueOf(arg1.split(",")[1]);
				return arg1Tally - arg0Tally;
			}
			
		});
		
		// now add in the ref
		for (String s : params) {
			String sRef = s.substring(0, s.indexOf(","));
			if (sRef.equals(ref)) {
				sizeSortedList.add(s);
			}
		}
		
		//stringify
		StringBuilder sb = new StringBuilder();
		for (String s : sizeSortedList) {
			if (sb.length() > 0) {
				sb.append(";");
			}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	/**
	 * Returns a representation of the supplied position as seen in the supplied amplicon and bins in the following format:
	 * base,count, ampliconId(total reads in amplicon),binId1(count),binId2(count).....
	 * @param cp
	 * @param overlappingProbes
	 * @param probeBinDist
	 * @return
	 */
	public static String getAmpliconDistribution(VcfRecord vcf, List<Probe> overlappingProbes, 
			Map<Probe, List<Bin>> probeBinDist, int minBinSize) {
		return getAmpliconDistribution(vcf, overlappingProbes, probeBinDist, minBinSize, false);
	}
	
	public static String getAmpliconDistribution(VcfRecord vcf, List<Probe> overlappingProbes, 
			Map<Probe, List<Bin>> probeBinDist, int minBinSize, boolean diagnosticMode) {
		StringBuilder sb = new StringBuilder();
		
		Map<String, List<Pair<Probe, Bin>>> baseDist = new HashMap<>();
		
		for (Probe amplicon : overlappingProbes) {
			
			List<Bin> bins = probeBinDist.get(amplicon);
			
			for (Bin b : bins) {
				/*
				 * only deal with bins that have >= minBinSize read in them
				 */
				if (b.getRecordCount() >= minBinSize) {
					String binSeq = getBaseAtPosition(vcf, amplicon, b);
					if (null != binSeq) {
						List<Pair<Probe, Bin>> probeBinPair = baseDist.get(binSeq);
						if (null == probeBinPair) {
							probeBinPair = new ArrayList<>();
							baseDist.put(binSeq, probeBinPair);
						}
						probeBinPair.add(new Pair<Probe, Bin>(amplicon, b));
					}
				}
			}
		}
		
		
		// convert map to sb
		for (Entry<String, List<Pair<Probe, Bin>>> entry : baseDist.entrySet()) {
			String bases = entry.getKey();
			List <Pair<Probe, Bin>> probeBinList = entry.getValue();
			int tally = 0;
			for (Pair<Probe, Bin> pair : probeBinList) {
				tally += pair.getRight().getRecordCount();
			}
			
			StringBuilder s = new StringBuilder(bases);
			s.append(Constants.COMMA);
			s.append(tally);
			if (diagnosticMode) {
				for (Pair<Probe, Bin> pair : probeBinList) {
					s.append(Constants.COMMA);
					s.append(pair.getLeft().getId()).append("/");
					s.append(pair.getRight().getId()).append("(").append(pair.getRight().getRecordCount()).append(")");
				}
			} else {
				TIntSet probeSet = new TIntHashSet();
				TIntSet binSet = new TIntHashSet();
				
				for (Pair<Probe, Bin> pair : probeBinList) {
					probeSet.add(pair.getLeft().getId());
					binSet.add(pair.getRight().getId());
				}
				s.append(Constants.COMMA);
				s.append(probeSet.size()).append("/").append(binSet.size());
				
			}
			if (sb.length() > 0) {
				sb.append(Constants.SEMI_COLON);
			}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	public static int getZeroBasedPositionInString(int mutationStartPosition, int binStartPosition) {
		return mutationStartPosition - binStartPosition;
	}
	
	public static String getBaseAtPosition(VcfRecord vcf, Probe amplicon, Bin bin) {
		
		String [] smithWatermanDiffs = bin.getSmithWatermanDiffs();
		if (null == smithWatermanDiffs) {
			logger.warn("bin does not contain sw diffs!!! bin: " + bin.getId() + ", probe: " + amplicon.getId() + ", vcf cp: " + vcf.getChrPosition().toIGVString());
		}
		String probeRef = amplicon.getReferenceSequence();
		// remove any indel characters - only checking
		String swRef = smithWatermanDiffs[0].replace("-", "");
		int offset = probeRef.indexOf(swRef);
		
		if ( offset == -1) {
			logger.warn("probeRef.indexOf(swRef) == -1!!! probe (id:bin.id: " + amplicon.getId() + ":" + bin.getId() + ") , probe ref: " + probeRef + ", swRef: " + swRef);
		}
		
		int length = vcf.getChrPosition().getLength();
		int positionInString = getZeroBasedPositionInString(vcf.getChrPosition().getPosition(), amplicon.getCp().getPosition() + offset);
		
		if (amplicon.getId() == 241) {
			logger.info("positionInString: " + positionInString +", from offset: " + offset + ", vcf.getChrPosition().getPosition(): " + vcf.getChrPosition().getPosition() +", amplicon.getCp().getPosition(): " + amplicon.getCp().getPosition());
		}
		
		return getMutationString(positionInString, length, smithWatermanDiffs);
		
	}
	
	public static String getMutationString(final int positionInString, final int eventLength, String [] smithWatermanDiffs) {
		
		int expectedStart = positionInString;
		int expectedEnd = expectedStart + eventLength;
		int binSequenceLength = smithWatermanDiffs[2].length();
		
		if (expectedStart < 0) {
			/*
			 * This happens when the bin is shorter than (starts after) the amplicon.
			 * In this situation, we need to return null
			 */
//			logger.info("bin " + bin.getId() + ", in amplicon " + amplicon.getId() + " has an offset of " + offset + ", which means that it starts after the position we are interested in " + vcf.getChrPosition().toIGVString());
			return null;
		}
		if (expectedEnd > binSequenceLength) {
			logger.warn("Expected end: " + expectedEnd + ", is greater than the length of the bin sequence: " + binSequenceLength);
			return null;
		}
		
		int additionalOffset = 0;
		String swDiffRegion = smithWatermanDiffs[1].substring(0, expectedStart); 
		if (swDiffRegion.contains(" ")) {
			
			// keep track of number of insertions and deletions
			
			for (int i = 0 ; i < expectedStart ; i++) {
				char c = swDiffRegion.charAt(i);
				if (c == ' ') {
					if ('-' == smithWatermanDiffs[0].charAt(i)) {
						// insertion
						additionalOffset++;
//					} else if ('-' == smithWatermanDiffs[2].charAt(i)) {
//						//deleteion
//						additionalOffset--;
//					} else {
//						//!@##$!!! something else-tion
//						logger.warn("Should have found either an indel at position " + i + " in smithWatermanDiffs[0] or smithWatermanDiffs[2]");
					}
				}
			}
		} else {
			// no additional offset required 
		}
		
		
		/*
		 *  next get the span of the event
		 *  If we have an insertion after the position we are dealing with, get the type and length
		 */
		int span = 0;
		int i = 1;
		while (expectedStart + i < binSequenceLength && smithWatermanDiffs[1].charAt(expectedStart + i) == ' ') {
			if (smithWatermanDiffs[0].charAt(expectedStart + i) == '-') {
				span++;
			}
			i++;
		}
		
		int finalStart = expectedStart + additionalOffset;
		int finalEnd = expectedEnd + additionalOffset + span;
		
		if (finalStart >=finalEnd) {
			logger.warn("finalStart: " + finalStart + " is greater than finalEnd: " + finalEnd);
		} else if (finalStart < 0) {
			logger.warn("finalStart: " + finalStart + " is less than 0");
		} else if (finalEnd < 0) {
			logger.warn("finalEnd: " + finalEnd + " is less than 0");
		} else if (finalEnd > smithWatermanDiffs[1].length()) {
			logger.warn("finalEnd: " + finalEnd + " is > than smithWatermanDiffs[1].length(): " + smithWatermanDiffs[1].length());
		} else {
			return smithWatermanDiffs[2].substring(finalStart, finalEnd).replace("-", "");
		}
		
		return null;
	}
	
	
	
	public static int noOfSlidesToGetPerfectMatch(final String s1, final String t1) {
		// it is assumed that the 2 char sequences do not match as they are
		if (StringUtils.isEmpty(s1) || StringUtils.isEmpty(t1)) {
			throw new IllegalArgumentException("Null string passed to ClinVarUtil.noOfSlidesToGetPerfectMatch(). s: " + s1 + ", t: " + t1);
		}
		
		// s and t need to be the same length
		if (s1.length() != t1.length()) {
			throw new IllegalArgumentException("Strings passed to ClinVarUtil.noOfSlidesToGetPerfectMatch() are not the same length s: " + s1 + ", t: " + t1);
		}
		String s = s1;
		String t = t1;
		
		int initialLength = s1.length();
		int noOfSlides = 0;
		// left slide first
		while (noOfSlides < initialLength &&  ! t.equals(s)) {
			noOfSlides++;
			t = t.substring(1);
			s = s.substring(0, s.length() -1);
		}
		
		
		// need a reliable check to see if noOfSlides is sufficiently large to trigger a RHS slide
		if (noOfSlides >= initialLength -1) {
			//perform a RHS slide
			s = s1;
			t = t1;
			noOfSlides = 0;
			// left slide first
			while (noOfSlides < initialLength &&  ! t.equals(s)) {
				noOfSlides--;
				s = s.substring(1);
				t = t.substring(0, t.length() -1);
			}
		}
		
		return noOfSlides;
	}
	
	
	public static int[] getBasicAndLevenshteinEditDistances(CharSequence s, CharSequence t) {
		// do exact match first
		if (s.equals(t)) {
			return new int[]{0,0};			// we have a match - set both distances to 0
		} else {
			// do basic distancing next
			int ed = ClinVarUtil.getBasicEditDistance(s, t);
			return ed > 1 ? new int[]{ ed, StringUtils.getLevenshteinDistance(s,t) } :  new int[]{ ed, ed };
			
		}
	}
}
