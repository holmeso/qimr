package com.genomiqa.tools;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.model.ChrRangePosition;
import org.qcmg.common.util.ChrPositionUtils;
import org.qcmg.common.util.TabTokenizer;

public class CallingRegionStats {
	
	String regionsFile;
	String geneModelFile;
	List<ChrPosition> regions = new ArrayList<>();
	Map<String, ChrPosition> geneCPMap;
	
	private void engage() throws IOException {
		loadRegionsFile();
		
		loadGeneModel();
		
		examineGeneModel();
	}
	
	/*
	 * TODO must add test, and soon
	 */
	private void examineGeneModel() {
		
		/*
		 * flag the genes that are not covered by the regions
		 */
		int unmatchedTally = 0;
		int overlapTally = 0;
		for (Entry<String, ChrPosition> entry : geneCPMap.entrySet()) {
			ChrPosition cp = entry.getValue();
			
			boolean foundMatch = false;
			boolean overlap = false;
			for (ChrPosition regionCP : regions) {
				if (ChrPositionUtils.isChrPositionContained(regionCP, cp)) {
					foundMatch = true;
					break;
				} else if (ChrPositionUtils.doChrPositionsOverlap(regionCP,  cp)) {
					overlap = true;
					System.out.println("overlap: " + regionCP.toIGVString() + " and " + cp.toIGVString());
				}
			}
			
			if ( ! foundMatch) {
				System.out.println("Could not find region that contains gene " + entry.getKey() + " : " + cp.toIGVString());
				unmatchedTally++;
				if (overlap) {
					overlapTally++;
				}
			}
		}
		
		System.out.println("total number of genes that don't match: " + unmatchedTally + " of which " + overlapTally + " overlap a region");
		
		
	}
	
	private void loadGeneModel() throws IOException {
		System.out.println("about to load gene model data");
		List<String> lines = Files.readAllLines(Paths.get(geneModelFile));
		System.out.println("number of entries in the gene model file: " + lines.size());
		
		Map<Object, List<String>> geneNameMap = lines.stream().filter(s ->  ! s.startsWith("#")).collect(Collectors.groupingBy(s -> s.substring(s.indexOf("gene_id") + 9, s.indexOf("\";"))));
		System.out.println("number of entries in the gene model map: " + geneNameMap.size());
		lines.clear();
		lines = null;
		
		geneCPMap = new HashMap<>();
		for (Entry<Object, List<String>> entry : geneNameMap.entrySet()) {
			
			/*
			 * loop over list, tokenize, and get start, stop and contig
			 */
			List<Integer> numbers = new ArrayList<>();
			String contig = null;
			for (String s : entry.getValue()) {
				String [] array = TabTokenizer.tokenize(s);
				int start = Integer.parseInt(array[3]);
				int end = Integer.parseInt(array[4]);
				if (null == contig) {
					contig = array[0];
				} else {
					/*
					 * check that contigs are the same...
					 */
					if ( ! contig.equals(array[0])) {
						System.out.println("found different contigs for gene: " + entry.getKey());
					}
				}
				numbers.add(start);
				numbers.add(end);
			}
			
			numbers.sort(null);
			// System.out.println("gene " + entry.getKey() + " is in " + contig + " with min start " + numbers.get(0) + " and max " + numbers.get(numbers.size()-1));
			if ( ! contig.startsWith("chr")) {
				contig = "chr" + contig;
			}
			
			geneCPMap.put(entry.getKey().toString(), new ChrRangePosition(contig, numbers.get(0), numbers.get(numbers.size() - 1)));
		}
		
		
	}
	
	private void loadRegionsFile() throws IOException {
		List<String> lines = Files.readAllLines(Paths.get(regionsFile));
		
		/*
		 * go through the headers, adding the contig lengths on the way, and get list of ChrPosition for regions
		 */
		long l = lines.stream()
				.filter(s -> s.startsWith("@SQ")).mapToLong(s -> Long.parseLong(s.substring(s.indexOf("LN:") + 3, s.indexOf("\t", s.indexOf("LN:")))))
				.sum();
		regions = lines.stream()
				.filter(s -> ! s.startsWith("@"))
				.map(s -> {
					String[]array = TabTokenizer.tokenize(s);
					return new ChrRangePosition(array[0], Integer.parseInt(array[1]), Integer.parseInt(array[2]));
					})
				.collect(Collectors.toList());
		
		System.out.println("genome length: " + l);
		System.out.println("number of regions: " + regions.size());
		long regionsLength = regions.stream().mapToLong(cp -> cp.getLength()).sum();
		System.out.println("number of bases covered by regions: " + regionsLength + " which is " + ((regionsLength * 100) / l) + "%");
	}
	
	
	public static void main(String[] args) {
		CallingRegionStats crs = new CallingRegionStats();
		crs.regionsFile = args[0];
		crs.geneModelFile = args[1];
		try {
			crs.engage();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
