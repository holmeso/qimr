package au.edu.qimr.vcftools.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Pair;
import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.string.StringUtils;
import org.qcmg.common.util.Constants;
import org.qcmg.common.util.SnpUtils;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.VcfUtils;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.common.vcf.header.VcfHeader.FormattedRecord;
import org.qcmg.common.vcf.header.VcfHeader.Record;
import org.qcmg.common.vcf.header.VcfHeaderUtils;

import au.edu.qimr.vcftools.Rule;

public class MergeUtils {
	
	private static final String GENOTYPE = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	public static final String FORMAT = "FORMAT";
	private static final QLogger logger = QLoggerFactory.getLogger(MergeUtils.class);

	public static List<String> mergeOtherHeaderRecords(List<Record> ...  loloRecs) {
		if (null == loloRecs || loloRecs.length == 0) {
			return Collections.emptyList();
		}
		
		AtomicInteger prefix = new AtomicInteger(1);
		List<String> mergedRecs = new ArrayList<>();
		Arrays.stream(loloRecs)
			.filter(list -> list != null && ! list.isEmpty())
			.forEach(list -> {
				mergedRecs.addAll(list.stream()
						.filter(r -> r != null && r.getData() != null && ! r.getData().equals(VcfHeaderUtils.BLANK_HEADER_LINE))
						.map(r -> r.getData().replaceAll(Constants.DOUBLE_HASH, Constants.DOUBLE_HASH + prefix.get() + Constants.COLON))
						.collect(Collectors.toList()));
				prefix.incrementAndGet();
			});
			
		return mergedRecs;
	}
	
	public static Map<Integer, Map<String, String>> getRules(Map<String, FormattedRecord> ...  loMaRecs) {
		return getHeaderAndRules(Arrays.asList(loMaRecs)).getRight();
	}
	
	public static Pair<VcfHeader, Rule> getMergedHeaderAndRules(VcfHeader ... headers) {
		
		if (canMergeBePerformed(headers)) {
			VcfHeader mergedHeader = new VcfHeader();
			List<Map<String, FormattedRecord>> infoHeaders = new ArrayList<>(headers.length);
			List<Map<String, FormattedRecord>> filterHeaders = new ArrayList<>(headers.length);
			List<Map<String, FormattedRecord>> formatHeaders = new ArrayList<>(headers.length);
			List<List<Record>> otherHeaders = new ArrayList<>(headers.length);
			for (VcfHeader h : headers) {
				infoHeaders.add(h.getInfoRecords());
				filterHeaders.add(h.getFilterRecords());
				formatHeaders.add(h.getFormatRecords());
				otherHeaders.add(h.getNonStandardRecords());
			}
			
			Pair<List<FormattedRecord>, Map<Integer, Map<String, String>>> infoPair = getHeaderAndRules(infoHeaders);
			for (FormattedRecord fr : infoPair.getLeft()) {
				mergedHeader.addInfo(fr);
			}
			Pair<List<FormattedRecord>, Map<Integer, Map<String, String>>> filterPair = getHeaderAndRules( filterHeaders);
			for (FormattedRecord fr : filterPair.getLeft()) {
				mergedHeader.addFilter(fr);
			}
			Pair<List<FormattedRecord>, Map<Integer, Map<String, String>>> formatPair = getHeaderAndRules(formatHeaders);
			for (FormattedRecord fr : formatPair.getLeft()) {
				mergedHeader.addFormat(fr);
			}
			List<String> mergedOtherRecords = mergeOtherHeaderRecords(otherHeaders.toArray(new List[]{}));
			for (String s : mergedOtherRecords) {
				mergedHeader.parseHeaderLine(s);
			}
			
			/*
			 * add IN= to 
			 */
			logger.info("checking that no existing info lines start with "+Constants.VCF_MERGE_INFO+"= ");
			if (infoPair.getLeft().contains(Constants.VCF_MERGE_INFO)) {
				logger.warn("Can't use "+Constants.VCF_MERGE_INFO+"= to mark records as having come from a particular input file - "+Constants.VCF_MERGE_INFO+"= is already in use!");
			}
			
			mergedHeader.addInfoLine(Constants.VCF_MERGE_INFO, ".","Integer", VcfHeaderUtils.DESCRITPION_MERGE_IN);
			mergedHeader.addFormatLine(VcfHeaderUtils.FORMAT_INFO, ".","String", VcfHeaderUtils.FORMAT_INFO_DESCRIPTION);
			//		 "Indicates which INput file this vcf record came from. Multiple values are allowed which indicate that the record has been merged from more than 1 input file");
			
			/*
			 * make sure SOMATIC has been added, and add the _n entry
			 */
			if ( ! infoPair.getLeft().contains(VcfHeaderUtils.INFO_SOMATIC)) {
				mergedHeader.addInfoLine(VcfHeaderUtils.INFO_SOMATIC, "0", "Flag", VcfHeaderUtils.INFO_SOMATIC_DESC);
			}
//			for (int i = 1 ; i <= headers.length ; i++) {
//				mergedHeader.addInfoLine(VcfHeaderUtils.INFO_SOMATIC + "_" + i, "0", "Flag", "Indicates that the nth input file considered this record to be somatic. Multiple values are allowed which indicate that more than 1 input file consider this record to be somatic");
//			}
			/*
			 * Build up new CHROM line with the samples updated with the input number appended to them
			 */
			StringBuilder sb = new StringBuilder();
			String [] firstChrLine = headers[0].getChrom().toString().split(Constants.TAB_STRING);
			/*
			 * add the first 9 columns to sb
			 */
			for (int i = 0 ; i < 9 ; i++) {
				StringUtils.updateStringBuilder(sb, firstChrLine[i], Constants.TAB);
			}
			for (int i = 0 ; i < headers.length ; i++) {
				String [] array = headers[i].getChrom().toString().split(Constants.TAB_STRING);
				/*
				 * add every column after FORMAT to sb, appending the numeric identifier
				 */
				boolean go = false;
				for (String s : array) {
					if (go) {
						StringUtils.updateStringBuilder(sb, s + "_" + (i+1), Constants.TAB);
					}
					if (FORMAT.equals(s)) {
						go = true;
					}
				}
			}
			mergedHeader.parseHeaderLine(sb.toString());
//			mergedHeader.parseHeaderLine(headers[0].getChrom().toString());
			
			Rule r = new Rule(headers.length);
			for (int i = 0 ; i < headers.length ; i++) {
				r.setRules(i, filterPair.getRight().get(i), infoPair.getRight().get(i), formatPair.getRight().get(i));
			}
			
			return Pair.of(mergedHeader,  r);
			
		} else {
			logger.warn("Unable to perform merge - please check that the vcf headers contain the same samples in the same order");
		}
		return null;
		
		
	}
	
	public static List<FormattedRecord> mergeHeaderRecords(Map<String, FormattedRecord> ...  loMaRecs) {
		return getHeaderAndRules(Arrays.asList(loMaRecs)).getLeft();
	}
	
	public static Pair<List<FormattedRecord>, Map<Integer, Map<String, String>>> getHeaderAndRules(List<Map<String, FormattedRecord> >  loMaRecs) {
		if (null == loMaRecs || loMaRecs.isEmpty()) {
			return Pair.of(Collections.emptyList(), Collections.emptyMap());
		}
		
		Map<Integer,Map<String, String>> replacementIds = new HashMap<>(4);
		
		Map<String, FormattedRecord> mergedRecsMap = loMaRecs.get(0);
		for (int i = 1 ; i < loMaRecs.size() ; i++) {
			Map<String, FormattedRecord> map = loMaRecs.get(i);
			for (Entry<String, FormattedRecord> entry : map.entrySet()) {
				FormattedRecord mergedRec = mergedRecsMap.get(entry.getKey());
				if (null == mergedRec) {
					// add
					mergedRecsMap.put(entry.getKey(), entry.getValue());
				} else {
					if (mergedRec.equals(entry.getValue())) {
						logger.info("Found identical header entry: " + mergedRec.toString());
					} else if (mergedRec.getData().substring(mergedRec.getData().indexOf("Description="))
							.equals(entry.getValue().getData().substring(entry.getValue().getData().indexOf("Description=")))) {
						/*
						 * Just match on description - ignore type and number for now (ever??)...
						 */
						logger.info("Found identical header entry (apart from type, and number): " + mergedRec.toString());
					} else if (mergedRec.getData().equals(GENOTYPE) || entry.getValue().getData().equals(GENOTYPE)) {
						/*
						 * Just match on description - ignore type and number for now (ever??)...
						 */
						logger.info("Found identical header entry (apart from type, and number): " + mergedRec.toString());
					} else {
						String newId = entry.getKey() + i;
						String existingRec = entry.getValue().getData();
						
						FormattedRecord updatedRec = existingRec.startsWith(VcfHeaderUtils.HEADER_LINE_INFO) ? new VcfHeader.InfoRecord(entry.getValue().getData().replace(entry.getKey(), newId))
						: existingRec.startsWith(VcfHeaderUtils.HEADER_LINE_FILTER) ? new VcfHeader.FilterRecord(entry.getValue().getData().replace(entry.getKey(), newId))
						: existingRec.startsWith(VcfHeaderUtils.HEADER_LINE_FORMAT) ? new VcfHeader.FormatRecord(entry.getValue().getData().replace(entry.getKey(), newId))
						: null;
						
						logger.info("bumping id from " + entry.getKey() + " to " + newId + " and adding to map. orig: " + mergedRec + ", updated: " + updatedRec);
						mergedRecsMap.put(newId, updatedRec);
						Map<String, String> replacement = replacementIds.get(i);
						if (null == replacement) {
							replacement = new HashMap<>();
							replacementIds.put(i, replacement);
						}
						replacement.put(entry.getKey(), newId);
					}
				}
			}
		}
		
		List<FormattedRecord> mergedRecs = new ArrayList<>(mergedRecsMap.values());
		return  Pair.of(mergedRecs, replacementIds);
	}
	
	public static Optional<String> getCombinedAlt(VcfRecord ... recs) {
		return Optional.ofNullable ( 
				Arrays.stream(
					Arrays.stream(recs)
					.map(VcfRecord::getAlt)
					.collect(Collectors.joining(Constants.COMMA_STRING))
				.split(Constants.COMMA_STRING))
				.distinct()
				.collect(Collectors.joining(Constants.COMMA_STRING)));
	}
	
	public static String getGT(String combinedAlts, String myAlts, String currentGT) {
		if ("0/0".equals(currentGT)) {
			return currentGT;
		}
		if (combinedAlts.equals(myAlts)) {
			return currentGT;
		}
		int index = currentGT.indexOf(Constants.SLASH_STRING);
		int gt1 = Integer.parseInt(currentGT.substring(0, index));
		int gt2 = Integer.parseInt(currentGT.substring(index + 1));
		
		String [] myAltsArray = myAlts.split(Constants.COMMA_STRING);
		String [] combinedAltsArray = combinedAlts.split(Constants.COMMA_STRING);
		
		int newGT1 = getNewPosition(combinedAltsArray, myAltsArray, gt1) ;
		int newGT2 = getNewPosition(combinedAltsArray, myAltsArray, gt2) ;
		
		return newGT1 + Constants.SLASH_STRING + newGT2;
	}
	
	public static int getNewPosition(String[] newAlts, String[] oldAlts, int oldPos) {
		if (oldPos == 0) return oldPos;
		
		String alt = oldAlts[oldPos - 1];
		int newPos = -1;
		
		for (int i = 1 ; i <= newAlts.length ; i++) {
			if (alt.equals(newAlts[i - 1])) {
				newPos = i;
				break;
			}
		}
		return newPos;
	}
	
	/**
	 * Moves any data that we don't like in the filter and info fields into specific fields in the format column.
	 * Modifies the VcfRecord supplied as an argument in place, and so is not referentially transparent
	 */
	public static void moveDataToFormat(VcfRecord r, String combinedAlt) {
		/*
		 * for info field, we are just looking for SOMATIC
		 * If this is present, move to the FORMAT-INF field for all samples (assuming that the format
		 */
		List<String> ff = r.getFormatFields();
		Map<String, String[]> ffMap = VcfUtils.getFormatFieldsAsMap(ff);
		int sampleCount = ff.size() -1;
		if (sampleCount > 0) {
//			if (VcfUtils.isRecordSomatic(r)) {
//				/*
//				 * Check to see if we have an INF entry in the map
//				 */
//				String [] infArr = ffMap.computeIfAbsent(VcfHeaderUtils.FORMAT_INFO, k -> new String[sampleCount]);
//				for (int i = 0 ; i < infArr.length ; i++) {
//					String s = infArr[i];
//					if (StringUtils.isNullOrEmptyOrMissingData(s)) {
//						infArr[i] = VcfHeaderUtils.INFO_SOMATIC;
//					} else {
//						infArr[i] += Constants.SEMI_COLON + VcfHeaderUtils.INFO_SOMATIC;
//					}
//				}
//			}
			
			/*
			 * filter field
			 */
			String filter = r.getFilter();
			if ( ! StringUtils.isNullOrEmptyOrMissingData(filter)) {
				String [] infArr = ffMap.computeIfAbsent(VcfHeaderUtils.FORMAT_FILTER, k -> new String[sampleCount]);
				for (int i = 0 ; i < infArr.length ; i++) {
					String s = infArr[i];
					if (StringUtils.isNullOrEmptyOrMissingData(s)) {
						infArr[i] = filter;
					} else {
						infArr[i] += Constants.SEMI_COLON + filter;
					}
				}
			}
		}
		
		/*
		 * if the combinedAlts is different from the alt for this record, we may need to update the GT field
		 */
		String aa = r.getAlt();
		if ( ! aa.equals(combinedAlt)) {
			
			String[] existingGTs = ffMap.get(VcfHeaderUtils.FORMAT_GENOTYPE);
			/*
			 * we could have an updated gts
			 */
			for (int z = 0 ; z < existingGTs.length ; z++) {
				String gt = existingGTs[z];
				if ( ! gt.equals(Constants.MISSING_DATA_STRING)) {
					String newGT = getGT(combinedAlt, aa, gt);
					if ( ! newGT.equals(gt)) {
						existingGTs[z] = newGT;
					}
				}
			}
			ffMap.put(VcfHeaderUtils.FORMAT_GENOTYPE, existingGTs);
		}
		
		/*
		 * update record
		 */
		if ( ! ffMap.isEmpty()) {
			r.setFormatFields(VcfUtils.convertFFMapToList(ffMap));
		}
	}
	
	/**
	 * returns true if the supplied string only contains the missing data symbol '.' and the colon delimiter ':'
	 * eg. .:.:.:.:.:.:.:. would return true, .:.:.:SOMATIC:.:.:. would return false.
	 * 
	 * @param f
	 * @return
	 */
	public static boolean isFormatSampleEmpty(String f) {
		if (StringUtils.isNullOrEmptyOrMissingData(f)) {
			return true;
		}
//		f.chars().forEach(System.out::println);
		return f.chars().allMatch(i -> i == 46 || i == 58);
	}
	
	/**
	 * Returns a new VcfRecord that takes the positional information from the suppliued VCfRecord, and uses the supplied alt for the alt...
	 * @param r
	 * @param combinedAlt
	 * @return
	 */
	public static VcfRecord getBaseVcfRecordDetails(VcfRecord r, String combinedAlt) {
		return new VcfRecord.Builder(r.getChrPosition(), r.getRef()).id(r.getId()).allele(combinedAlt).build();
	}
	
	public static VcfRecord mergeRecords(Map<Integer, Map<String, String>> rules, VcfRecord caller1, VcfRecord caller2) {
		VcfRecord mr = null;
		
		
		if (null == caller1 && null == caller2) {
			throw new IllegalArgumentException("MergeUtils.mergeRecords called will no records!!!");
		} else if (null == caller1) {
			Map<String, String> thisRecordsRules = null != rules ? rules.get(1) : null;
			/*
			 * just got caller 2 data
			 */
			mr = getBaseVcfRecordDetails(caller2, caller2.getAlt());
			mr.setInfo(caller2.getInfo());
			moveDataToFormat(caller2, mr.getAlt());
			
			/*
			 * add format fields from caller2 to mr
			 */
			List<String> ff = caller2.getFormatFields();
			mr.setFormatFields(ff);
			
			/*
			 * add format columns for missing caller 
			 */
			int noOfSamples = ff.size() -1;
			for (int i = 0 ; i < noOfSamples ; i++) {
				VcfUtils.addMissingDataToFormatFields(mr, 1);
			}
			
		} else if (null == caller2) {
			Map<String, String> thisRecordsRules = null != rules ? rules.get(0) : null;
			/*
			 * just got caller 1 data
			 */
			mr = getBaseVcfRecordDetails(caller1, caller1.getAlt());
			moveDataToFormat(caller2, mr.getAlt());
			/*
			 * add format fields from caller1 to mr
			 */
			List<String> ff = caller1.getFormatFields();
			mr.setFormatFields(ff);
			
			/*
			 * add format columns for missing caller 
			 */
			int noOfSamples = ff.size() -1;
			for (int i = 0 ; i < noOfSamples ; i++) {
				VcfUtils.addMissingDataToFormatFields(mr, noOfSamples);
			}
			
			
		} else {
			/*
			 * got data from both callers
			 */
			
			Optional<String> combinedAlts = getCombinedAlt(caller1, caller2);
			mr = getBaseVcfRecordDetails(caller1, combinedAlts.get());
			
			Map<String, String> caller1Rules = null != rules ? rules.get(0) : null;
			if (null != caller1Rules && ! caller1Rules.isEmpty()) {
				
				/*
				 * INFO
				 */
				for (String s : caller1.getInfo().split(Constants.SEMI_COLON_STRING)) {
					int equalsIndex = s.indexOf(Constants.EQ);
					String key = equalsIndex > -1 ? s.substring(0, equalsIndex) : s;
					String replacementKey = caller1Rules.get(key);
					if (null == replacementKey) {
						mr.getInfoRecord().appendField(key, (equalsIndex > -1) ? s.substring(equalsIndex+1) : s);
					} else {
						mr.getInfoRecord().addField(replacementKey, (equalsIndex > -1) ? s.substring(equalsIndex+1) : s);
					}
				}
				
			} else {
				mr.appendInfo(caller1.getInfo(), false);
			}
			Map<String, String> caller2Rules = null != rules ? rules.get(1) : null;
			if (null != caller2Rules && ! caller2Rules.isEmpty()) {
				
				/*
				 * INFO
				 */
				for (String s : caller2.getInfo().split(Constants.SEMI_COLON_STRING)) {
					int equalsIndex = s.indexOf(Constants.EQ);
					String key = equalsIndex > -1 ? s.substring(0, equalsIndex) : s;
					String replacementKey = caller2Rules.get(key);
					if (null == replacementKey) {
						mr.getInfoRecord().appendField(key, (equalsIndex > -1) ? s.substring(equalsIndex+1) : s);
					} else {
						mr.getInfoRecord().addField(replacementKey, (equalsIndex > -1) ? s.substring(equalsIndex+1) : s);
					}
				}
				
			} else {
				mr.appendInfo(caller2.getInfo(), false);
			}
			
			moveDataToFormat(caller1, mr.getAlt());
			moveDataToFormat(caller2, mr.getAlt());
			
			/*
			 * add format fields from caller1 to mr
			 */
			
			Map<String, String[]> caller1FFMap = VcfUtils.getFormatFieldsAsMap(caller1.getFormatFields());
			Map<String, String[]> caller2FFMap = VcfUtils.getFormatFieldsAsMap(caller2.getFormatFields());
			
			if (null != caller1Rules && ! caller1Rules.isEmpty()) {
				/*
				 * update any keys depending on rules
				 */
				for (Entry<String, String> entry : caller1Rules.entrySet()) {
					if (caller1FFMap.containsKey(entry.getKey())) {
						String [] values = caller1FFMap.remove(entry.getKey());
						caller1FFMap.put(entry.getValue(), values);
					}
				}
			}
			if (null != caller2Rules &&  ! caller2Rules.isEmpty()) {
				/*
				 * update any keys depending on rules
				 */
				for (Entry<String, String> entry : caller2Rules.entrySet()) {
					if (caller2FFMap.containsKey(entry.getKey())) {
						String [] values = caller2FFMap.remove(entry.getKey());
						caller2FFMap.put(entry.getValue(), values);
					}
				}
			}
			
			
			List<String> ff1 = VcfUtils.convertFFMapToList(caller1FFMap);
			if (ff1.size() > 1) {
				mr.setFormatFields(ff1);
			}
			List<String> ff2 =VcfUtils.convertFFMapToList(caller2FFMap);
			if (ff2.size() > 1) {
				VcfUtils.addAdditionalSamplesToFormatField(mr, ff2);
			}
			
			/*
			 * remove somatic from info field should it be there
			 */
			mr.getInfoRecord().removeField(SnpUtils.SOMATIC);
			
		}
		return mr;
	}
	
	/**
	 * Merge vcf records that have same position, ref, alt and number of samples into a single record
	 */
	public static VcfRecord mergeRecords(Map<Integer, Map<String, String>> rules, VcfRecord ... records) {
		if (null == records || records.length == 0) {
			throw new IllegalArgumentException("MergeUtils.mergeRecords called will no records!!!");
		}
		
		/*
		 * build up alt string
		 */
		
		Optional<String> combinedAlt = getCombinedAlt(records);
//		assert combinedAlt.isPresent() : "Null returned from getCombinedAlt";
		
		/*
		 * Get common values from 1st record
		 */ 
//		VcfRecord mergedRecord =  VcfUtils.createVcfRecord(records[0].getChrPosition(),   records[0].getId(), records[0].getRef(), records[0].getAlt());
		VcfRecord mergedRecord =  new VcfRecord.Builder(records[0].getChrPosition(), records[0].getRef())
									.id(records[0].getId()).allele(combinedAlt.get()).build();
				
		
 
		
		/*
		 * Update id, info, filter, and format fields
		 */
		
		for (int i = 0 ; i < records.length ; i++) {
			VcfRecord r = records[i];
			Map<String, String> thisRecordsRules = null != rules ? rules.get(i) : null;
			String suffix = "_" + (i+1);
			
			mergedRecord.appendId(r.getId());
			List<String> rFF =  r.getFormatFields() ;
			
			/*
			 * if the combinedAlts is different from the alt for this record, we may need to update the GT field
			 */
			String aa = r.getAlt();
			if ( ! aa.equals(combinedAlt.get())) {
				
				Map<String, String[]> ffMap = VcfUtils.getFormatFieldsAsMap(rFF);
				String[] existingGTs = ffMap.get(VcfHeaderUtils.FORMAT_GENOTYPE);
				/*
				 * we could have an updated gts
				 */
				boolean needToUpdate = false;
				for (int z = 0 ; z < existingGTs.length ; z++) {
					String gt = existingGTs[z];
					if ( ! gt.equals(Constants.MISSING_DATA_STRING)) {
						String newGT = getGT(combinedAlt.get(), aa, gt);
						if ( ! newGT.equals(gt)) {
							existingGTs[z] = newGT;
							needToUpdate = true;
						}
					}
				}
				
				
				/*
				 * update map
				 */
				if (needToUpdate) {
					ffMap.put(VcfHeaderUtils.FORMAT_GENOTYPE, existingGTs);
					rFF = VcfUtils.convertFFMapToList(ffMap);
					r.setFormatFields(rFF);
				}
				
			}
			/*
			 * remove SOMATIC from info field, and add to format INF subfield
			 */
			List<String> formatInfo = new ArrayList<>(rFF.size());
			formatInfo.add(VcfHeaderUtils.FORMAT_INFO);
			boolean isSomatic = VcfUtils.isRecordSomatic(r);
			
			for (int j = 1 ; j < rFF.size(); j++) {
				formatInfo.add(isSomatic && ! rFF.get(j).startsWith(Constants.MISSING_DATA_STRING) ? SnpUtils.SOMATIC : Constants.MISSING_DATA_STRING);
			}
//			formatInfo.add(isSomatic && ! rFF.get(1).startsWith(Constants.MISSING_DATA_STRING) ? SnpUtils.SOMATIC : Constants.MISSING_DATA_STRING);
//			formatInfo.add(isSomatic && ! rFF.get(2).startsWith(Constants.MISSING_DATA_STRING)? SnpUtils.SOMATIC : Constants.MISSING_DATA_STRING);
			VcfUtils.addFormatFieldsToVcf(r,formatInfo);
			if (isSomatic) {
				r.getInfoRecord().removeField(SnpUtils.SOMATIC);
			}
			
			
			if (null != thisRecordsRules && ! thisRecordsRules.isEmpty()) {
				
				/*
				 * INFO
				 */
				for (String s : r.getInfo().split(Constants.SEMI_COLON_STRING)) {
					int equalsIndex = s.indexOf(Constants.EQ);
					String key = equalsIndex > -1 ? s.substring(0, equalsIndex) : s;
					String replacementKey = thisRecordsRules.get(key);
					if (null == replacementKey) {
						mergedRecord.getInfoRecord().appendField(key, (equalsIndex > -1) ? s.substring(equalsIndex+1) : s);
					} else {
						mergedRecord.getInfoRecord().addField(replacementKey, (equalsIndex > -1) ? s.substring(equalsIndex+1) : s);
					}
				}
				
			} else {
				mergedRecord.appendInfo(r.getInfo(), false);
			}
			
			/*
			 * FORMAT
			 */
			if (null != rFF &&  ! rFF.isEmpty()) {
				if (null == thisRecordsRules) {
					if (i == 0) {
						VcfUtils.addFormatFieldsToVcf(mergedRecord, r.getFormatFields());
					} else {
						VcfUtils.addAdditionalSamplesToFormatField(mergedRecord, r.getFormatFields());
					}
//					VcfUtils.addFormatFieldsToVcf(mergedRecord, r.getFormatFields(), true);
				} else {
					/*
					 * create new header string, substituting in the new values should any be present in the rules map
					 */
					String [] formatHeadersArray = rFF.get(0).split(Constants.COLON_STRING);
					StringBuilder newHeader = new StringBuilder();
					for (String s : formatHeadersArray) {
						if (newHeader.length() > 0) {
							newHeader.append(Constants.COLON);
						}
						String replacementHeaderAttribute = thisRecordsRules.get(s);
						newHeader.append(replacementHeaderAttribute == null ? s : replacementHeaderAttribute);
					}
					if ( ! newHeader.toString().equals(rFF.get(0))) {
						rFF.set(0, newHeader.toString());
					}
					if (i == 0) {
						VcfUtils.addFormatFieldsToVcf(mergedRecord,rFF);
					} else {
						VcfUtils.addAdditionalSamplesToFormatField(mergedRecord, rFF);
					}
//					VcfUtils.addFormatFieldsToVcf(mergedRecord, rFF, true);
					
				}
			}
			
			
			/*
			 * FILTER
			 */
			if ( ! StringUtils.isNullOrEmptyOrMissingData(r.getFilter())) {
				for (String s : r.getFilter().split(Constants.SEMI_COLON_STRING)) {
					String replacementKey = (null == thisRecordsRules) ? null : thisRecordsRules.get(s);
					if (null == replacementKey) {
						mergedRecord.addFilter(s + suffix);
					} else {
						mergedRecord.addFilter(replacementKey + suffix);
					}
				}
			}
		}
		return mergedRecord;
	}
	
	/*
	 * Need to be same ChrPosition, ref and alt
	 */
//	public static boolean areRecordsEligibleForMerge(VcfRecord ... records) {
//		if (null == records || records.length == 0) {
//			return false;
//		}
//	
//		VcfRecord r1 = records[0];
//		return Arrays.stream(records)
//			.allMatch(r -> r.equals(r1));
//		
//	}
	
	/**
	 * Same sample merge test
	 * Need headers to contain the same samples in the same order
	 * Samples can't be null
	 * 
	 * Returns true if all headers contain the same samples in the same order, false otherwise
	 * 
	 * @param headers
	 * @return
	 */
	public static boolean canMergeBePerformed(VcfHeader ... headers) {
		
		if (null == headers || headers.length == 0) {
			return false;
		}
		String [] firstIds = headers[0].getSampleId();
		if (null == firstIds || firstIds.length == 0) {
			return false;
		}
		
		/*
		 * Get sample ids for each header and check that they are the same for each (number and order)
		 */
//		boolean doSampleIdsMatch = 
		return Arrays.stream(headers)
			.map(header -> header.getSampleId())
			.allMatch(array -> {
				return Arrays.deepEquals(array, firstIds);
			});
	}

}
