package au.edu.qimr.indel.pileup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.string.StringUtils;
import org.qcmg.common.util.Constants;
import org.qcmg.common.util.IndelUtils;
import org.qcmg.common.util.IndelUtils.SVTYPE;
import org.qcmg.common.vcf.VcfFormatFieldRecord;
import org.qcmg.common.vcf.VcfInfoFieldRecord;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.VcfUtils;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.vcf.VCFFileReader;

public class ReadGatkIndels extends ReadIndels{
	
	static final String GATK_INFO_DP = "DP";
	static final String GATK_INFO_AN = "AN";	
	static final String GATK_INFO_AC = "AC";
	static final String GATK_INFO_AF = "AF";
	static final String GATK_INFO_MLEAC = "MLEAC";
	static final String GATK_INFO_MLEAF = "MLEAF";	 
		
	static final String GATK_FORMAT_GT = VcfHeaderUtils.FORMAT_GENOTYPE;   //"GT";
	static final String GATK_FORMAT_GD = VcfHeaderUtils.FORMAT_GENOTYPE_DETAILS; //GD
	static final String GATK_FORMAT_AD = VcfHeaderUtils.FORMAT_ALLELIC_DEPTHS ;  //"AD";
	static final String GATK_FORMAT_DP = VcfHeaderUtils.FORMAT_READ_DEPTH;    //"DP";
	static final String GATK_FORMAT_PL = VcfHeaderUtils.FORMAT_POSSIBLE_GENOTYPES; //"PL";
	static final String GATK_FORMAT_GQ = VcfHeaderUtils.FORMAT_GENOTYPE_QUALITY; //"GQ";
	
	
	public ReadGatkIndels(QLogger logger) { super(logger); }
		
	/**
	 * 
	 * @param f: input vcf file, only keep the first sample column
	 * @param isTestFile: if true, it will be merged to control vcf and it's first sample column will be append to second column;
	 * Otherwise keep the control vcf and with only the first sample column. 
	 * @throws IOException
	 */
	public void readFile(File f, boolean isTest) throws IOException{
		int indelNew = 0;
		int indelDup = 0;
		int indelMultiAltNo = 0;
    	int inLines = 0;
    	int inMultiAltNo = 0;
    	int errVcf = 0; 
	
        try (VCFFileReader reader = new VCFFileReader(f)) {        	 
        	header = (!isTest)? reader.getHeader(): VcfHeaderUtils.mergeHeaders(header, reader.getHeader(), false);;	        	
         	//no chr in front of position
			for (final VcfRecord re : reader) {					
			   inLines ++;				 
	 	       addGenoDetails(re);	
	 	       String[] alleles = re.getAlt().split(",");
	 	       if(alleles.length > 1) inMultiAltNo ++; //multi allele input variants
	 	        
	 	       //keep the first sample column only
 				for(int i = 0; i < alleles.length; i++){						 					 
					SVTYPE type = IndelUtils.getVariantType(re.getRef(), alleles[i]);
					if( !type.equals(SVTYPE.DEL) && !type.equals(SVTYPE.INS)) continue; 
					
					//retrive this single allele indel or split alleles, 
					VcfRecord  vcf1 = ( alleles.length == 1)? re: decomposeVcf(re, alleles[i], i) ;						
					vcf1.setFilter(Constants.MISSING_DATA_STRING); //reset filter
					
					//only keep first sample column of second vcf
					VcfFormatFieldRecord fre = vcf1.getSampleFormatRecord(1);
					if(fre == null && (errVcf ++) < errRecordLimit){
						logger.warn("input vcf record missing format column:" + vcf1.toSimpleString());
						continue; 
					}else	
						vcf1.setFormatFields( fre.toStringList());	
					
					if(!isTest){
						if(positionRecordMap.containsKey(vcf1) && (indelDup ++) < errRecordLimit)					
	     					logger.warn("duplicate variants exsits and only keep last one:" + vcf1.toSimpleString());	     							 
						else
							positionRecordMap.put(vcf1, vcf1);  
					}else
						mergeIndel(positionRecordMap.get(vcf1), vcf1) ;
				}   
    		}
 	      }
                      
    	//counts from each input
		logger.info(indelNew + " indels are found from first input file: " + f.getCanonicalPath());
		logger.info(indelMultiAltNo + " indels are split from multi alleles in first input file.");
		logger.info(inLines  + " variants record exsits inside first input file.");
		logger.info(inMultiAltNo + " variants record with multi Alleles exsits inside first input file.");						
	}
	
	/**
	 * 
	 * @param re: input vcf record with only one sample column
	 * @param isTestFile: if true, shif the sample column to second column, otherwise to first
	 * @return true if this vcf is first time appear
	 */
	private void mergeIndel(  VcfRecord control, VcfRecord test){
		
		//a somatic vcf
		if(control == null){			 
			VcfUtils.addMissingDataToFormatFields(test, 1);			
			positionRecordMap.put(test, test);								
			return ; 
		}		
		
		//merge format, here control vcf have only one sample column			
		VcfUtils.addAdditionalSampleToFormatField(control, test.getSampleFormatRecord(1).toStringList());

		//merge INFO AD 		 
		VcfInfoFieldRecord infoControl = control.getInfoRecord();					
		for( String key : new String[]{GATK_INFO_AC, GATK_INFO_AN,GATK_INFO_MLEAC, GATK_INFO_DP } ){
			String value1 =  infoControl.getField(key);
			String value2 = test.getInfoRecord().getField(key);			
			if(value1 == null || value2 == null) continue; //missing sub-field			
			int join = Integer.getInteger(value1) + Integer.getInteger(value2);				
			infoControl.addField(key, join+"");			 
		}
		
	}	
	
	/**
	 * only for GATK vcf record, split multi alleles to mutli vcf record
	 * @param re: original vcf record
	 * @param alleles: an array of splited alleles, maximum number is 3
	 * @param altOrder: current allele order, maximum number is 3
	 * @param gd: genotype details, eg. A/T
	 * @return new vcf record conatin only current allele. only restore current allele's AC, AF, MLEAC, MLEAF value, only have first sample column
	 * GT:GQ:DP:PL:AD are updated as well, add GD value to sample column 
	 */
	private VcfRecord decomposeVcf(VcfRecord re, String alt, int altOrder){
		
		//set current allele
		VcfRecord re1 = new VcfRecord.Builder(re.getChrPosition(), re.getRef(), alt)
				.id(re.getId()).qualString(re.getQualString()).build();
		
		//reset info column 
		VcfInfoFieldRecord infoRe = re.getInfoRecord();			
		//replace with splited AC and AF for each allel
		for(String key : new String[]{GATK_INFO_AC, GATK_INFO_AF,GATK_INFO_MLEAC, GATK_INFO_MLEAF  }){
			String value =  infoRe.getField(key);
			if(value == null ) continue; //missing sub-field
			String[] acs = infoRe.getField(key).split(",");
			if(altOrder > acs.length)
				logger.warn(key + " value in INFO sub-field is invalid: " + re.toString());
			else
				infoRe.setField(key, acs[altOrder]) ;
		}
		re1.setInfo(infoRe.toString()); 
		
		List<String> format = re.getFormatFields();
		
		//format column will be null if original missing data  
		if(format == null && format.size() < 2) return re1; 
		
		//reset format here the list is 2 element: formart and sample column1			
		VcfFormatFieldRecord formatRe = new  VcfFormatFieldRecord(format.get(0), re.getFormatFields().get(1));
		String[] sFormat = new String[6];
		String[] sKeys = new String[]{ VcfHeaderUtils.FORMAT_GENOTYPE, VcfHeaderUtils.FORMAT_GENOTYPE_DETAILS , VcfHeaderUtils.FORMAT_ALLELIC_DEPTHS,
				VcfHeaderUtils.FORMAT_READ_DEPTH, VcfHeaderUtils.FORMAT_POSSIBLE_GENOTYPES, VcfHeaderUtils.FORMAT_GENOTYPE_QUALITY } ;
				 				
		for( int i = 0; i < 6; i ++){
			String ff = formatRe.getField( sKeys[i] );
			sFormat[i] = (StringUtils.isNullOrEmptyOrMissingData(ff)) ? Constants.MISSING_DATA_STRING : ff; 
		}
	 
		String[] ads = sFormat[2].split(",");
		String[] pls = sFormat[5].split(",");
						 			 	
		//only deal with 0/1 1,1 and 1/2
		if(altOrder == 0){     //alt1
			sFormat[0].replace("2",".").replace("3", ".");
			sFormat[2] = ( ads.length > 1)?  ads[0] + "," + ads[1] : ".,.";		//ref and alt1	
			sFormat[5] = ( pls.length > 2)?  pls[0] + "," + pls[1] + "," + pls[2] : ".,.,."; // 00, 01, 11					
		}else if(altOrder == 1){  //alt2
			sFormat[0].replace("2","1").replace("3", ".").replace("1", "."); 
			sFormat[2] =  ( ads.length > 2)? ads[0] + "," + ads[2] : ".,.";		//ref and alt2		
			sFormat[5] = (pls.length > 5)? pls[0] + "," + pls[3] + "," + pls[5] : ".,.,."; // 00, 01, 11								
		}else  if(altOrder == 2){ //alt3
			sFormat[0].replace("3","1").replace("2", ".").replace("1", ".");
			sFormat[2] =  (ads.length > 3)? ads[0] + "," + ads[3] : ".,.";		//ref and alt2		
			sFormat[5] = (pls.length > 9)? pls[0] + "," + pls[6] + "," + pls[9] : ".,.,."; // 00, 01, 11
																	
		}else{ //can't deal with alleles more than 3
			sFormat[0] = (sFormat[0].contains("/") ) ? "./." : ".|."; 
			sFormat[2] = ".,.";
			sFormat[5] =  ".,.,.";			
		}
		//GT:GD:AD:DP:GQ:PL  original DP for total coverage; original GQ for confidence . for confident unsure 
 		
		format.clear();
		format.add(String.join(":", sKeys));
		format.add(String.join(":", sFormat));		
		re1.setFormatFields(format);	
		
		return re1; 
	}
		
	/**
	 * 
	 * @param vcf: add Genotyp Details into sample column.
	 *  eg. ...A T,C.. GT 0/1 1/2, become ...A T,C.. G:GD 0/1:A/T 1/1:T/C
	*/
	private void addGenoDetails ( VcfRecord vcf ){	
		List<String> format = vcf.getFormatFields();
		if(format == null || format.size() < 2) return ; 
				
		//add GD to second field 
		List<String> gdFormat = new ArrayList<String>();		 
		for(int i = 1; i < format.size(); i ++){			
			VcfFormatFieldRecord re = new  VcfFormatFieldRecord(format.get(0), format.get(i));
			String gd = IndelUtils.getGenotypeDetails(re, vcf.getRef(), vcf.getAlt() );
			re.setField(1, VcfHeaderUtils.FORMAT_GENOTYPE_DETAILS, gd);			
			if(i == 1) gdFormat.add(re.getFormatColumnString());
			gdFormat.add(re.getSampleColumnString());
		}
				
		vcf.setFormatFields(gdFormat);
	}

	
//	public void readFirstFile(File first ) throws IOException{		
//	int indelNew = 0;
//	int indelDup = 0;
//	int indelMultiAltNo = 0;
//	int inLines = 0;
//	int inMultiAltNo = 0;
//
//    try (VCFFileReader reader = new VCFFileReader(first)) {        	 
//    	header = reader.getHeader();	
//     	//no chr in front of position
//
//		for (final VcfRecord re : reader) {					
//		   inLines ++;			
//		   //format data from control, set default as germline
// 	       re.setFilter(Constants.MISSING_DATA_STRING);
// 	       addGenoDetails(re);	 	        
// 	       String[] alleles = re.getAlt().split(",");
// 	       if(alleles.length > 1) inMultiAltNo ++; //multi allele input variants
// 	        
//				for(int i = 0; i < alleles.length; i++){						 					 
//				SVTYPE type = IndelUtils.getVariantType(re.getRef(), alleles[i]);
//				if( !type.equals(SVTYPE.DEL) && !type.equals(SVTYPE.INS)) continue; 
//				
//				//retrive this single allele indel or split alleles
//				VcfRecord  vcf1 = ( alleles.length == 1)? re: decomposeVcf(re, alleles[i], i) ;					
//					if(positionRecordMap.containsKey(vcf1) && (indelDup ++) < errRecordLimit){						
//						logger.warn("same variants already exsits, this one will be discard:\n" + positionRecordMap.get(vcf1).toString() );
//						continue; //no overwrite but just warning
//					}
//					positionRecordMap.put(vcf1, vcf1);  
//					indelNew ++;   	 	        	 
//			}   
//		}
//	      }
//           
//	//counts from each input
//	logger.info(indelNew + " indels are found from first input file: " + first.getCanonicalPath());
//	logger.info(indelMultiAltNo + " indels are split from multi alleles in first input file.");
//	logger.info(inLines  + " variants record exsits inside first input file.");
//	logger.info(inMultiAltNo + " variants record with multi Alleles exsits inside first input file.");							    		                
//}
	
}
