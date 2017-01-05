package au.edu.qimr.indel.pileup;

import java.io.File;
import java.io.IOException;
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
		
	static final String GATK_FORMAT_GT = "GT";
	static final String GATK_FORMAT_AD = "AD";
	static final String GATK_FORMAT_DP = "DP";
	static final String GATK_FORMAT_PL = "PL";
	static final String GATK_FORMAT_GQ = "GQ";
	
	public ReadGatkIndels(QLogger logger) {
		super(logger);
		 
	}
	
	public void readFirstFile(File first ) throws IOException{		
		int indelNew = 0;
		int indelDup = 0;
		int indelMultiAltNo = 0;
    	int inLines = 0;
    	int inMultiAltNo = 0;
	
        try (VCFFileReader reader = new VCFFileReader(first)) {
        	 
        	header = reader.getHeader();	
         	//no chr in front of position
			for (final VcfRecord re : reader) {					
				inLines ++;			
				//format data from control, set default as germline
	 	        re.setFilter(Constants.MISSING_DATA_STRING);
	 	        		
	 	        		
				//add genotype details GD to sample format
				String gd = IndelUtils.getGenotypeDetails( new VcfFormatFieldRecord(re.getFormatFields().get(0),re.getFormatFields().get(1)),
						re.getRef(), re.getAlt() );				
				
				
			//	resetGenotype(re);
    			String[] alleles = re.getAlt().split(",");
    			if(alleles.length > 1) inMultiAltNo ++; //multi allele input variants
    			
				 
				for(int i = 0; i < alleles.length; i++){				
					SVTYPE type = IndelUtils.getVariantType(re.getRef(), alleles[i]);
					if( !type.equals(SVTYPE.DEL) && !type.equals(SVTYPE.INS)) continue; 
					 
					VcfRecord vcf1 = getSplitVcf(re, alleles, i, gd);
//					VcfRecord vcf1 = re; 	 	        	   	
// 	        		if(alleles.length > 1){	
// 	        			vcf1 = VcfUtils.resetAllel(re, alleles[i]); //reset allele column
// 	        			indelMultiAltNo ++; //mutli allele indel
// 	        		}
 	        		
 					if(positionRecordMap.containsKey(vcf1) && (indelDup ++) < errRecordLimit){						
 						logger.warn("same variants already exsits, this one will be discard:\n" + positionRecordMap.get(vcf1).toString() );
 						continue; //no overwrite but just warning
 					}
 					positionRecordMap.put(vcf1, vcf1);  
 					indelNew ++;
   	 	        	 
				}  //for each allele		 
    		}
 	      }
               
    	//counts from each input
		logger.info(indelNew + " indels are found from first input file: " + first.getCanonicalPath());
		logger.info(indelMultiAltNo + " indels are split from multi alleles in first input file.");
		logger.info(inLines  + " variants record exsits inside first input file.");
		logger.info(inMultiAltNo + " variants record with multi Alleles exsits inside first input file.");							
    		                
	}	
	/**
	 * only for GATK vcf record, split multi alleles to mutli vcf record
	 * @param re: original vcf record
	 * @param alleles: an array of splited alleles, maximum number is 3
	 * @param altOrder: current allele order, maximum number is 3
	 * @param gd: genotype details, eg. A/T
	 * @return new vcf record conatin only current allele. only restore current allele's AC, AF, MLEAC, MLEAF value; 
	 * GT:GQ:DP:PL:AD are updated as well, add GD value to sample column 
	 */
	private VcfRecord getSplitVcf(VcfRecord re, String[] alleles, int altOrder, String gd){
		
		VcfRecord re1 =  re; 
		
		if(alleles.length > 1)  {		
			//set current allele
			re1  =  new VcfRecord.Builder(re.getChrPosition(), re.getRef(), alleles[altOrder])
					.id(re.getId()).qualString(re.getQualString()).filter(re.getFilter()).build();
			
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
			
		
			//reset format here the list is 2 element: formart and sample column1
			List<String> format = re.getFormatFields();
			//do nothing for null
			if(format != null && format.size() > 1){				
				VcfFormatFieldRecord formatRe = new  VcfFormatFieldRecord(format.get(0), re.getFormatFields().get(1));
				String gt = formatRe.getField(GATK_FORMAT_GT);
				
				//only deal with 0/1 1,1 and 1/2
				//first allele
				if(altOrder == 0){
					gt.replace("2",".").replace("3",".");
				}
				
				//re1.setFormatFields(re.getFormatFields());
/*
 * 	static final String GATK_FORMAT_GT = "GT";
	static final String GATK_FORMAT_AD = "AD";
	static final String GATK_FORMAT_DP = "DP";
	static final String GATK_FORMAT_PL = "PL";
	static final String GATK_FORMAT_GQ = "GQ";
			
 */
			}
		}
		
		
		return re1; 
	}
	

}
