package au.edu.qimr.indel.pileup;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.model.ChrRangePosition;
import org.qcmg.common.string.StringUtils;
import org.qcmg.common.util.Constants;
import org.qcmg.common.util.IndelUtils;
import org.qcmg.common.util.IndelUtils.SVTYPE;
import org.qcmg.common.vcf.VcfFormatFieldRecord;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.VcfUtils;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.vcf.VCFFileReader;


public class ReadIndels {
	QLogger logger; 
	private VcfHeader header; 
	
	private final int errRecordLimit = 100;
	private int overwriteNo = 0;
	private int errGTNo = 0;
 
	//here key will be uniq for indel: chr, start, end, allel 
//	private final  Map<Integer, VcfRecord> positionRecordMap = new  ConcurrentHashMap<Integer, VcfRecord>();
	private final  Map<VcfRecord, VcfRecord> positionRecordMap = new  ConcurrentHashMap<VcfRecord, VcfRecord>();
	
//	private final HashSet<VcfRecord> positionRecordMap = new  HashSet<VcfRecord>();

	public ReadIndels( QLogger logger){
		this.logger = logger; 		
	}
	
	/**
	 * merge first sample column of input to existing variants which is stored on hash map
	 * @param f: input of vcf file
	 * @throws IOException
	 */
	public void appendIndels(File f) throws IOException{
				
		//only keep the first sample column, and put "." to second column
		for(VcfRecord vcf : positionRecordMap.values()){
			List<String> format = vcf.getFormatFields();
			if(format != null){
				while(format.size() > 2)
					format.remove(format.size()-1);   
				vcf.setFormatFields(format);			 
				VcfUtils.addMissingDataToFormatFields(vcf, 2);
			}			
		}	
		logger.info("only keep first sample column of tumour input vcf." );
		
    	//merge variants
       	int inLines = 0;
    	int inVariants = 0;   	
    	int indelnewCount = 0;
    	int overlapCount= 0;
    	int multiAltNo = 0;
	
        try (VCFFileReader reader = new VCFFileReader(f)) {
        	header = VcfHeaderUtils.mergeHeaders(header, reader.getHeader(), false);
			for (final VcfRecord re : reader) {	
				inLines ++;
				resetGenotype(re);
    			String[] StrAlts = re.getAlt().split(",");    
    			if(StrAlts.length > 1) multiAltNo ++;
				for(String alt : StrAlts){
					inVariants ++;
					SVTYPE type = IndelUtils.getVariantType(re.getRef(), alt);
	 	        	if(type.equals(SVTYPE.DEL) ||type.equals(SVTYPE.INS) ){	 
	 	        		VcfRecord vcf1 = re; 
	 	        		//reset allele column
	 	        		if(StrAlts.length > 1){
	 	        			vcf1 = VcfUtils.resetAllel(re, alt);
//	 	        					new VcfRecord(re.toString().trim().split("\t"));
//	 	        			vcf1.setAlt(alt);
	 	        		}
     					if(!mergeIndel(vcf1))
     						indelnewCount ++;
     					else overlapCount ++; 				 
	 	        	}
				}  	    			
			}
		}
		logger.info(String.format("Find %d indels from %d variants within file: %s", inVariants, inLines, f));
		logger.info("the number of reads contains multi Aleles from second file is " + multiAltNo);
		logger.info(indelnewCount + " indel variants are only appeared in second vcf! ");
		logger.info(overlapCount + " indel variants are appeared in both input vcf! ");			
		logger.info(positionRecordMap.size() + " indel variants position are selected into output! ");        
	}	
	
	
	/**
	 * Add this vcf record if not exists on RAM, move  the original first sample column to second column and mark missing data '.' on the first column; 
	 * Or merge this vcf record into existed variants: replace the second sample column of exist variants with first sample column of new one.
	 * @param secVcf: a vcf record
	 * @return true if same variants exist and merge them; otherwise return false by adding this new variants 
	 */
	public boolean mergeIndel(  VcfRecord secVcf){
	    //get vcf with same ChrPointPosition, ref and alt
		VcfRecord existingvcf = positionRecordMap.get(secVcf);
				
		//only keep first sample column of second vcf
		List<String> secformat = secVcf.getFormatFields();
		while(secformat != null && secformat.size() > 2 )
			secformat.remove(secformat.size()-1); 

		//copy secVcf with sample column "FORMAT <missing data> oriSample" only 
		if(existingvcf == null){
			existingvcf = new VcfRecord.Builder(secVcf.getChrPosition(), secVcf.getRef())
				.id(secVcf.getId()).allele(secVcf.getAlt()).filter(secVcf.getFilter()).build();
			existingvcf.setInfo(secVcf.getInfo());	
			
			//insert . to first sample column and then shift original sample to second 
			if(secformat != null){
				existingvcf.setFormatFields(secformat);				
				VcfUtils.addMissingDataToFormatFields(existingvcf, 1);
			}
			positionRecordMap.put(existingvcf, existingvcf);		
						
			return false; 
		}else{
			//only keep first sample column of exsiting vcf
			List<String> format1 = existingvcf.getFormatFields();
			
			if(format1 != null){
				while( format1.size() > 2 )
					format1.remove(format1.size()-1); 			
				existingvcf.setFormatFields(format1);
				//merge exiting and second vcf format, the exsitingvcf already inside map	
				VcfUtils.addAdditionalSampleToFormatField(existingvcf,  secformat) ;				
			}
			return true; 			
		}		 
	}	
	
	/**
	 * load variants to hash map
	 * @param f: vcf input file
	 * @param sampleCode: replace sample code inside the input vcf file 
	 * @throws IOException
	 */
	public void LoadIndels(File f) throws IOException{
    	int inVariants = 0;
    	int inLines = 0;
    	int multiAltNo = 0;
	
        try (VCFFileReader reader = new VCFFileReader(f)) {
        	if(header == null)
        		header = reader.getHeader();	
        	else
        		header = VcfHeaderUtils.mergeHeaders(header, reader.getHeader(), false);
        	//no chr in front of position
			for (final VcfRecord re : reader) {	
				inLines ++;
				resetGenotype(re);
    			String[] alleles = re.getAlt().split(",");
    			if(alleles.length > 1) multiAltNo ++;
    			
				for(String alt : alleles){
					inVariants ++;
					SVTYPE type = IndelUtils.getVariantType(re.getRef(), alt);
	 	        	if(type.equals(SVTYPE.DEL) ||type.equals(SVTYPE.INS) ){
	 	        		VcfRecord vcf1 = VcfUtils.resetAllel(re, alt);
     					if(positionRecordMap.containsKey(vcf1) && (overwriteNo ++) < errRecordLimit){ 						
     						logger.warn("overwriting same variants:\n" + positionRecordMap.get(vcf1).toString() );
    					}
     					positionRecordMap.put(vcf1, vcf1);    					
	 	        	}
				}  		 
    		}
 	      }
		  logger.info(String.format("Find %d indels from %d variants (%d records lines) within file: %s",
						positionRecordMap.size(), inVariants, inLines, f.getAbsoluteFile()));
		  if(overwriteNo >= errRecordLimit)
			  logger.warn("there are big number (" + overwriteNo + ") of same variant but maybe with different annotation are overwrited/removed");		
		  if(errGTNo >= errRecordLimit)
			  logger.warn("there are big number (" + overwriteNo + ") of input vcf contains invalid GT field on Format column");				  
		  
		  logger.info(String.format("There are %d input record contains multi alleles within file: %s", multiAltNo, f.getAbsoluteFile()));
	}
	
	/**
	 * change the input vcf by putting '.' on "GT" field if there are multi Alleles existis;
	 * do nothing if not multi Alleles or not GT field on format column
	 * @param vcf: input vcf record
	 */
	public void resetGenotype(VcfRecord vcf){
//		if(!vcf.getAlt().contains(","))
//			return;
		
		if(vcf.getFormatFields() == null)
			return; 
		
		//add GD to second field 
		List<String> format = vcf.getFormatFields();		
		String[] details = new String[format.size()];
		
		//store the merged field GT:GD
		details[0] = VcfHeaderUtils.FORMAT_GENOTYPE + ":" + VcfHeaderUtils.FORMAT_GENOTYPE_DETAILS;
		
		for(int i = 1; i < format.size(); i ++){
			VcfFormatFieldRecord re = new  VcfFormatFieldRecord(format.get(0), format.get(i));
			//if "GT" field is not exist, do nothing
			String sgt = re.getField(VcfHeaderUtils.FORMAT_GENOTYPE) ;
			details[i] = (vcf.getAlt().contains(Constants.COMMA_STRING))? Constants.MISSING_DATA_STRING : sgt; 
			details[i] += Constants.COLON_STRING;
			
			
			if(StringUtils.isNullOrEmpty(sgt) || sgt.equals(Constants.MISSING_DATA_STRING)){
				details[i] += Constants.MISSING_DATA_STRING;
				continue; 
			}else{
				boolean isbar = sgt.contains("\\|");
				String[] gts = isbar? sgt.split("\\|") : sgt.split("\\/"); 
				if(gts.length != 2 && (errGTNo++) < errRecordLimit){					
					logger.warn("invalid GT field in format column: " + vcf.toString());
					details[i] += Constants.MISSING_DATA_STRING;
					continue;
				}
					
				String[] sgd = new String[2];
				String[] alts = vcf.getAlt().split(",");
				//only allow three multi-allele, otherwise too confuse
				for(int j = 0; j < 2; j ++) 
					switch (gts[j].trim()){
						case "0": sgd[j] = vcf.getRef(); 
						break;
						case "1": sgd[j] = (alts.length >= 1) ? alts[0] : Constants.MISSING_DATA_STRING; 
						break;
						case "2": sgd[j] = (alts.length >= 2) ? alts[1] : Constants.MISSING_DATA_STRING; 
						break;
						case "3": sgd[j] = (alts.length >= 3) ? alts[2] : Constants.MISSING_DATA_STRING; 
						break;
						default:
							sgd[j] = Constants.MISSING_DATA_STRING; 
					}
				
				details[i] += isbar? String.join("|", sgd) : String.join("/", sgd);
//				newformat.add(re.toString());
			}
		}
		
		//replace GT with GT:GD actually insert GD
		List<String> newformat = new ArrayList<String>();
		for(int i = 0; i < format.size(); i ++){
			List<String> ff = Arrays.asList(format.get(i).split(Constants.COLON_STRING));
			ff.set(0, details[i]);			
			newformat.add(String.join(Constants.COLON_STRING,ff));			
		}
		
		vcf.setFormatFields(newformat);
	}
	
	/**
	 * 
	 * @return a map of, key is the indel position, value is the list a vcf record on that position. 
	 * @throws Exception
	 */
	public Map<ChrRangePosition, IndelPosition> getIndelMap() throws Exception{	
		
		Map<ChrRangePosition,IndelPosition> indelPositionMap = new  ConcurrentHashMap<ChrRangePosition,IndelPosition>();
		for(VcfRecord vcf : positionRecordMap.values()){
			ChrRangePosition indelPos = ChrRangePosition.valueOf(vcf.getChrPosition()); 
			if(indelPositionMap.containsKey(indelPos)) 
				indelPositionMap.get(indelPos).addVcf( vcf );
 			  else 
				indelPositionMap.put(indelPos, new IndelPosition(vcf));
			
//			//debug
//			System.out.println("getIndelMap():: " + indelPos.toIGVString());
//			for(VcfRecord vcfxx : indelPositionMap.get(indelPos).getIndelVcfs())
//				System.out.println(vcfxx.toString());
		}	

		return indelPositionMap;
	}
	
	public VcfHeader getVcfHeader(){ return header;	}
}
