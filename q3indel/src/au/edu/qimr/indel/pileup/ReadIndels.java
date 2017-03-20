package au.edu.qimr.indel.pileup;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.ChrRangePosition;

import org.qcmg.common.util.IndelUtils;
import org.qcmg.common.util.IndelUtils.SVTYPE;
import org.qcmg.common.vcf.VcfFormatFieldRecord;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.vcf.VCFFileReader;




public class ReadIndels {
	static final String FILTER_SOMATIC = "SOMATIC";
	
	QLogger logger =  QLoggerFactory.getLogger(ReadIndels.class);; 
	protected VcfHeader header; 
	
	protected final int errRecordLimit = 100;	
	//counts from each input, {No of new indel, No of overlap indel, No of indels with mult Allel, No of inputs variants, No of input variants with multi Allel}
	protected final int[] counts = {0,0, 0,0,0}; 
 
	//here key will be uniq for indel: chr, start, end, allel 
	protected final  Map<VcfRecord, VcfRecord> positionRecordMap = new  ConcurrentHashMap<VcfRecord, VcfRecord>();	
	
	public ReadIndels(  File[] files) throws IOException{ 
		if(files == null) return; 
		for(File f : files) {
			logger.info("loading vcf file: " + f.getAbsolutePath());
			readFile(f);			
		}		
	}
	
	/**
	 * load variants to hash map
	 * @param f: vcf input file
	 * @param sampleCode: replace sample code inside the input vcf file 
	 * @throws IOException
	 */
	private void readFile(File f ) throws IOException{		
		int indelNew = 0;
		int indelOverlap = 0;
		int indelMultiAltNo = 0;
    	int inLines = 0;
    	int inMultiAltNo = 0;
	
        try (VCFFileReader reader = new VCFFileReader(f)) {        	
        	header = (header == null)? reader.getHeader() : VcfHeaderUtils.mergeHeaders(header, reader.getHeader(), false);
        	
        	//no chr in front of position
			for (final VcfRecord re : reader) {	
				inLines ++;
				addGenoDetails(re);	
				
				String[] alleles =  re.getAlt().contains(",") ? re.getAlt().split(","): new String[] {re.getAlt()};				
    			if(alleles.length > 1) inMultiAltNo ++; //multi allele input variants
    			
				//for(String alt : alleles){
    			for( int i = 0; i < alleles.length; i++ ){						
					SVTYPE type = IndelUtils.getVariantType(re.getRef(), alleles[i]);
	 	        	if(type.equals(SVTYPE.DEL) ||type.equals(SVTYPE.INS) ){	
	 	        		VcfRecord vcf1 = re; 
	 	        		if(alleles.length > 1){ 	        		
		 	        		vcf1 = decomposeVcf(re, alleles[i], i );
		 	        		indelMultiAltNo ++; //mutli allele indel
	 	        		} 
	 	        		
	 	        		//somatic or germline	 	        		
     					if(positionRecordMap.containsKey(vcf1) && (indelOverlap ++) < errRecordLimit){						
 //    						logger.warn("same variants already exsits, this one will be discard:\n" + positionRecordMap.get(vcf1).toString() );
     						logger.warn("same variants already exsits, this one will be discard:\n" +  vcf1.toString() );

     						continue; //no overwrite but just warning
     					}
     					positionRecordMap.put(vcf1, vcf1);
     					indelNew ++;
   	 	        	}//done for current indel allele
				}  		 
    		}
 	      }
               
    	//counts from each input
    	counts[0] = indelNew;
    	counts[1] = indelOverlap;
    	counts[2] = indelMultiAltNo;
    	counts[3] = inLines;
    	counts[4] = inMultiAltNo;
    		                
	}
	
	private VcfRecord decomposeVcf(VcfRecord re, String alt, int altOrder){
		//set current allele
		VcfRecord re1 = new VcfRecord.Builder(re.getChrPosition(), re.getRef(), alt)
				.id(re.getId()).qualString(re.getQualString()).build();

		List<String> sformat = re.getFormatFields();
		List<String> field = new ArrayList<String>();
		field.add(sformat.get(0));
		for(int i = 1; i < sformat.size(); i ++){
			VcfFormatFieldRecord formatRe = new  VcfFormatFieldRecord(sformat.get(0), sformat.get(i));
			String gt = formatRe.getField(VcfHeaderUtils.FORMAT_GENOTYPE);
			if(altOrder == 0) gt = gt.replace("2",".").replace("3", ".");
			else if(altOrder == 1) gt = gt.replace("3", ".").replace("1", ".").replace("2","1"); 
			else  if(altOrder == 2) gt = gt.replace("2", ".").replace("1", ".").replace("3","1");
			else gt = gt.contains("/")  ? "./." : ".|."; 			
			formatRe.setField(VcfHeaderUtils.FORMAT_GENOTYPE, gt);
			field.add(formatRe.getSampleColumnString());
		}
		
		re1.setFormatFields(field);
				
		return re1; 
	}
	
	
	/**
	 * 
	 * @param vcf: add Genotyp Details into sample column.
	 *  eg. ...A T,C.. GT 0/1 1/2, become ...A T,C.. G:GD 0/1:A/T 1/1:T/C
	*/
	protected void addGenoDetails ( VcfRecord vcf ){	
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
	
	/**
	 * 
	 * @return a map of, key is the indel position, value is the list a vcf record on that position. 
	 * @throws Exception
	 */
	public Map<ChrRangePosition, IndelPosition> getIndelMap() throws Exception{	
	
		Map<ChrRangePosition,IndelPosition> indelPositionMap = new  ConcurrentHashMap<ChrRangePosition,IndelPosition>();
		for(VcfRecord vcf : positionRecordMap.values()){			
			ChrRangePosition indelPos = new ChrRangePosition(vcf.getChrPosition(), vcf.getChrPosition().getEndPosition()); 
			if(indelPositionMap.containsKey(indelPos)) 
				indelPositionMap.get(indelPos).addVcf( vcf );
				  else 
				indelPositionMap.put(indelPos, new IndelPosition(vcf));
		}	
	
		return indelPositionMap;
	}
	
	public VcfHeader getVcfHeader(){ return header;	}
	
}
