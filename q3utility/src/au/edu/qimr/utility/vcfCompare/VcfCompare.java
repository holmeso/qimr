package au.edu.qimr.utility.vcfCompare;



import htsjdk.tribble.readers.TabixReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.qcmg.common.log.*;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.vcf.VCFFileReader;
import org.qcmg.vcf.VCFFileWriter;
import org.qcmg.common.vcf.header.*;
public class VcfCompare {
	
	protected final Map<ChrPosition,VcfRecord> positionRecordMap = new HashMap<ChrPosition,VcfRecord>();
//	protected final Map<ChrPosition, Integer> positionRecordMap = new HashMap<ChrPosition,Integer>();
	protected VcfHeader header  = null;
//	protected VcfHeader header2 = null;
	
	VcfCompare(File primary, File secondary) throws Exception {
		
		//read record into RAM, meanwhile wipe off the ID field value;
        try ( VCFFileReader reader1 = new VCFFileReader(primary);
        		VCFFileReader reader2 = new VCFFileReader(secondary)) {
        	String[] firstSamples = new String[] {null}; 
        	String[] secondSamples = new String[] {null}; 
        	
        	if(reader1.getHeader().getSampleId() == null)
        		firstSamples = reader1.getHeader().getSampleId();
        	
           	if(reader2.getHeader().getSampleId() == null)
        		secondSamples = reader2.getHeader().getSampleId();
           	          	
    		if(firstSamples.length != secondSamples.length)
    			throw new Exception();
        		       		
			for(int i = 0; i < firstSamples.length; i ++)
				if(! firstSamples[i].toLowerCase().equals(secondSamples[i].toLowerCase()))
					throw new Exception();
        		        	
        	//merge header
        	header = VcfHeaderUtils.mergeHeaders(reader1.getHeader(), reader2.getHeader(), false);
        	
        	//load first file
        	for (final VcfRecord re : reader1){         		
        		ChrPosition pos = re.getChrPosition();
        		re.appendInfo("vcfFrom=f1");
        		positionRecordMap.put(new ChrPosition(pos.getChromosome(), pos.getPosition(), pos.getEndPosition(), re.getAlt()), re);	 
        	}
        	
        	//compare to second file
        	for (final VcfRecord re : reader2){
        		ChrPosition pos = re.getChrPosition();
        		pos = new ChrPosition(pos.getChromosome(), pos.getPosition(), pos.getEndPosition(), re.getAlt());        		
        		if(positionRecordMap.get(pos) == null){  
        			re.appendInfo("vcfFrom=");
        			re.setInfo(re.getInfo() + "vcfFrom=f2");
        			positionRecordMap.put(pos, re); //second file
        		}else{ 
        			re.setInfo(re.getInfo() + "vcfFrom=both");
        			positionRecordMap.put(pos, re); //both file       	
        		}	
        	}         	
		} 
//        	TabixReader  t1 = new TabixReader( primary.getAbsolutePath() );
//			TabixReader  t2 = new TabixReader( secondary.getAbsolutePath()) ; 
	}
	
	void VcfMerge( File output) throws IOException{
		//rehead
		//merge
		
		
		final List<ChrPosition> orderedList = new ArrayList<ChrPosition>(positionRecordMap.keySet());
		Collections.sort(orderedList);
		
		try(VCFFileWriter writer = new VCFFileWriter( output)) {			
			for(final VcfHeader.Record record: header)  {
				writer.addHeader(record.toString());
			}
			long count = 0; 
			for (final ChrPosition position : orderedList) {				
				VcfRecord record = positionRecordMap.get(position); 
				writer.add( record );	
				count ++;
			}
		}  
	}
	
			
	public static void main(String[] args) throws Exception{
		//check arguments
		Options options = new Options( args);	
		if(! options.commandCheck()){ 	System.exit(1);	}		
		
		QLogger logger =  options.getLogger(args);		
		try{     	
			File foutput = new File(options.getIO("output"));
			File primary = new File(options.getIO("primaryVcf"));			
			File addition = new File(options.getIO("additionalVcf"));			 
			
			VcfCompare compare = new VcfCompare( primary, addition);
			logger.info("Merging VCF ...");	 		
			compare.VcfMerge( foutput);
			
			logger.info("outputed VCF   " );
			//count reads number in each window and output  
        	logger.logFinalExecutionStats(0);
        	System.exit(0);
        }catch(Exception e){ 
        	logger.error(e.toString());
        	logger.logFinalExecutionStats(1);
        	System.err.println(e.toString());       	 
        	System.exit(1);
        }		
	}
	
}
