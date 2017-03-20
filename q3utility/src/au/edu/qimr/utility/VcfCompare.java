package au.edu.qimr.utility;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import org.qcmg.common.log.*;
import org.qcmg.common.meta.QExec;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.model.ChrPositionName;
import org.qcmg.common.util.Constants;
import org.qcmg.common.vcf.*;
import org.qcmg.common.vcf.header.*;
import org.qcmg.vcf.*;

public class VcfCompare {
	
	
	final Map<ChrPositionName,VcfRecord> positionRecordMap = new HashMap<ChrPositionName,VcfRecord>();
	VcfHeader header  = null;	
	private long noPrimary = 0;
	private long noAdditional = 0;
	private long noBoth = 0; 
	private long noOutput = 0;
	 	
	public VcfCompare(File primary, File secondary) throws IOException,RuntimeException  {		
		//read record into RAM, meanwhile wipe off the ID field value;
        try ( VCFFileReader reader1 = new VCFFileReader(primary);
        		VCFFileReader reader2 = new VCFFileReader(secondary)) {
        	String[] firstSamples = new String[] {null}; 
        	String[] secondSamples = new String[] {null};         	
       	
        	if(reader1.getHeader().getSampleId() != null)
        		firstSamples = reader1.getHeader().getSampleId();
        	
           	if(reader2.getHeader().getSampleId() != null)
        		secondSamples = reader2.getHeader().getSampleId();
          	          	
    		if(firstSamples.length != secondSamples.length)
    			throw new RuntimeException(firstSamples.length + " != "  + secondSamples.length);
    		        		       		
			for(int i = 0; i < firstSamples.length; i ++)
				if( firstSamples[i] != null &&  secondSamples[i] != null 
					&& ! firstSamples[i].toLowerCase().equals(secondSamples[i].toLowerCase()))
					throw new RuntimeException(firstSamples[i] + "!=" + secondSamples[i]);
        		        	
         	//load first file
        	for (final VcfRecord re : reader1){     
        		noPrimary ++; 
        		ChrPosition pos = re.getChrPosition();
        		re.appendInfo(VcfCompareOptions.Info_From + "=1");        		
        		positionRecordMap.put(new ChrPositionName(pos.getChromosome(),pos.getStartPosition(), pos.getEndPosition(), re.getAlt()), re);	 
        	}
        	
        	//compare to second file
        	for (final VcfRecord re : reader2){
        		noAdditional ++;
        		ChrPositionName pos =  new ChrPositionName(re.getChromosome(),re.getPosition(), re.getChrPosition().getEndPosition(), re.getAlt());        		
        		if(positionRecordMap.get(pos) == null)  //second file
        			re.appendInfo(VcfCompareOptions.Info_From + "=2"); 
        		else{ //both file
        			noBoth ++; 
        			re.appendInfo(VcfCompareOptions.Info_From + "=0"); 
        		}
        		positionRecordMap.put(pos, re); //replace
        	}         	
		} 
	}
	
	public void reheader( String cmd, VcfCompareOptions options) throws Exception {	
		DateFormat df = new SimpleDateFormat("yyyyMMdd");
		String version = VcfCompare.class.getPackage().getImplementationVersion();
		String pg = VcfCompare.class.getPackage().getImplementationTitle();
		final String fileDate = df.format(Calendar.getInstance().getTime());
		final String uuid = QExec.createUUid();		

		try ( VCFFileReader reader1 = new VCFFileReader(options.getIO(VcfCompareOptions.primaryInput));
				VCFFileReader reader2 = new VCFFileReader(options.getIO(VcfCompareOptions.additionalInput))) {
		   
			VcfHeader h2 = reader2.getHeader();
			String inputUuid = (h2.getUUID() == null)? null: h2.getUUID().getMetaValue();   
		//	h2.replace("##" + Options.additionalInput + "=" + inputUuid + ":" + options.getIO(Options.additionalInput));
			h2.addOrReplace( VcfHeaderUtils.STANDARD_INPUT_LINE  + "=" + inputUuid + ":" + options.getIO(VcfCompareOptions.additionalInput));
			
			VcfHeader h1 = reader1.getHeader();	    	   
			inputUuid = (h1.getUUID() == null)? null: h1.getUUID().getMetaValue();  
			h1.addOrReplace( VcfHeaderUtils.STANDARD_INPUT_LINE  + "=" + inputUuid + ":" + options.getIO(VcfCompareOptions.primaryInput));
			h1.addOrReplace(VcfHeaderUtils.STANDARD_FILE_DATE + "=" + fileDate);
			h1.addOrReplace(VcfHeaderUtils.STANDARD_UUID_LINE + "=" + uuid);
			h1.addOrReplace(VcfHeaderUtils.STANDARD_SOURCE_LINE + "=" + pg+"-"+version);	
			
			header = VcfHeaderUtils.mergeHeaders(reader1.getHeader(), reader2.getHeader(), false);
		}
	 		
		if(version == null) version = Constants.NULL_STRING_UPPER_CASE;
	    if(pg == null ) pg = Constants.NULL_STRING_UPPER_CASE;
	    if(cmd == null) cmd = Constants.NULL_STRING_UPPER_CASE;
		VcfHeaderUtils.addQPGLineToHeader(header, pg, version, cmd);
		header.addInfo(VcfCompareOptions.Info_From, "1", "Integer", VcfCompareOptions.Info_From_Description);
	
	}
	
	public void VcfMerge( File output) throws IOException{
		final List<ChrPositionName> orderedList = new ArrayList<ChrPositionName>(positionRecordMap.keySet());
		Collections.sort(orderedList);
		
		try(VCFFileWriter writer = new VCFFileWriter( output)) {			
			for(final VcfHeaderRecord record: header)  {
				writer.addHeader(record.toString());
			}
		
			for (final ChrPosition position : orderedList) {				
				VcfRecord record = positionRecordMap.get(position); 
				writer.add( record );	
				noOutput ++; 
			}
		}  
	}
	
	public long getCountPrimaryOnly() {return noPrimary - noBoth; }
	public long getCountAddtionalOnly(){return noAdditional - noBoth; }
	public long getCountBoth(){return noBoth; }
	public long getCountOutput(){return noOutput;}
	
	public static void main(final String[] args) throws Exception {		
	
		//check arguments
		VcfCompareOptions options = new VcfCompareOptions( args);	
		if(! options.commandCheck()){ System.out.println("command check failed! ");	System.exit(1);	}		
		
		 
		QLogger logger =  options.getLogger(args);		
		try{    
			File foutput = new File(options.getIO(VcfCompareOptions.output));
			File primary = new File(options.getIO("primaryInput"));			
			File addition = new File(options.getIO("additionalInput"));		
						
			VcfCompare compare = new VcfCompare( primary, addition);
			logger.info("Merging VCF ...");	 	
			compare.reheader( Messages.reconstructCommandLine(args), options);

			compare.VcfMerge( foutput);
			logger.info("Total outputed variant record is " + compare.getCountOutput() );
			logger.info("Including variant record appeared from both input: " + compare.getCountBoth() );
			logger.info("Including variant record appeared from primary input: " + compare.getCountPrimaryOnly() );
			logger.info("Including variant record appeared from additional input: " + compare.getCountAddtionalOnly() );
			
			if(compare.getCountOutput() != 
					(compare.getCountBoth() + compare.getCountPrimaryOnly() + compare.getCountAddtionalOnly()) )
				logger.warn("The total output variants number is not same as total number from both inputs! ");			 
			
			//count reads number in each window and output  
			logger.info(  getUsedMemory());
        	logger.logFinalExecutionStats(0);
        	System.exit(0);
        }catch(Exception e){ 
        	logger.error(e.toString());
        	logger.logFinalExecutionStats(1);
        	System.err.println(e.toString());       	 
        	System.exit(1);
        }		
	}
	
	static String getUsedMemory(){
		final long MEGABYTE = 1024L * 1024L;

		Runtime runtime = Runtime.getRuntime();
		runtime.gc();	
		long MTotal =  runtime.totalMemory() / MEGABYTE;
//		long MUsed = MTotal - runtime.freeMemory() / MEGABYTE;
		long MMax = runtime.maxMemory() / MEGABYTE;		 
		
		return String.format("ResourcesUsed: mem=%dM; vmem=%dM.", MTotal,MMax);
	}
	
}
