package au.edu.qimr.utility;

import java.io.File;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.meta.*;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.tab.TabbedFileWriter;
import org.qcmg.vcf.VCFFileReader;
 
public class Vcf2DCC1 {

	
	public Vcf2DCC1(Vcf2Dcc1Options option, QLogger logger) throws Exception {
		File output = new File(option.getOutputFileName());		 
		File input = new File(option.getInputFileName());		

		
		try (VCFFileReader reader = new VCFFileReader(input);
				TabbedFileWriter writer = new TabbedFileWriter(output);) {
			VcfHeader head = reader.getHeader();
			
			//read header
			 QExec exec = new  QExec(option.getPGName(), option.getVersion(), (String[]) null, option.getCommandLine(), null);
//			 public QDccMeta(String uuid, String analyzed_sample_id, String matched_sample_id, 
//						String assembly_version, String platform, String experimental_protocol, String base_calling_algorithm,
//						String alignment_algorithm, String variation_calling_algorithm, String donor_id) {
//			 QDccMeta meta = new QDccMeta();
			 
			//convert variants
			 
			 
		}
			
		
		
	}
	
	public static void main(final String[] args) throws Exception {		
		
		//check arguments
		Vcf2Dcc1Options options = new Vcf2Dcc1Options( args);	
		if(! options.commandCheck()){ System.out.println("command check failed! ");	System.exit(1);	}		
		
		 
		QLogger logger =  options.getLogger(args, Vcf2DCC1.class);	
		
		try{    
			Vcf2DCC1 convert = new Vcf2DCC1(options, logger);		
			
		
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
