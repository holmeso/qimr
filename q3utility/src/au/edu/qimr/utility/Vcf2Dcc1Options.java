package au.edu.qimr.utility;

import java.util.List;


public class Vcf2Dcc1Options extends Options{
	private final String input = "VcfInput";
	private final String output = "DccOutput";
	
//	private final QLogger logger ;
	
	private static final String OUTPUT_DESCRIPTION = Messages.getMessage("DCC_OUTPUT_OPTION_DESCRIPTION");
	private static final String INTPUT_DESCRIPTION = Messages.getMessage("VCF_INPUT_OPTION_DESCRIPTION");
 	
	Vcf2Dcc1Options(final String[] args) throws Exception{
		
		parser.accepts(output, OUTPUT_DESCRIPTION).withRequiredArg().ofType(String.class).describedAs("dcc output file");
		parser.accepts(input, INTPUT_DESCRIPTION).withRequiredArg().ofType(String.class).describedAs("vcf input vcf");
				
		parseArgs( args, au.edu.qimr.utility.Vcf2DCC1.class);
		commandCheck();
	}
	 
	 @Override
	public boolean commandCheck() throws Exception{
			//quit system after provide help or version info
			if( hasHelp() || hasVersion() ){ 
				 System.exit(0); 
			} 		 	
			
			if (options.nonOptionArguments().size() > 0) {	
				@SuppressWarnings("unchecked")
				List<String> nonoptions = (List<String>) options.nonOptionArguments();

				for(String str : nonoptions){
					System.err.println("INVALID OPTION: " + str);
				} 				
				System.exit(1); 
			}
		
			if( ! checkIO(getUniqIO(input), getUniqIO(output))){
				System.exit(1);
			}
			
			return true;				
		}
	 
	 public String getInputFileName() throws Exception{
		 
		 return getUniqIO(input); 
	 }
	 
	 public String getOutputFileName() throws Exception{
		 
		 return getUniqIO(output); 
	 }
 
 
}
