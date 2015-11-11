package au.edu.qimr.qannotate.options;

import static java.util.Arrays.asList;
import joptsimple.OptionSet;
import au.edu.qimr.qannotate.Messages;

/*
 * parse command line to options. 
 */
public class ConfidenceOptions extends Options {
 	

    /**
     * check command line and store arguments and option information
     */
 	
 	// public Options(){  parser = new OptionParser(); this.Mode = null; } 
 	public ConfidenceOptions(){ super(Options.MODE.confidence);   }
 	
    @Override
    public boolean parseArgs(final String[] args) throws Exception{  	
    	 
        parser.acceptsAll( asList("h", "help"), Messages.getMessage("HELP_OPTION_DESCRIPTION"));
        parser.acceptsAll( asList("i", "input"), Messages.getMessage("INPUT_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("input vcf");
        parser.acceptsAll( asList("o", "output"), Messages.getMessage("OUTPUT_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("output vcf"); 
        parser.accepts("mode", "run confident mode").withRequiredArg().ofType(String.class).describedAs("dbsnp");
        parser.accepts("log", LOG_DESCRIPTION).withRequiredArg().ofType(String.class);
        parser.accepts("loglevel",  LOG_LEVEL_OPTION_DESCRIPTION).withRequiredArg().ofType(String.class);
        
        parser.accepts(test,  Messages.getMessage("TUMOUR_SAMPLEID_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("testSample");
        parser.accepts(control, Messages.getMessage("NORMAL_SAMPLEID_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("controlSample");	        

        
        final OptionSet options = parser.parse(args);   
        
        if(options.has("h") || options.has("help")){
        	displayHelp(Messages.getMessage("CONFIDENCE_USAGE"));
            return false;
        }
        
        
        if( !options.has("log")){
            System.out.println(Messages.getMessage("LOG_OPTION_DESCRIPTION"));            
            return false;
        } else{  
        	logFileName = (String) options.valueOf("log");  	
        	logLevel = (String) options.valueOf("loglevel");
        }   
        
        commandLine = Messages.reconstructCommandLine(args) ;
        
        testSample = (options.has(test))? (String)options.valueOf(test) : null;
        controlSample = (options.has(control))? (String)options.valueOf(control) : null;
 
        
        //check IO
        inputFileName = (String) options.valueOf("i") ;      	 
        outputFileName = (String) options.valueOf("o") ; 
//        databaseFileName = (String) options.valueOf("d") ; 

        final String[] inputs= (databaseFileName != null) ? new String[]{ inputFileName,databaseFileName}: new String[]{ inputFileName} ;
        final String[] outputs = new String[]{outputFileName};
        final String [] ios = new String[inputs.length + outputs.length];
        System.arraycopy(inputs, 0, ios, 0, inputs.length);
        System.arraycopy(outputs, 0, ios, inputs.length, outputs.length);
        return checkInputs(inputs )  && checkOutputs(outputs ) && checkUnique(ios);

    } 

   


   
}