package au.edu.qimr.utility;

import static java.util.Arrays.asList;

import java.io.File;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;

import au.edu.qimr.utility.Messages;
/*
 * parse command line to options. 
 */
public abstract class Options {
	
   protected static final String VERSION_DESCRIPTION = Messages.getMessage("VERSION_OPTION_DESCRIPTION");	 
   protected static final String HELP_DESCRIPTION = Messages.getMessage("HELP_OPTION_DESCRIPTION");  
    
   protected static final String LOG_DESCRIPTION = Messages.getMessage("LOG_OPTION_DESCRIPTION");
   protected static final String LOG_LEVEL_OPTION_DESCRIPTION = Messages.getMessage("LOG_LEVEL_OPTION_DESCRIPTION");
   protected static final String USAGE = Messages.getMessage("USAGE");

   protected static final String version = au.edu.qimr.utility.Main.class.getPackage().getImplementationVersion();	

   protected final String input = "input";
   protected final String output = "output";

  
   
//   protected static final String test = "test";
//   protected static final String control = "control";
   protected final OptionParser parser = new OptionParser();  
	protected OptionSet options;	

	protected  boolean commandCheck = false;
    protected String commandLine;

		
	protected String outputFileName = null;
	protected String inputFileName = null;
    
//	protected String logFileName = null;
//	protected  String logLevel;  
	protected QLogger logger = null;
	public abstract boolean commandCheck() throws Exception;
    
    public void parseArgs(final String[] args, Class myclass) throws Exception{ 

       	parser.allowsUnrecognizedOptions(); 
        parser.acceptsAll( asList("h", "help"), HELP_DESCRIPTION);
        parser.acceptsAll( asList("v", "version"), VERSION_DESCRIPTION);
        options  = parser.parse(args);   
        commandLine = Messages.reconstructCommandLine(args) ;
        
//        if(options.has("v") || options.has("version")){
//            System.out.println( "Current version is " + getVersion());
//            return false;
//        }
               
//		logger = getLogger(args,myclass); 				         
//        commandCheck = commandCheck(); 
//       
      
	} 
    
    public void displayHelp() throws Exception {
		    System.out.println(Messages.getMessage("USAGE"));  
		    Thread.sleep(1);
		    parser.printHelpOn(System.err);	        
    }
    
    public String getVersion(){
    	return version;
    }
    
    public String getPGName(){
        return Messages.getProgramName();
    }
    

//    protected boolean checkUnique(String[] ios ){   
//    	
//    	final File[] fios = new File[ios.length];
// //   	Path[] pios = new Path[ios.length];
//    	
//    	for (int i = 0; i < ios.length; i ++)
//    		fios[i] = new File(ios[i]); 
//    	   	
//    	for(int  i = ios.length -1; i > 0; i --)
//    		for (int j = i-1; j >= 0; j -- )
//				try {
//					//if( Files.isSameFile(pios[i], pios[j]))
//					if(fios[i].getCanonicalFile().equals(fios[j].getCanonicalFile()))
//						throw new Exception( "below command line values are point to same file: \n\t" + fios[i] + "\n\t" + fios[j]   ) ;
//				} catch (final Exception e) {
//					//e.printStackTrace();
//					System.err.println(e.getMessage());
//					return false;
//				}
//
//    	return true;    	
//    }	 
    
//    protected boolean checkOutputs( String[]  outputs){   
//    	for(int i = 0; i < outputs.length; i ++){
//	        final File out = new File(outputs[i] );
//	        if((out.exists() && !out.canWrite()) || !out.getParentFile().canWrite() ){
//	        	System.err.println( Messages.getMessage("OUTPUT_ERR_DESCRIPTION",out.getName()) );
//	        	 return false;
//	        }
//	        
//    	}
//    	return true;
//    }
   
	public String getCommandLine() {	return commandLine; }		
	
	 boolean hasHelp() throws Exception{
		 if(options.has("h") || options.has("help")){    
        	System.out.println(USAGE);			 
        	System.out.println(HELP_DESCRIPTION);
            parser.printHelpOn(System.err);
            return true;       
        }
		return false;
	 }
	 
	 boolean hasVersion()throws Exception{
		 if(options.has("v") || options.has("version")){
			System.out.println(VERSION_DESCRIPTION);
			System.err.println(version);    	
			return true;       
	    }
		return false;
	 }
	 
	 
    /**
     * check input and output files
     * @return true if input file readable and output file writable
     */
    protected boolean checkInputs( String[] inputs ){   
        String errMessage = null;
        
        for(int i = 0; i < inputs.length; i ++){
        	final File in = new File(inputs[i] );
	        if(!in.exists()) 
	            errMessage = Messages.getMessage("NONEXIST_INPUT_FILE", in.getPath());
	        else if(!in.isFile())       
	            errMessage = Messages.getMessage("FILE_NOT_DIRECTORY", in.getPath());
	         else if(!in.canRead())
	        	errMessage = Messages.getMessage("UNREAD_INPUT_FILE",in.getPath());     
	           
	        if(errMessage != null){
	        	System.err.println(errMessage);
	        	 return false;
	        }
	       }  	
	        	
	       return true;
    }	
    
    protected boolean checkOutputs( String[] outputs ){
       String errMessage = null;
        
        for(int i = 0; i < outputs.length; i ++){
        	final File in = new File(outputs[i] );
	        if(!in.exists()) 
	            errMessage = Messages.getMessage("NONEXIST_INPUT_FILE", in.getPath());
	        else if(!in.isFile())       
	            errMessage = Messages.getMessage("FILE_NOT_DIRECTORY", in.getPath());
	         else if(!in.canRead())
	        	errMessage = Messages.getMessage("UNREAD_INPUT_FILE",in.getPath());     
	           
	        if(errMessage != null){
	        	System.err.println(errMessage);
	        	 return false;
	        }
	       }  	
	        	
	       return true;    	
    }
    
    
    protected boolean checkUnique(String[] ios ){   
		
		final File[] fios = new File[ios.length];
		
		for (int i = 0; i < ios.length; i ++)
			fios[i] = new File(ios[i]); 
		   	
		for(int  i = ios.length -1; i > 0; i --)
			for (int j = i-1; j >= 0; j -- )
				try {
					//if( Files.isSameFile(pios[i], pios[j]))
					if(fios[i].getCanonicalFile().equals(fios[j].getCanonicalFile()))
						throw new Exception( "below command line values are point to same file: \n\t" + fios[i] + "\n\t" + fios[j]   ) ;
					} catch (final Exception e) {
						//e.printStackTrace();
					System.err.println(e.getMessage());
					return false;
				}
	
		return true;    	
	}	    
    
    /**
     * check input and output files
     * @return true if input file readable and output file writable
     */
    protected boolean checkIO(String input, String output ){   	
    	final File in = new File(input );
        final File out = new File(output  );
               
		if(! in.getAbsolutePath().equals(out.getAbsoluteFile())){
			System.err.println(Messages.getMessage("SAME_FILES", input, output));		
			return false; 
		}
                
        return   checkInputs(new String[] {input})  && checkOutputs(new String[] {input});        
    }
	 
    //IO parameters
	String getUniqIO(String io) throws Exception{		

		int size = options.valuesOf(io).size();
		if( size > 1){
			throw new Exception("multiple "+ io + " files specified" );
		}
		else if( size < 1 ){
			throw new Exception(" missing or invalid IO option specified: " + io );
		}
		
		return options.valueOf(io).toString();		 
	}
	
	QLogger getLogger(String[] args, Class myclass) throws Exception{
		
		// configure logging
		if(logger != null)
			return logger; 
		
		String logLevel = (String) options.valueOf("loglevel");
		String logFile;
		if(options.has("log")){
			logFile = options.valueOf("log").toString();			 
		}else{
			logFile = options.valueOf(output) + ".log";			 
		}		
		
		logger = QLoggerFactory.getLogger( myclass, logFile,logLevel);
		logger.logInitialExecutionStats(myclass.toString(), version, args);	
		return logger;
	}
}
