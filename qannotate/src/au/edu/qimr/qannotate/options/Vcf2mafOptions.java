package au.edu.qimr.qannotate.options;

import static java.util.Arrays.asList;
import joptsimple.OptionSet;
import au.edu.qimr.qannotate.Messages;


public class Vcf2mafOptions extends Options {
	 public Vcf2mafOptions( ) {  super(Options.MODE.vcf2maf);	  }
	 
	 public static final String unkown = "unkown";
	 public static final String null_String = "null";
	 public static final String default_center = "QIMR_Berghofer"; 
	 
	 String center ; 
	 String sequencer; 
	 String outputDir;
	 String donorId; 
	 boolean appendSample = false; 
	 boolean lowMaf = false; 

	 
	 
	 final String  Description_sequencer = "eg.  <Illumina GAIIx, Illumina HiSeq,SOLID,454, ABI 3730xl, Ion Torrent PGM,Ion Torrent Proton,PacBio RS, Illumina MiSeq,Illumina HiSeq 2500,454 GS FLX Titanium,AB SOLiD 4 System>";
	 
	 @Override
	    public boolean parseArgs(final String[] args) throws Exception{  	
	    	 
	        parser.acceptsAll( asList("h", "help"), Messages.getMessage("HELP_OPTION_DESCRIPTION"));
	        parser.acceptsAll( asList("i", "input"), Messages.getMessage("INPUT_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("input vcf");
	        parser.acceptsAll( asList("o", "output"),  "output maf file with full path").withRequiredArg().ofType(String.class).describedAs("output maf"); 
	        parser.accepts("outdir", Messages.getMessage("MAF_OUTPUT_DIRECTORY_OPTION_DESCRIPTION") ).withRequiredArg().ofType(String.class).describedAs("output file location");
	        parser.accepts("donor",  Messages.getMessage("DONOR_ID_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("donor id");
	        
	        parser.accepts(test,  Messages.getMessage("TUMOUR_SAMPLEID_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("testSample");
	        parser.accepts(control, Messages.getMessage("NORMAL_SAMPLEID_DESCRIPTION")).withRequiredArg().ofType(String.class).describedAs("controlSample");	        
	        parser.accepts("center", "Genome sequencing center").withRequiredArg().ofType(String.class).describedAs("center");
	        parser.accepts("sequencer", Description_sequencer).withRequiredArg().ofType(String.class).describedAs("Sequencer");	 
	        parser.accepts("lowMaf", Messages.getMessage("LOW_MAF_DESCRIPTION"));
	        
	        parser.accepts("mode", "run vcf2maf").withRequiredArg().ofType(String.class).describedAs("vcf2maf");
	       // "(compulsary) database location"
 	        parser.accepts("log", LOG_DESCRIPTION).withRequiredArg().ofType(String.class);
	        parser.accepts("loglevel",  LOG_LEVEL_OPTION_DESCRIPTION).withRequiredArg().ofType(String.class);
//	        parser.accepts("sampleFlag",  Messages.getMessage("SAMPLE_FLAG_DESCRIPTION")).withRequiredArg().ofType(String.class);
	        final OptionSet options = parser.parse(args);   
	        
	        if(options.has("h") || options.has("help")){
	        	displayHelp(Messages.getMessage("VCF2MAF_USAGE"));
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
	        //check IO
	        inputFileName = (String) options.valueOf("i") ;      	 
	        outputFileName = (String) options.valueOf("o") ; 
	        outputDir = (options.has("outdir"))? (String)options.valueOf("outdir") : null;
	        center = (options.has("center"))? (String)options.valueOf("center") : null;
	        sequencer = (options.has("sequencer"))? (String)options.valueOf("sequencer") : null;	 
	        
	        testSample = (options.has(test))? (String)options.valueOf(test) : null;
	        controlSample = (options.has(control))? (String)options.valueOf(control) : null;
	        donorId = (options.has("donor"))? (String)options.valueOf("donor") : null;
	        
	        if( options.has("lowMaf"))	lowMaf = true;
	        
  	        return true;
	     } 
	 
	 public String getCenter(){  return center; }
	 public String getSequencer(){  return sequencer; }
	 public String getOutputDir(){return outputDir;}
	 public String getDonorId(){return donorId; }
	 public boolean doOutputLowMaf(){ return lowMaf; }
	 
}