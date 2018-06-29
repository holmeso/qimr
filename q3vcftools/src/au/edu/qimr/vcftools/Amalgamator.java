package au.edu.qimr.vcftools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.qcmg.common.model.ChrPositionName;
import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.meta.QExec;
import org.qcmg.common.util.Constants;
import org.qcmg.common.util.FileUtils;
import org.qcmg.common.util.LoadReferencedClasses;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.VcfUtils;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.vcf.VCFFileReader;


public class Amalgamator {
	
	private static QLogger logger;
	private QExec exec;
	private String[] vcfFiles;
	private String outputFileName;
	private String version;
	private String logFile;
	private int exitStatus;
	
	private final Map<ChrPositionName, String[][]> positions = new HashMap<>();
	
	
	protected int engage() throws IOException {
			
		logger.info("about to load vcf files");
		loadVcfs();
		
		writeOutput();
		return exitStatus;
	}
	
	private void writeOutput() throws FileNotFoundException {
		List<ChrPositionName> recs = new ArrayList<>(positions.keySet());
		Collections.sort(recs);
		
		logger.info("writing output");
		
		
		try (PrintStream ps = new PrintStream(new FileOutputStream(new File(outputFileName)))) {
			/*
			 * put input files along with their positions into the file header
			 */
			int j = 1;
			for (String s : vcfFiles) {
				ps.println("##" + j++ + ": vcf file: " + s);
			}
			
			
			String header = "#chr\tposition\tref\talt";
			for (int i = 1 ; i <= vcfFiles.length ; i++) {
				header += "\tGT:" + i;
			}
			for (int i = 1 ; i <= vcfFiles.length ; i++) {
				header += "\tAC:" + i;
			}
			ps.println(header);
			for (ChrPositionName cpn : recs) {
				String[][] arrs = positions.get(cpn);
				String gts = "";
				String covs = "";
				for(String[] arr : arrs) {
					gts += "\t" + (null != arr[0] ? arr[0] : "./.");
					covs += "\t" + (null != arr[1] ? arr[1] : ".");
				}
				ps.println(cpn.getChromosome() + "\t" + cpn.getStartPosition() + "\t" + cpn.getName() + gts + covs);
			}
		}
		
		logger.info("writing output- DONE");
		
	}
	
	private void loadVcfs() throws IOException {
		int i = 0;
		int fileCount = vcfFiles.length;
		int index = 0;
		for (String s : vcfFiles) {
			
			try (VCFFileReader reader = new VCFFileReader(new File(s))) {
				for (VcfRecord rec : reader) {
					if (++ i % 1000000 == 0) {
						logger.info("hit " + i + " entries");
					}
					
					/*
					 * we want HC or PASS
					 */
					if (isRecordHighConfOrPass(rec)) {
						
						/*
						 * we want germline (for now - need to optionalise)
						 */
						if ( ! isRecordSomatic(rec)) {
							/*
							 * only deal with snps and compounds for now
							 */
							if (VcfUtils.isRecordASnpOrMnp(rec)) {
								
								String ref = rec.getRef();
								String alt = rec.getAlt();
								
								if (ref.length() > 1) {
									/*
									 * split cs into constituent snps
									 */
									for (int z = 0 ; z < ref.length() ; z++) {
										
										ChrPositionName cpn  = new ChrPositionName(rec.getChrPosition().getChromosome(), rec.getChrPosition().getStartPosition() + z,  rec.getChrPosition().getStartPosition() + z , ref.charAt(z) + "\t" + alt.charAt(z));
										String [][] arr = positions.computeIfAbsent(cpn, v -> new String[fileCount][2]);
										
										/*
										 * most likely don't have GT and AC info for cs
										 */
										arr[index][0] =  "C/S";
										arr[index][1] =  "C,S,0";
									}
								} else {
								
								
									ChrPositionName cpn  = new ChrPositionName(rec.getChrPosition().getChromosome(), rec.getChrPosition().getStartPosition(),  rec.getChrPosition().getStartPosition() + ref.length() - 1 , ref + "\t" + alt);
									String [][] arr = positions.computeIfAbsent(cpn, v -> new String[fileCount][2]);
									
									List<String> ffList = rec.getFormatFields();
									/*
									 * get position of GT from first entry, and then get second (germline and qsnp)
									 */
									String [] formatHeaders = ffList.get(0).split(":");
									int j = 0;
									int gtPosition = 0;
									int oabsPosition = 0;
									
									for (String h : formatHeaders) {
										if (VcfHeaderUtils.FORMAT_GENOTYPE.equals(h)) {
											gtPosition = j;
										} else  if (VcfHeaderUtils.FORMAT_ALLELE_COUNT.equals(h)) {
											oabsPosition = j;
										}
										j++;
									}
									String [] params = ffList.get(1).split(":"); 
									String gts = params[gtPosition];
									/*
									 * this could contain the ampesand - if so, get first (qsnp) element
									 */
									int ampesandIndex = gts.indexOf(Constants.VCF_MERGE_DELIM);
									if (ampesandIndex > -1) {
										gts = gts.substring(0, ampesandIndex);
									}
									arr[index][0] =  gts;
									/*
									 * get allele dist next
									 */
									String oabs = params[oabsPosition];
									ampesandIndex = oabs.indexOf(Constants.VCF_MERGE_DELIM);
									if (ampesandIndex > -1) {
										oabs = oabs.substring(0, ampesandIndex);
									}
//									logger.info("attempting to get allele dist for ac: " + oabs);
									Map<String, Integer> alleleDist = VcfUtils.getAllelicCoverageFromAC(oabs);
									arr[index][1] =  getAllelicDistFromMap(alleleDist, ref, alt);
									
									
								}
							}
						}
					}
				}
			}
			logger.info("input: " + (index+1) + " has " + i + " entries");
			logger.info("positions size: " + positions.size());
			i = 0;
			index++;
		}
		logger.info("Number of positions to be reported upon: " + positions.size());
	}
	
	/*
	 * Returns a string representing the coverage counts of ref, alt and rest based on entries in the supplied map
	 */
	static String getAllelicDistFromMap(Map<String, Integer> map, String ref, String alt) {
		int refCount = map.computeIfAbsent(ref, v -> Integer.valueOf(0)).intValue();
		int altCount = map.computeIfAbsent(alt, v -> Integer.valueOf(0)).intValue();
		int restCount = map.entrySet().stream().filter(kv -> ! kv.getKey().equals(ref) && ! kv.getKey().equals(alt)).mapToInt(kv -> kv.getValue().intValue()).sum();
		
		return refCount + Constants.COMMA_STRING + altCount + Constants.COMMA_STRING + restCount; 
	}
	
	public static final boolean isRecordHighConfOrPass(VcfRecord r) {
		if (r.getInfo().contains("CONF=HIGH_1,HIGH_2")) {
			return true;
		}
		List<String> l = r.getFormatFields();
		/*
		 * all fields in list ( apart from first) should have PASS
		 */
		long numberOfPasses = l.stream().filter(s -> s.contains("PASS")).count();
		
		return numberOfPasses == l.size() - 1;
	}
	
	public static final boolean isRecordSomatic(VcfRecord r) {
		if (VcfUtils.isRecordSomatic(r)) {
			return true;
		}
		List<String> l = r.getFormatFields();
		/*
		 * all fields in list ( apart from first) should have PASS
		 */
		long numberOfSomatics = l.stream().filter(s -> s.contains("SOMATIC")).count();
		
		return numberOfSomatics == (l.size() - 1) / 2;
	}
	
	public static void main(String[] args) throws Exception {
		//loads all classes in referenced jars into memory to avoid nightly build sheninegans
		LoadReferencedClasses.loadClasses(Amalgamator.class);
		
		Amalgamator qp = new Amalgamator();
		int exitStatus = qp.setup(args);
		if (null != logger) {
			logger.logFinalExecutionStats(exitStatus);
		} else {
			System.err.println("Exit status: " + exitStatus);
		}		
		System.exit( exitStatus );
	}

	private int setup(String[] args) throws Exception {
		int returnStatus = 1;
		Options options = new Options(args);

		if (options.hasHelpOption()) {
			System.err.println(Messages.USAGE);
			options.displayHelp();
			returnStatus = 0;
		} else if (options.hasVersionOption()) {
			System.err.println(Messages.getVersionMessage());
			returnStatus = 0;
		} else if (options.getVcfs().length < 1) {
			System.err.println(Messages.USAGE);
		} else {
			// configure logging
			logFile = options.getLog();
			version = Amalgamator.class.getPackage().getImplementationVersion();
			if (null == version) {	version = "local"; }
			logger = QLoggerFactory.getLogger(Amalgamator.class, logFile, options.getLogLevel());
			exec = logger.logInitialExecutionStats("q3vcftools Amalgamator", version, args);
			
			// get list of file names
			vcfFiles = options.getVcfs();
			if (vcfFiles.length < 1) {
				throw new Exception("INSUFFICIENT_ARGUMENTS");
			} else {
				// loop through supplied files - check they can be read
				for (int i = 0 ; i < vcfFiles.length ; i++ ) {
					if ( ! FileUtils.canFileBeRead(vcfFiles[i])) {
						throw new Exception("INPUT_FILE_ERROR: "  +  vcfFiles[i]);
					}
				}
			}
			
			// set outputfile - if supplied, check that it can be written to
			if (null != options.getOutputFileName()) {
				String optionsOutputFile = options.getOutputFileName();
				if (FileUtils.canFileBeWrittenTo(optionsOutputFile)) {
					outputFileName = optionsOutputFile;
				} else {
					throw new Exception("OUTPUT_FILE_WRITE_ERROR");
				}
			}
			
			logger.info("vcf input files: " + Arrays.deepToString(vcfFiles));
			logger.info("outputFile: " + outputFileName);
			
			return engage();
		}
		return returnStatus;
	}

}