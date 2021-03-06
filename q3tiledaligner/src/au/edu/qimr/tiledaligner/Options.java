/**
 * © Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
/**
 * All source code distributed as part of the AdamaJava project is released
 * under the GNU GENERAL PUBLIC LICENSE Version 3, a copy of which is
 * included in this distribution as gplv3.txt.
 */
package au.edu.qimr.tiledaligner;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

final class Options {

	private static final String HELP_DESCRIPTION = Messages.getMessage("HELP_OPTION_DESCRIPTION");
	private static final String VERSION_DESCRIPTION = Messages.getMessage("VERSION_OPTION_DESCRIPTION");
	private static final String NO_OF_THREADS_OPTION_DESCRIPTION = Messages.getMessage("NO_OF_THREADS_OPTION_DESCRIPTION");
	private static final String INCLUDE_OPTION_DESCRIPTION = Messages.getMessage("INCLUDE_OPTION_DESCRIPTION");
	private static final String TAGS_OPTION_DESCRIPTION = Messages.getMessage("TAGS_OPTION_DESCRIPTION");
	private static final String TAGS_INT_OPTION_DESCRIPTION = Messages.getMessage("TAGS_INT_OPTION_DESCRIPTION");
	private static final String TAGS_CHAR_OPTION_DESCRIPTION = Messages.getMessage("TAGS_CHAR_OPTION_DESCRIPTION");
	private static final String MAX_RECORDS_OPTION_DESCRIPTION = Messages.getMessage("MAX_RECORDS_OPTION_DESCRIPTION");
	private static final String LOG_OPTION_DESCRIPTION = Messages.getMessage("LOG_OPTION_DESCRIPTION");
	private static final String MIN_BIN_SIZE_OPTION_DESCRIPTION = Messages.getMessage("MIN_BIN_SIZE_OPTION_DESCRIPTION");
	private static final String LOG_LEVEL_OPTION_DESCRIPTION = Messages.getMessage("LOG_LEVEL_OPTION_DESCRIPTION");
	private static final String OUTPUT_FILE_DESCRIPTION = Messages.getMessage("OUTPUT_FILE_DESCRIPTION");
	private static final String INPUT_FILE_DESCRIPTION = Messages.getMessage("INPUT_FILE_DESCRIPTION");
	private static final String VALIDATION_STRINGENCY_OPTION_DESCRIPTION = Messages.getMessage("VALIDATION_STRINGENCY_DESCRIPTION");
	private static final String NO_HTML_DESCRIPTION = Messages.getMessage("NO_HTML_OPTION_DESCRIPTION");
	
	private final OptionParser parser = new OptionParser();
	private final OptionSet options;
//	private final List<String> fileNames;
//	private final String[] fastqs;
	private final String outputFileName;
	private final String log;
	private final Integer minBinSize;
	private final String input;
	private final String logLevel;

	@SuppressWarnings("unchecked")
	Options(final String[] args) throws Exception {
//		parser.acceptsAll(asList("h", "help"), HELP_DESCRIPTION);
//		parser.acceptsAll(asList("v", "version"), VERSION_DESCRIPTION);
		parser.accepts("help", HELP_DESCRIPTION);
		parser.accepts("version", VERSION_DESCRIPTION);
		
		parser.accepts("input", INPUT_FILE_DESCRIPTION).withRequiredArg().ofType(String.class);
		parser.accepts("output", OUTPUT_FILE_DESCRIPTION).withRequiredArg().ofType(String.class);
//		parser.accepts("exclude", EXCLUDES_OPTION_DESCRIPTION).withRequiredArg().ofType(String.class).withValuesSeparatedBy(',');
		parser.accepts("minBinSize", MIN_BIN_SIZE_OPTION_DESCRIPTION).withRequiredArg().ofType(Integer.class);
		parser.accepts("log", LOG_OPTION_DESCRIPTION).withRequiredArg().ofType(String.class);
		parser.accepts("loglevel", LOG_LEVEL_OPTION_DESCRIPTION).withRequiredArg().ofType(String.class);
		parser.accepts("fastqs", TAGS_CHAR_OPTION_DESCRIPTION).withRequiredArg().ofType(String.class).withValuesSeparatedBy(',');
		parser.posixlyCorrect(true);

		options = parser.parse(args);
		
		// log file
		log = (String) options.valueOf("log");
		// logLevel
		logLevel = (String) options.valueOf("loglevel");
		
		minBinSize =  (Integer) options.valueOf("minBinSize");
		input = (String) options.valueOf("input");
		
		// inputs
//		List<String> inputs = (List<String>) options.valuesOf("fastqs");
//		fastqs = new String[inputs.size()];
//		inputs.toArray(fastqs);
		
		
		// output
		outputFileName = (String) options.valueOf("output");
		
		if ( ! options.nonOptionArguments().isEmpty()) {}
//			throw new IllegalArgumentException(Messages.getMessage("USAGE"));
		
//		List<String> nonoptions = options.nonOptionArguments();
//		fileNames = new String[nonoptions.size()];
//		nonoptions.toArray(fileNames);
	}

	boolean hasVersionOption() {
		return options.has("version");
	}

	boolean hasHelpOption() {
		return options.has("help");
	}
	boolean hasMinBinSizeOption() {
		return options.has("minBinSize");
	}
	
	boolean hasNoHtmlOption() {
		return options.has("nohtml");
	}
	
	boolean hasLogOption() {
		return options.has("log");
	}
	
	boolean hasLogLevelOption() {
		return options.has("loglevel");
	}
	
	Integer getMinBinSize() {
		return minBinSize;
	}
	
	String getInput() {
		return input;
	}

//	String[] getFastqs() {
//		return fastqs;
//	}
//	List<String> getFileNames() {
//		return fileNames;
//	}
	
	
	String getLog() {
		return log;
	}
	String getLogLevel() {
		return logLevel;
	}
	public String getOutputFileName() {
		return outputFileName;
	}

	void displayHelp() throws Exception {
		parser.printHelpOn(System.err);
	}
	
	String getValidation() {
		if (options.has("validation")) {
			return (String) options.valueOf("validation");
		} else return null;
	}

}
