package au.edu.qimr.indel.pileup;


import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.AbstractQueue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMRecord;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.meta.QBamId;
import org.qcmg.common.meta.QExec;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.string.StringUtils;
import org.qcmg.common.util.IndelUtils;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.picard.SAMFileReaderFactory;
import org.qcmg.picard.util.QBamIdFactory;
import org.qcmg.qbamfilter.query.QueryExecutor;
import org.qcmg.vcf.VCFFileWriter;

import au.edu.qimr.indel.Options;
import au.edu.qimr.indel.pileup.Homopolymer;
import au.edu.qimr.indel.pileup.IndelPileup;
import au.edu.qimr.indel.pileup.IndelPosition;
import au.edu.qimr.indel.pileup.ReadIndels;


public class IndelMT {
	public static final int  MAXRAMREADS = 1500; //maximum number of total reads in RAM
	
 
	class contigPileup implements Runnable {

		private final AbstractQueue<IndelPosition> qIn;
		private final AbstractQueue<IndelPileup> qOut;
		private final Thread mainThread;
		final CountDownLatch pLatch;
		private SAMSequenceRecord contig;
		private File bam; 
		private QueryExecutor exec;
		
		//unit Test only
		contigPileup(){
			this.qIn = null;
			this.qOut = null;
			this.mainThread = null;
			this.pLatch = new CountDownLatch(2); // testing only
		}

		/**
		 * 
		 * @param qIn : store SAM record from input file
		 * @param qOutGood : store unmatched record based on query
		 * @param qOutBad: store unmatched record based on query (null is allowed)
		 * @param query : query string
		 * @param maxRecords : queue size
		 * @param mainThread : parent thread
		 * @param rLatch : the counter for reading thread
		 * @param fLatch : the counter for filtering thread (current type)
		 */
		contigPileup(SAMSequenceRecord contig,  AbstractQueue<IndelPosition> qIn, File bam, QueryExecutor exec,
				AbstractQueue<IndelPileup> qOut, Thread mainThread, CountDownLatch latch)  {
			this.qIn = qIn;
			this.qOut = qOut;
			this.mainThread = mainThread;
			this.pLatch = latch;
//			this.wLatch = wLatch;
			this.bam = bam;
			this.contig = contig;
			this.exec = exec; 
 		}

		@Override
		public void run() {
		 	List<SAMRecord> current_pool = new ArrayList<SAMRecord>();
		 	List<SAMRecord> next_pool = new ArrayList<SAMRecord>(); 
		 	int size = qIn.size();	
			if (size <= 0) {
		 		logger.info("There is no indel fallen in contig: " + contig.getSequenceName() );		 		
		 		return;
		 	}
			
			IndelPosition topPos= qIn.poll();
			File index = new File(bam.getAbsolutePath() + ".bai");
			
//			try (SAMFileReader Breader = new SAMFileReader(bam )){	
			try (SamReader Breader =  SAMFileReaderFactory.createSAMFileReader(bam, index); ){		
				SAMRecordIterator ite = Breader.query(contig.getSequenceName(), 0, contig.getSequenceLength(),false);		
			 	while (ite.hasNext()) {	
			 		SAMRecord re = ite.next(); 
			 		//bam file already sorted, skip non-indel region record
			 		if(re.getAlignmentEnd() < topPos.getStart()) continue; 
			 		
			 		//only interested pass filter record
			 		boolean passFilter; 	
			 		if(exec != null )
						passFilter = exec.Execute(re);
					else
						passFilter = !re.getReadUnmappedFlag() && (!re.getDuplicateReadFlag() || !options.excludeDuplicates());
			 		if(! passFilter ) continue; 			 					 	 
			 			 		
			 		//whether in current indel region
			 		if(re.getAlignmentStart() <= topPos.getEnd() && current_pool.size() < MAXRAMREADS )
			 			current_pool.add(re);
			 		else{
			 			next_pool.add(re); 
			 			//pileup
			 			IndelPileup pileup= new IndelPileup(topPos, options.getSoftClipWindow(), options.getNearbyIndelWindow());
			 			pileup.pileup(current_pool);
			 			qOut.add(pileup);
			 			
			 			//prepare for next indel position
			 		//	topPos = qIn.poll();
			 			if( (topPos = qIn.poll()) == null) break; 
			 			
			 			resetPool(topPos,  current_pool, next_pool); 	 
//			 			//debug
//			 			if(current_pool.size() > 1000 ||next_pool.size()>1000 )
//			 			System.out.println("debug: " + topPos.getChrPosition().toIGVString() + " : " + current_pool.size() + " (current, next) " +next_pool.size());
			 			
			 		}
			 	}	 	
			
			 	//after loop check all pool
			 	do{			
			 		//check whether previous loop also used up all indel position
			 		if(topPos == null) break; 
			 			
		 			IndelPileup pileup= new IndelPileup(topPos, options.getSoftClipWindow(), options.getNearbyIndelWindow());
		 			pileup.pileup(current_pool);
		 			qOut.add(pileup);
					
					if( (topPos = qIn.poll()) == null) break; 

					resetPool(topPos,  current_pool, next_pool); 							
//		 			//debug
//					if(current_pool.size() > 1000 ||next_pool.size()>1000 )
//		 			System.out.println(topPos.getChrPosition().toIGVString() + " (final): " + current_pool.size() + " (current, next) " +next_pool.size());
			 	}while( true ) ;					 					 
			} catch (Exception e) {
				logger.error("Exception caught in pileup thread", e);
				mainThread.interrupt();
			} finally {
				pLatch.countDown();
				logger.info( size + " indels is completed pileup from " + contig.getSequenceName() + " on " + bam.getName());

 			}			
		}
	
		
		/**
		 * it swap SAMRecord between currentPool and nextPool. After then, the currentPool will contain all SAMRecord overlapping topPos position, 
		 * the nextPool will contain all SAMRecord start after topPos position.  All SAMRecord end before topPos position will be remvoved from both pool. 
		 * @param topPos:   pileup position
		 * @param currentPool: a list of SAMRecord overlapped previous pileup Position
		 * @param nextPool: a list of SAMRecord behind previous pileup Position
		 */
		 void resetPool( IndelPosition topPos, List<SAMRecord> currentPool, List<SAMRecord> nextPool){
			
			
				List<SAMRecord> tmp_current_pool = new ArrayList<SAMRecord>();							
				List<SAMRecord> tmp_pool = new ArrayList<SAMRecord>();	
				
//				if(nextPool.size() > MAXRAMREADS ) nextPool.clear();				
//				if(currentPool.size() >  MAXRAMREADS ) currentPool.clear();
//				//debug
//				if(currentPool.size() > MAXRAMREADS || nextPool.size() > MAXRAMREADS ){
//					System.out.println(topPos.getChrPosition().toIGVString() + ": bf moving, next pool size : " + nextPool.size() + ", current pool size: " + currentPool.size());
//				}				
				tmp_pool.addAll(nextPool);
				
				//check read record behind on current position			
				for( SAMRecord  re : tmp_pool ){
					//aligned position before indel
					if(re.getAlignmentEnd() < topPos.getStart())
						nextPool.remove(re);
					//aligned position cross indel
					else if(re.getAlignmentStart() <= topPos.getEnd()){ 	 					 
						tmp_current_pool.add(re);	
						nextPool.remove(re);
					}			 
				}	
				
				tmp_pool.clear();
				tmp_pool.addAll(currentPool);
				//check already read record  for previous pileup
				for( SAMRecord  re1 : tmp_pool ){
					//aligned position before indel
					if(re1.getAlignmentEnd() < topPos.getStart())
						currentPool.remove(re1);
					//aligned position after indel
					else if(re1.getAlignmentStart() > topPos.getEnd()){
						nextPool.add(re1);
						currentPool.remove(re1);
					}
				}
				//debug
				//if(topPos.getChrPosition().getPosition() ==148848273){
//				if(currentPool.size() > MAXRAMREADS || nextPool.size() > MAXRAMREADS ){
//					System.out.println("current pool reduced to " + currentPool.size());
//					System.out.println("then adding from next pool : " + tmp_current_pool.size());
//				}	
				
				//merge samrecord
				currentPool.addAll(tmp_current_pool);
		}		
		
	}
		
	
	class homopoPileup implements Runnable {
		private final AbstractQueue<IndelPosition> qIn;
		private final AbstractQueue<Homopolymer> qOut;
		private final Thread mainThread;
		private File referenceFile; 
		private String contig; 
		private final byte[] referenceBase;
		private int window; 
		final CountDownLatch pLatch;
		
		homopoPileup(String contig,   AbstractQueue<IndelPosition> qIn, File reference,  
				AbstractQueue<Homopolymer> qOut, int window, Thread mainThread, CountDownLatch latch) {
			this.qIn = qIn;
			this.qOut = qOut;
			this.mainThread = mainThread;
			this.referenceFile = reference; 
			this.window = window;
			this.pLatch = latch; 
			this.contig = contig; 
			
			IndexedFastaSequenceFile indexedFasta = Homopolymer.getIndexedFastaFile(referenceFile);
			this.referenceBase = indexedFasta.getSequence(contig).getBases();
 	
 		}
		 
 
		@Override
		public void run() {
			int size = qIn.size();
			try {
  
				try {
					IndelPosition pos;					
					while ((pos = qIn.poll()) != null)  
						qOut.add(new Homopolymer(pos, referenceBase, window));
				} finally {					 
					pLatch.countDown();
				}
			}catch(Exception e){
				logger.error("Exception caught in homopolymer thread", e);
				mainThread.interrupt();

			} finally {
				pLatch.countDown();
				logger.info(size  + " indels had been checked homopolymer from " + contig);
			}
		}
		
	}
	
	Options options; 
	QLogger logger; 
	ReadIndels indelload;
		
	private final List<SAMSequenceRecord> sortedContigs = new ArrayList<SAMSequenceRecord>();
	private Map<ChrPosition, IndelPosition> positionRecordMap ;
	
	//unit test purpose
	@Deprecated
	IndelMT(){}
	public IndelMT(Options options, QLogger logger) throws IOException  {		
		this.options = options;	
		this.logger = logger; 
		
		SamReader reader = SAMFileReaderFactory.createSAMFileReader(options.getTestBam());  
		for (final SAMSequenceRecord contig : reader.getFileHeader().getSequenceDictionary().getSequences())  
			sortedContigs.add(contig);
		reader.close(); 
		
		//loading indels 
		this.indelload = new ReadIndels(logger);		
		if(options.getRunMode().equalsIgnoreCase("gatk")){	
			//first load control			
			indelload.LoadIndels(options.getControlInputVcf());		
			//then test
			indelload.appendIndels(options.getTestInputVcf());
		}else if(options.getRunMode().equalsIgnoreCase("pindel")){	
			for(int i = 0; i < options.getInputVcfs().size(); i ++)
				indelload.LoadIndels(options.getInputVcfs().get(i));		
		}
	}
	
	
//	public IndelMT(File inputVcf, Options options, QLogger logger) throws IOException  {		
//		this.options = options;	
//		this.logger = logger; 	
//	}
//
//	public IndelMT(File inputTumourVcf, File inputNormalVcf, Options options, QLogger logger) throws Exception {
//		this(inputTumourVcf, options, logger); 
//		indelload.appendIndels(inputNormalVcf); 		
//	}
//	
	public int process(final int threadNo) throws Exception {
		positionRecordMap = indelload.getIndelMap();
		if(positionRecordMap == null || positionRecordMap.size() == 0){
			logger.info("Exit program since there is no indels loaded from inputs");
			return 0; 
		}			
		
        final CountDownLatch pileupLatch = new CountDownLatch(sortedContigs.size() * 2); // filtering thread               
        
        final AbstractQueue<IndelPileup> tumourQueue = new ConcurrentLinkedQueue<IndelPileup>();
        final AbstractQueue<Homopolymer> homopoQueue = new ConcurrentLinkedQueue<Homopolymer>();
        final AbstractQueue<IndelPileup> normalQueue = new ConcurrentLinkedQueue<IndelPileup>();
        // set up executor services
        ExecutorService pileupThreads = Executors.newFixedThreadPool(threadNo);    	
    	
    	//each time only throw threadNo thread, the loop finish untill the last threadNo                    	
    	for(SAMSequenceRecord contig : sortedContigs ){      		
    		pileupThreads.execute(new contigPileup(contig, getIndelList(contig), options.getTestBam() ,null ,
    				 tumourQueue, Thread.currentThread() ,pileupLatch));
    		       		
    		pileupThreads.execute(new contigPileup(contig, getIndelList(contig), options.getControlBam(),null ,
    				normalQueue, Thread.currentThread(),pileupLatch ));
    		
    		pileupThreads.execute(new homopoPileup(contig.getSequenceName(), getIndelList(contig), options.getReference(),
    				homopoQueue, options.nearbyHomopolymer, Thread.currentThread(),pileupLatch));    		
    	}
    	pileupThreads.shutdown();
    	
		// wait for threads to complete
		try {
			logger.info("waiting for  threads to finish (max wait will be 20 hours)");
			pileupThreads.awaitTermination(20, TimeUnit.HOURS);
			logger.info("All threads finished");
			
			 writeVCF( tumourQueue, normalQueue, homopoQueue,options.getOutput(),indelload.getVcfHeader());			
			
		} catch (Exception e) {
			logger.error("Exception caught whilst waiting for threads to finish: " + e.getMessage(), e);
			throw e;
		} finally {
            pileupThreads.shutdownNow();
		}
        
		return 0; 
	}
	

	private void writeVCF(AbstractQueue<IndelPileup> tumourQueue, AbstractQueue<IndelPileup> normalQueue, AbstractQueue<Homopolymer> homopoQueue, File output, VcfHeader header ) throws Exception{
		
		IndelPileup pileup;
		if(positionRecordMap == null ){
			logger.warn("the indel map: positionRecordMap point to nothing");
			return; 		
		}
			
		while((pileup = tumourQueue.poll()) != null ){
			ChrPosition pos = pileup.getChrPosition();
			IndelPosition indel = positionRecordMap.get(pos);
			indel.setPileup(true, pileup);			
		}
		while((pileup = normalQueue.poll()) != null ){
			ChrPosition pos = pileup.getChrPosition();
			IndelPosition indel = positionRecordMap.get(pos);
			indel.setPileup(false, pileup);			
		}
		
		Homopolymer homopo;
		while((homopo = homopoQueue.poll()) != null ){
			ChrPosition pos = homopo.getChrPosition();
			IndelPosition indel = positionRecordMap.get(pos);
			indel.setHomopolymer(homopo);
		}
		
		final AbstractQueue<IndelPosition> orderedList = getIndelList(null);
		logger.info("reading indel position:  " + orderedList.size() );
		try(VCFFileWriter writer = new VCFFileWriter( output)) {	
			
			List<File> inputs = new ArrayList<File>();
			if(options.getRunMode().equalsIgnoreCase("gatk")){
				inputs.add(options.getTestInputVcf());
				inputs.add(options.getControlInputVcf());				
			}else
				inputs = options.getInputVcfs();
							
			//reheader
			getHeaderForIndel(header);	
		      //debug
		    for(final VcfHeader.Record record: header)  
		    	System.out.println( "header: " + record.toString());
		 	
        	for(final VcfHeader.Record record: header)  
        		writer.addHeader(record.toString());
			 
        	//adding indels
			long count = 0;
			IndelPosition indel; 
			while( (indel = orderedList.poll()) != null)
				for(int i = 0; i < indel.getMotifs().size(); i++){					
					writer.add( indel.getPileupedVcf(i, options.getMinGematicNovelStart(), options.getMinGematicSupportOfInformative()) );	
					count ++;
				}
						
			logger.info("outputed VCF record:  " + count);	
		}
		
	}
	
	 	private void getHeaderForIndel(VcfHeader header ) throws Exception{

//		VcfHeader header = existHeader; // new VcfHeader();
		 QExec qexec = options.getQExec();
		 
		final DateFormat df = new SimpleDateFormat("yyyyMMdd");
 		header.parseHeaderLine(VcfHeaderUtils.CURRENT_FILE_VERSION);		
		header.parseHeaderLine(VcfHeaderUtils.STANDARD_FILE_DATE + "=" + df.format(Calendar.getInstance().getTime()));		
		header.parseHeaderLine(VcfHeaderUtils.STANDARD_UUID_LINE + "=" + qexec.getUuid().getValue());
		header.parseHeaderLine(VcfHeaderUtils.STANDARD_SOURCE_LINE + "=" + qexec.getToolName().getValue() + " v" + qexec.getToolVersion().getValue());
		
		header.parseHeaderLine(VcfHeaderUtils.STANDARD_DONOR_ID + "=" + options.getDonorId());
		header.parseHeaderLine(VcfHeaderUtils.STANDARD_CONTROL_SAMPLE + "=" + options.getControlSample());		
		header.parseHeaderLine(VcfHeaderUtils.STANDARD_TEST_SAMPLE + "=" + options.getTestSample());		
		
//		List<File> inputs = new ArrayList<File>();
		if(options.getRunMode().equalsIgnoreCase("gatk")){
			header.parseHeaderLine(VcfHeaderUtils.STANDARD_INPUT_LINE + "_GATK_TEST=" + options.getTestInputVcf().getAbsolutePath());
			header.parseHeaderLine(VcfHeaderUtils.STANDARD_INPUT_LINE + "_GATK_CONTROL=" + options.getControlInputVcf().getAbsolutePath());
 		}else
	 		for(int i = 0; i < options.getInputVcfs().size(); i ++)
	 			header.parseHeaderLine(VcfHeaderUtils.STANDARD_INPUT_LINE + "_PINDEL=" + options.getInputVcfs().get(i).getAbsolutePath());

		String normalBamName = options.getControlBam().getAbsolutePath();			
		header.parseHeaderLine( VcfHeaderUtils.STANDARD_CONTROL_BAM  + "=" + normalBamName);
		header.parseHeaderLine( "##qControlBamUUID=" + QBamIdFactory.getBamId(normalBamName));
		String tumourBamName = options.getTestBam().getAbsolutePath();
		header.parseHeaderLine( VcfHeaderUtils.STANDARD_CONTROL_BAM  + "=" + tumourBamName);
		header.parseHeaderLine( "##qControlBamUUID=" + QBamIdFactory.getBamId(tumourBamName));		 			 			 	
		header.parseHeaderLine( "##qAnalysisId=" + options.getAnalysisId() );
				
		//add filter
        header.addFilterLine(IndelUtils.FILTER_COVN12, IndelUtils.DESCRITPION_FILTER_COVN12 );
        header.addFilterLine(IndelUtils.FILTER_COVN8,  IndelUtils.DESCRITPION_FILTER_COVN8 );
        header.addFilterLine(IndelUtils.FILTER_COVT,  IndelUtils.DESCRITPION_FILTER_COVT );
        header.addFilterLine(IndelUtils.FILTER_HCOVN,  IndelUtils.DESCRITPION_FILTER_HCOVN );
        header.addFilterLine(IndelUtils.FILTER_HCOVT,  IndelUtils.DESCRITPION_FILTER_HCOVT );
        header.addFilterLine(IndelUtils.FILTER_MIN,  IndelUtils.DESCRITPION_FILTER_MIN );
        header.addFilterLine(IndelUtils.FILTER_NNS,  IndelUtils.DESCRITPION_FILTER_NNS );
        header.addFilterLine(IndelUtils.FILTER_TPART,  IndelUtils.DESCRITPION_FILTER_TPART );
        header.addFilterLine(IndelUtils.FILTER_NPART,  IndelUtils.DESCRITPION_FILTER_NPART );
        header.addFilterLine(IndelUtils.FILTER_TBIAS,  IndelUtils.DESCRITPION_FILTER_TBIAS );
        header.addFilterLine(IndelUtils.FILTER_NBIAS,  IndelUtils.DESCRITPION_FILTER_NBIAS );
        
		final String SOMATIC_DESCRIPTION = String.format("There are more than %d novel starts  or "
				+ "more than %.2f soi (number of supporting informative reads /number of informative reads) on control BAM",
				options.getMinGematicNovelStart(), options.getMinGematicSupportOfInformative());

		header.addInfoLine(VcfHeaderUtils.INFO_SOMATIC, "1", "String", SOMATIC_DESCRIPTION);
		header.addInfoLine(IndelUtils.INFO_NIOC, "1", "String", IndelUtils.DESCRITPION_INFO_NIOC);
		header.addInfoLine(IndelUtils.INFO_HOMCNTXT, "1", "String", IndelUtils.DESCRITPION_INFO_HOMCNTXT); 					  			
		header.addFormatLine(IndelUtils.INFO_ACINDEL, "1", "String", IndelUtils.DESCRITPION_INFO_ACINDEL);

		VcfHeaderUtils.addQPGLineToHeader(header, qexec.getToolName().getValue(), qexec.getToolVersion().getValue(), qexec.getCommandLine().getValue() 
				+  " [runMode: " + options.getRunMode() + "]");        
        		
		//header will automatic append CHROM line if get from input file	
		VcfHeaderUtils.addSampleId(header, options.getControlSample(), 1);
		VcfHeaderUtils.addSampleId(header, options.getTestSample(), 2);
		 
	}
	 
 

	/**
	 * 
	 * @param contig: contig name or null for whole reference
	 * @return a sorted list of IndelPotion on this contig; return whole reference indels if contig is null
	 */
	private  AbstractQueue<IndelPosition>  getIndelList( SAMSequenceRecord contig ){	  
		if (positionRecordMap == null || positionRecordMap.size() == 0)
			return new ConcurrentLinkedQueue<IndelPosition>(); 			  
		  
		List<IndelPosition> list = new ArrayList<IndelPosition> ();	
		if(contig == null){ //get whole reference
			list.addAll(positionRecordMap.values());	
		}else{	  //get all chrPosition on specified contig	
		for(ChrPosition pos : positionRecordMap.keySet())
			if(pos.getChromosome().equals(contig.getSequenceName()))
				list.add(positionRecordMap.get(pos));	 
		}
  
		final Comparator<String> chrComparator = new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				return sortedContigs.indexOf(o1) - sortedContigs.indexOf(o2);
			}
		};
		Collections.sort(list, new Comparator<IndelPosition>() {
			@Override
			public int compare(IndelPosition o1, IndelPosition o2) {
				final int diff = chrComparator.compare(o1.getChrPosition().getChromosome(), o2.getChrPosition().getChromosome());
				if (diff != 0) return diff;
				return o1.getChrPosition().getPosition() - o2.getChrPosition().getPosition();
			}
		});
		
		return new ConcurrentLinkedQueue<IndelPosition>(list);
  }
	 
		
}