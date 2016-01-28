package au.edu.qimr.indel.pileup;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.util.IndelUtils.SVTYPE;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.picard.SAMFileReaderFactory;

public class IndelPileupTest {
//	static final String inputIndel = "indel.vcf"; 
	static final String inputBam = "tumor.sam"; 
	QLogger logger = QLoggerFactory.getLogger(IndelPileupTest.class);
	
	@BeforeClass
	public static void createInput() {	
		CreateSam();
	}
	
	@AfterClass
	public static void deleteInput() {	
		new File(inputBam).delete();
	}	
		
	@Test
	public void insertTest() throws Exception{
								
		//get pool		
		SamReader inreader =  SAMFileReaderFactory.createSAMFileReader(new File(inputBam));
		List<SAMRecord> pool = new ArrayList<SAMRecord>();
        for(SAMRecord record : inreader)
        	pool.add(record);       
        inreader.close();
		 
        //pileup
		List<VcfRecord> vcfs = new ArrayList<VcfRecord>();
		vcfs.add(new VcfRecord("chr1",	183014,	null, "G",	"GTT"));
		vcfs.add(new VcfRecord("chr1",	183014,	null, "G",	"GT"));
		IndelPosition indel = new IndelPosition (vcfs, SVTYPE.INS);
		IndelPileup  pileup = new IndelPileup( indel, 13, 3); 		
		pileup.pileup(pool);
	 		
		//assert first insertion vcf
        assertTrue(pileup.getInformativeCount() == 3);
        assertTrue(pileup.getsuportReadCount(0) == 1); 
        assertTrue(pileup.getnovelStartReadCount(0)== 1);
        assertTrue(pileup.getparticalReadCount(0) == 1); 	        
        
        //assert second insertion vcf
        assertTrue(pileup.getInformativeCount() == 3);
        assertTrue(pileup.getsuportReadCount(1) == 1); 
        assertTrue(pileup.getnovelStartReadCount(1)== 1);
        assertTrue(pileup.getparticalReadCount(1) == 1); 	
	}	
	
	@Test
	public void deleteTest() throws Exception{
 		//get delete indel
		VcfRecord vs = new VcfRecord("chr1", 197, null, "CAG", "C");				 
		IndelPosition indel = new IndelPosition (vs);
		
		//make pool
		List<SAMRecord> pool = new ArrayList<SAMRecord>();				
		try(SamReader inreader =  SAMFileReaderFactory.createSAMFileReader(new File(inputBam));){
	        for(SAMRecord record : inreader){
	        	if(record.getAlignmentStart() <= indel.getEnd())
	        	pool.add(record);
	        }
		}
        IndelPileup pileup = new IndelPileup( indel, 13, 3); 	
        pileup.pileup(pool);
                
        assertTrue(pileup.getmotif(0).equals("AG"));
        assertTrue(pileup.getInformativeCount() == 2);
        assertTrue(pileup.getsuportReadCount(0) == 2); 
        assertTrue(pileup.getnovelStartReadCount(0)== 2);
        assertTrue(pileup.getparticalReadCount(0) == 0); 		
	}

		
    public static void CreateSam(){
        List<String> data = new ArrayList<String>();
        data.add("@HD	VN:1.0	SO:coordinate");
        data.add("@RG	ID:20140717025441134	SM:eBeads_20091110_CD	DS:rl=50");
        data.add("@PG	ID:qtest::Test	VN:0.2pre");
        data.add("@SQ	SN:chr1	LN:249250621");
        data.add("@SQ	SN:chr11	LN:243199373");
        data.add("@CO	create by qcmg.qbamfilter.filter::TestFile");
        data.add("1997_1173_1256	99	chr1	183011	60	112M1D39M	=	183351	491	" + 
        		"TATGTTTTTTAGTAGAGACAGGGTCTCACTGTGTTGCCCAGGCTAGTCTCTAACTCCTGGGCTCAAATTATCCTCCCCACTTGGCCTCCCAAAAGGATTGGATTACAGGCATAAGCCACTGCCCCAAGCCTAAAATTTTTTAAAGTACCAT	FFFFFFJFJJJJFJJFJFJJFAJJJJJJJJJJFJJJJJJJFJJJJJJJFJFJJJJJJJJJFJAJJJJJJJJJJJJJJJJFJFJJJJFJJJFFFJJJJAFFJFJJJJJFJFJFJJJFJJJFFJJJAFFJJJJFAFJFFFJJAFFAJJJJFJ7	ZC:i:4	MD:Z:112^A39	PG:Z:MarkDuplicates.5	RG:Z:20140717025441134");	                                 
        data.add("1997_1173_1257	163	chr1	182999	60	16M1I105M29S	=	183397	540	" + 
        		"TACATTTAAAAATATGTTTTTTTAATAGAGACAGGGTCTCACTGTGTTGCCCAGGCTAGTCTCAAACTCCTGGGCTCAAATTATCCTCCCCACTTGGCCTCCCAAAAGGATTGGATTACAGGNANNNNNNNNNGNCNNNNNGCTAAAANTT	AAFFFJJJJJJJJJJJAJJJFJJJJJJ<J7FFJF-JFAFJJ<AJAFAJFF-FJJ-FJJAAFFJAFA-FFAFF<FFAFJAFJ<JJA7F-<-AJ<<J<F<FFJ-J<A7J-F-FJA7-<F<7J<<#-#########A#A#####AAFFFJA#--	ZC:i:4	MD:Z:23G38T58	PG:Z:MarkDuplicates.5	RG:Z:20140717025441134");
        data.add("1997_1173_1258	163	chr1	183001	9	14M2I106M29S	=	183287	440	" + 
        		"CATTTAAAAATATGTTTTTTTTAATAGAGACAGGGTCTCACTGTGTTGCCCAGGCTAGTCTCAAACTCCTGGGCTCAAATTATCCTCCCCACTTGGCCTCCCAAAAGGATGGGATTACAGGCNTNNNNNNNNNCNCNNNNNCTAAAATNTT	AFFFFJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJAJJJJJJJJJJAFJJJFJJJJJJJJJJJJJJJJJJJFJFFFAJJJJJJJJJJJJJJJJJJJJJ#A#########J#A#####JFJJJJF#JF	ZC:i:4	MD:Z:21G38T47T11	PG:Z:MarkDuplicates.5	RG:Z:20140717025441134");

       //deletion
       data.add("1997_1173_1267	113	chr1	64	60	134M2D17M	chr1	72846	0	" +
       "GAAAATACTAAACCACACCAGGTGTGGTGTCACATGCCTGTGGTCTCAGGTACTTGGGAGGCTGAGGTGGGAGGATCGCTTGAACCCAGGAAGTTGAGGCTGCAGTGAGTTGTGATTACACCAGCCTGGGTGACAGTGTCACCCTGTCTCA	JF7-<7--7-77-JJFFFJFJF<JJ<<JAFJAJJAAF<JJ<AFJ-JJAJJJAJAJJAAAJF-JFF7FFJFJAFAFAFJA<JF--FJA-F--JJAJJFJJ<FJJJJ<JJJJJJJJJJJJJJ<FAFJ<AA-JJJ<JJJJJJJJJFAFJAFFAA	ZC:i:5	MD:Z:12T121^TC17	PG:Z:MarkDuplicates.5	RG:Z:20140717025441134");
       data.add("1997_1173_1268	177	chr1	67	57	131M2D20M	chr1	72680	0	" +
       "AATACTAAAACACACCAGGTGTGGTGTCACATGCCTGTGGTCTCAGGNANTNGNGANGNTNAGGTGGGAGGATCGCTTGAACCCAGGAAGTTGAGGCTGCAGTGAGTTGTGATTACACCAGCCTGGGTGACAGTGTCACCCTGTCTCAAAA	JJJJFJJFFFJJFJJJFFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJ#J#J#J#JJ#J#F#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJFFFAA	ZC:i:4	MD:Z:47A1G1A1C2C1G1C70^TC20	PG:Z:MarkDuplicates.5	RG:Z:20140717025441134");        
        
        BufferedWriter out;
        try {
           out = new BufferedWriter(new FileWriter(inputBam ));
           for (String line : data) 
        	   out.write(line + "\n");          
           out.close();
        } catch (IOException e) {
            System.err.println("IOException caught whilst attempting to write to SAM test file: " + inputBam  + e);
        }              
    }


}
