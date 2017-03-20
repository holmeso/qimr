package au.edu.qimr.indel.pileup;

import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.After;

import org.junit.Before;

import org.junit.Test;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.model.ChrRangePosition;
import org.qcmg.common.vcf.VcfFormatFieldRecord;
import org.qcmg.common.vcf.header.*;
import org.qcmg.common.vcf.header.VcfHeaderRecord;
import org.qcmg.common.vcf.header.VcfHeaderUtils;

import au.edu.qimr.indel.Main;
import au.edu.qimr.indel.Q3IndelException;
import au.edu.qimr.indel.Support;

public class ReadIndelsTest {
	public final static String input1 = "input1.vcf";
	public final static String input2 = "input2.vcf";
	public final static String input3 = "input3.vcf";
	
	@Before 
	public void createInput() {	 
		createVcf();		
	}
	
	@After
	public void clearInput() {	 		
		new File(input1).delete();
		new File(input2).delete();
		new File(input3).delete();
	}
	
	@Test
	public void multiAltTest(){
		
		//only keep one variants since the allel CA is CTX
		//we only deal with indel which is single base on ref or allel
		//indel size can't be exceed 200
		try{
			ReadIndels indelload = new ReadIndels( new File[]{new File(input3)});		
			//indelload.LoadIndels(new File(input3), "");				
			Map<ChrRangePosition, IndelPosition> positionRecordMap = indelload.getIndelMap();
			assertTrue(positionRecordMap.size() == 1);
			 
			for( ChrPosition key : positionRecordMap.keySet()){
				IndelPosition indel = positionRecordMap.get(key);
				assertTrue(indel.getIndelVcfs().size() == 1);
				assertEquals("C", indel.getIndelVcf(0).getAlt());
				
				//check format GT field
				List<String> format = indel.getIndelVcfs().get(0).getFormatFields();
				
				VcfFormatFieldRecord record = new VcfFormatFieldRecord(format.get(0), format.get(1));
				assertTrue(record.getField("GT").equals("1/."));
				
				record = new VcfFormatFieldRecord(format.get(0), format.get(2));
				assertTrue(record.getField("GT").equals("1/1"));
			}

		}catch(Exception e){
			assertFalse(true);
		}		
	}

		
	@Test
	//test different records with same hascode, but hashCode() can change from Java versions
	public void LoadIndelsTest2()  {
		List<String> data = new ArrayList<String>();
		data.add("chr2\t71697867\t.\tTTTCC\tT\t115.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=0.361;ClippingRankSum=-0.361;DP=11;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQ0=0;MQRankSum=1.881;QD=10.52;ReadPosRankSum=-0.922;SOR=1.609\tGT:AD:DP:GQ:PL\t0/1:5,4:9:99:153,0,145");
		data.add("chr2\t241133989\t.\tGAGGTGGAGCGTAGTGTTGAAATTGCATCCATTGTGGGGCAGTGTTGGA\tG\t457.73\t.\t AC=1;AF=0.500;AN=2;BaseQRankSum=2.462;ClippingRankSum=-0.359;DP=26;FS=3.628;MLEAC=1;MLEAF=0.500;MQ=59.33;MQ0=0;MQRankSum=-2.103;QD=17.61;ReadPosRankSum=-0.449;SOR=0.991\tGT:AD:DP:GQ:PL\t0/1:13,13:26:99:495,0,2645");
			
		Support.createVcf(data, input1);
		
						
		try{
			ReadIndels indelload = new ReadIndels(new File[]{ new File(input1) } );
			//load single file
//			indelload.LoadIndels(new File(input1), "");	
			Map<ChrRangePosition, IndelPosition> positionRecordMap = indelload.getIndelMap();
			assertTrue(positionRecordMap.size() == 2);		
			
			int code1 = 0, code2 = 0;
			for( ChrPosition key : positionRecordMap.keySet()){
				IndelPosition indel = positionRecordMap.get(key);
	 
				if(indel.getStart() == 71697868 && indel.getFullChromosome().equals("chr2")){
					indel.getMotif(0).equals("TTCC");
					code1 = indel.getIndelVcf(0).hashCode();
				}else if(indel.getStart() == 241133990 && indel.getFullChromosome().equals("chr2")){
					indel.getMotif(0).equals("AGGTGGAGCGTAGTGTTGAAATTGCATCCATTGTGGGGCAGTGTTGGA");
					code2 = indel.getIndelVcf(0).hashCode();
				}else
					assertFalse(true);
				//can't run it for all environment 
//				assertTrue(code1 == code2);			 
 			}
	
		}catch(Exception e){
			System.out.println(Q3IndelException.getStrackTrace(e));
			assertFalse(true);
		}
		
	}		
			
	@Test
	//load inputs only in case of pindel
	public void LoadIndelsTest()  {
		
		try{			
			//load single file
			ReadIndels indelload = new ReadIndels( new File[]{ new File(input1) });	//indelload.LoadIndels(new File(input1),"");	
			assertTrue(getHeaderLineCounts(indelload.getVcfHeader()) == 7);
						
			//load two files eg pindel
			File[] fs = new File[]{new File(input1), new File(input2)};
			indelload = new ReadIndels( fs);	
			assertTrue(getHeaderLineCounts(indelload.getVcfHeader()) == 8);	
			
			Map<ChrRangePosition, IndelPosition> positionRecordMap = indelload.getIndelMap();
			assertTrue(positionRecordMap.size() == 3);		
			for( ChrPosition key : positionRecordMap.keySet()){
				IndelPosition indel = positionRecordMap.get(key);
				assertFalse(indel.getIndelVcf(0).getFormatFields().get(1).equals(indel.getIndelVcf(0).getFormatFields().get(2)));
				if(indel.getStart() == 59033423){					
					assertTrue( indel.getMotif(0).equals("C"));
					assertTrue( indel.getMotif(1).equals("CG"));
					 					
					//check GT:GD  from existing one,  no long overwrite existing one
					assertTrue(indel.getIndelVcf(0).getFormatFields().get(0).equals("GT:GD:AD:DP:GQ:PL"));
					assertTrue(indel.getIndelVcf(0).getFormatFields().get(1).equals("0/1:T/TC:7,4:11:99:257,0,348"));
					assertTrue(indel.getIndelVcf(0).getFormatFields().get(2).equals("0/1:T/TC:17,2:19:72:72,0,702"));	
					
					assertTrue(indel.getIndelVcf(1).getFormatFields().get(1).equals("0/.:T/A:7,5"));
					assertTrue(indel.getIndelVcf(1).getFormatFields().get(2).equals("./.:A/TC:9,9"));		
					
				}			
	 			 			 		
				//change inputs order
				fs = new File[]{new File(input2), new File(input1)};
				assertTrue(getHeaderLineCounts(indelload.getVcfHeader()) == 8);	
				
//				indelload = new ReadIndels( null ,fs  );	
//				indelload.LoadIndels(new File(input2),"");			
//				assertTrue(getHeaderLineCounts(indelload.getVcfHeader()) == 7);				
//				//load second file, in case of pindel
//				indelload.LoadIndels(new File(input1),"");
				
			}
		}catch(Exception e){
			System.out.println(Q3IndelException.getStrackTrace(e));
			assertFalse(true);
		}
	}
	
	static int getHeaderLineCounts(VcfHeader header){
		int no = 0; 
		for(final VcfHeaderRecord record: header )  
			no ++;		
		return no; 		
	}
	
	public void createVcf(){	
		
		List<String> head = new ArrayList<String>();
		head.add("##fileformat=VCFv4.1");
		head.add("##FILTER=<ID=PASS,Description=\"test\">");
		head.add("##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"test\">");
		head.add("##contig=<ID=chrMT,length=16569>");     
		head.add("##contig=<ID=chrY,length=59373566>");	  
				    
		List<String> head1 = new ArrayList<String>(head);
		head1.add("##INFO=<ID=SOMATIC1,Number=0,Type=Flag,Description=\"test1\">");       
		head1.add(VcfHeaderUtils.STANDARD_FINAL_HEADER_LINE + "\tFORMAT\tS1\tS2");
		List<String> data1 = new ArrayList<String>();
		data1.add("chrY	59033286	.	GT	G	724.73	PASS	SOMATIC1	GT:AD:DP:GQ:PL	0/0:0:0:0:0,0,0	0/1:80,17:97:99:368,0,3028");
		data1.add("chrY	59033423	.	T	TC	219.73	PASS	SOMATIC1	GT:AD:DP:GQ:PL	0/1:7,4:11:99:257,0,348	0/1:17,2:19:72:72,0,702"); 
		Support.createVcf(head1, data1, input1);
		       
		head1 = new ArrayList<String>(head);
		head1.add("##PG=\"creating second file\"");
		head1.add(VcfHeaderUtils.STANDARD_FINAL_HEADER_LINE + "\tFORMAT\tS1\tS2");      
		data1.clear();        
		data1.add("chrY	59033285	.	GGT	G	724.73	PASS	SOMATIC	GT:AD:DP:GQ:PL	0/1:131,31:162:99:762,0,4864	0/1:80,17:97:99:368,0,3028");
		data1.add("chrY	59033286	.	GT	G	724.73	PASS	SOMATIC	GT:AD:DP:GQ:PL	0/1:131,31:162:99:762,0,4864	0/1:80,17:97:99:368,0,3028");
		data1.add("chrY	59033423	.	T	A,TC,TCG	219.73	PASS	SOMATIC	GT:AD	0/1:7,5	1/2:9,9");            
		Support.createVcf(head1, data1, input2);
		
		data1.clear();
		data1.add("chrY	59033286	.	CAA	C,CA	724.73	PASS	.	GT:AD:DP:GQ:PL	1/2:14,38,25:77:99:1229,323,592,448,0,527	1/1:14,38,25:77:99:1229,323,592,448,0,527");
		data1.add("chrY	59033287	.	GTGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	G	724.73	PASS	.	GT:AD	1/2:14,38,25	1/1:14,38,25");
		data1.add("chrY	59033287	.	G	GTGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	724.73	PASS	.	GT:AD	1/2:14,38,25	1/1:14,38,25");

		Support.createVcf(head1, data1, input3);
            
	}
	
	
	
}
