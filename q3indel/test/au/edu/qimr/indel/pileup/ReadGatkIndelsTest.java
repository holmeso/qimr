package au.edu.qimr.indel.pileup;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.model.ChrRangePosition;

import au.edu.qimr.indel.Q3IndelException;
import au.edu.qimr.indel.Support;

public class ReadGatkIndelsTest {
	final String input1 = "input1.vcf";
	final String input2 = "input2.vcf";
//	final String input3 = "input3.vcf";
	@Before 
	public void createInput() {	createVcf(); }
	
	@After
	public void clearInput() {	 		
		new File(input1).delete();
		new File(input2).delete();
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
		head1.add("#CHROM	POS	ID      REF     ALT     QUAL	FILTER	INFO	FORMAT	S1	S2"); 		
		List<String> data1 = new ArrayList<String>();
		data1.add("chrY	59033286	.	GT	G	724.73	PASS	SOMATIC1	GT:AD:DP:GQ:PL	0/1:80,17:97:99:368,0,3028");
		data1.add("chrY	59033423	.	T	TC	219.73	PASS	SOMATIC1	GT:AD:DP:GQ:PL	0/1:7,4:11:99:257,0,348"); 
		Support.createVcf(head1, data1, input1);
		       
		head1 = new ArrayList<String>(head);
		head1.add("#CHROM	POS	ID      REF     ALT     QUAL	FILTER	INFO	FORMAT	S1	S2");        
		data1.clear();        
		data1.add("chrY	59033285	.	GGT	G	724.73	PASS	SOMATIC	GT:AD:DP:GQ:PL	0/1:131,31:162:99:762,0,4864");
		data1.add("chrY	59033286	.	GT	G	724.73	PASS	SOMATIC	GT:AD:DP:GQ:PL	1/1:131,31:162:99:762,0,4864");
		data1.add("chrY	59033423	.	T	A,TC,TCG	219.73	PASS	SOMATIC	GT:AD	2/3:7,5	1/3:9,9");            
		Support.createVcf(head1, data1, input2);
	}
	
	
	@Test
	public void mergeIndelsTest(){
	 
		try{
			File[] fs = new File[]{new File(input1), new File(input2)};
			ReadGatkIndels indelload = new ReadGatkIndels( fs);	
			assertTrue(ReadIndelsTest.getHeaderLineCounts(indelload.getVcfHeader()) == 7);				
						
			Map<ChrRangePosition, IndelPosition> positionRecordMap = indelload.getIndelMap();
			assertTrue(positionRecordMap.size() == 3);		
			for( ChrPosition key : positionRecordMap.keySet()){
				IndelPosition indel = positionRecordMap.get(key);
				if(indel.getStart() == 59033286){	//	 chrY	59033285	.	GGT	G

					assertTrue( indel.getIndelVcf(0).getFormatFields().get(1).equals(".:.:.:.:.:.") );
					assertTrue( indel.getIndelVcf(0).getFormatFields().get(2).equals("0/1:GGT/G:131,31:162:99:762,0,4864") );	
					assertTrue( indel.getIndelVcf(0).getInfo().equals("SOMATIC" ) ); //info column from second file  
				}else if(indel.getStart() == 59033287){
					//indels only appear on second file,  
					assertTrue(indel.getIndelVcf(0).getFormatFields().get(1).equals("0/1:GT/G:80,17:97:99:368,0,3028" )); 
					assertTrue(indel.getIndelVcf(0).getFormatFields().get(2).equals("1/1:G/G:131,31:162:99:762,0,4864" ));	
					assertTrue(indel.getIndelVcf(0).getInfo().equals("SOMATIC1" )); //info column from first file  
				}else if(indel.getStart() == 59033423){	
					//merge indels but split alleles
					assertTrue(indel.getIndelVcf(0).getFormatFields().get(1).equals("0/1:T/TC:7,4:11:99:257,0,348"));
 					assertTrue(indel.getIndelVcf(0).getFormatFields().get(2).equals("1/.:TC/TCG:.:.:.:."));
 					assertTrue(indel.getIndelVcf(0).getAlt().equals("TC"));
 					assertTrue(indel.getIndelVcf(0).getInfo().equals("SOMATIC1" )); //info column from first file  
 					
					assertTrue(indel.getIndelVcf(1).getFormatFields().get(1).equals(".:.:.:.:.:." ));
					assertTrue(indel.getIndelVcf(1).getFormatFields().get(2).equals("./1:TC/TCG:.:.:.:."));
					assertTrue(indel.getIndelVcf(1).getAlt().equals("TCG"));
					assertTrue(indel.getIndelVcf(1).getInfo().equals("SOMATIC" ));  //info column from second file
				}									
			}
		}catch(Exception e){
			System.err.println(Q3IndelException.getStrackTrace(e));
			assertFalse(true);
		}
		
	}	
	
	
	
}
