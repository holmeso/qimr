package au.edu.qimr.utility;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import org.junit.*;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.header.VcfHeader;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
import org.qcmg.vcf.VCFFileReader;
 
public class VcfCompareTest {
	public static String outputName = "output.vcf";
	public static String inputPrimaryName = "primary.vcf";
	public static String inputAdditionalName = "additional.vcf";
	public static String args[] = "--primaryInput primary.vcf --additionalInput additional.vcf --output merge.vcf --log output.log".split(" ");			
	
	@BeforeClass
	public static void createInput() throws Exception{	
		createVcf(inputPrimaryName);
		createVcf(inputAdditionalName);
		
        final List<String> data = new ArrayList<String>();
        data.add("chr1\t250065\t.\tG\tA\t309.32\tCOVN12\tSOMATIC\tGT:AD:NNS\t.:.:0\t0/1:3,10:8");
        try(BufferedWriter out = new BufferedWriter(new FileWriter(inputPrimaryName,true));) {          
            for (final String line : data)   out.write(line +"\n");                  
         }  
		
        data.clear();
        data.add("chr1\t260065\t.\tG\tT\t.\t.\t.\tGT:AD:NNS\t.:.:0\t0/1:3,10:8");
        try(BufferedWriter out = new BufferedWriter(new FileWriter(inputAdditionalName,true));) {          
            for (final String line : data)   out.write(line +"\n");                  
         } 
	}		
 
	 @AfterClass
	 public static void deleteIO(){ 
		 new File(inputPrimaryName).delete();
		 new File(inputAdditionalName).delete();		 
		 new File(outputName).delete();		 
	 }	 

		@Test
		public void checkHeader() throws Exception  {
			String inputAdditionalName = "additional1.vcf";
	        final List<String> data = new ArrayList<String>();
	        data.add("##fileformat=VCFv4.0");
	        data.add("##fileDate=20150819");
	        data.add(VcfHeaderUtils.STANDARD_UUID_LINE + "=abcd_12345678_xzy_999666333");	        
	        data.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXTERN-MELA-20140526-102");    
	        
	        try(BufferedWriter out = new BufferedWriter(new FileWriter(inputAdditionalName));) {          
	            for (final String line : data)   out.write(line +"\n");                  
	         } 			
			
			try {
				VcfCompare compare = new VcfCompare( new File(inputPrimaryName), new File(inputAdditionalName));				
			    fail( "My method didn't throw when I expected it to" );
			} catch (RuntimeException expectedException) {
				new File(inputAdditionalName).delete();	
			}
		}	
 
	@Test
	public void checkMerge() throws Exception {
 	   //call VcfCompare		
		VcfCompare compare = new VcfCompare( new File(inputPrimaryName), new File(inputAdditionalName));				
		compare.reheader( Messages.reconstructCommandLine(args ), new Options( args)	);
		compare.VcfMerge( new File(outputName));    
			
		try ( VCFFileReader reader = new VCFFileReader(new File(outputName)) ) {  
			VcfHeader header = reader.getHeader();		
 			assertTrue(header.getFileVersion().toString().equals(VcfHeaderUtils.STANDARD_FILE_VERSION + "=VCFv4.0"));		
			String fileDate = new SimpleDateFormat("yyyyMMdd").format(Calendar.getInstance().getTime());
			assertTrue(header.getFileDate().toString().equals(VcfHeaderUtils.STANDARD_FILE_DATE + "=" + fileDate));			
 			assertFalse(header.getUUID().toString().equals(VcfHeaderUtils.STANDARD_UUID_LINE + "=abcd_12345678_xzy_999666333"));
			assertTrue(header.getInfoRecords().containsKey(Options.Info_From));
		}
		
		//check counts
		assertTrue( compare.getCountBoth() == 3 );
		assertTrue( compare.getCountAddtionalOnly() == 1 );
		assertTrue( compare.getCountPrimaryOnly() == 1 );
		assertTrue( compare.getCountOutput() == 5 );
		
		//check Record
       try ( VCFFileReader reader = new VCFFileReader(new File(outputName)) ) {        		        	
         	//load first file
        	for (final VcfRecord re : reader){   
        		int pos = re.getPosition();	        	 
        		switch (pos){
	        		case 248845:		        			 
	        		case 248945:
	        			re.getInfo().equals(Options.Info_From +"=0");
	        		case 249065:	
	        			re.getInfoRecord().getField(Options.Info_From).equals("0");
	        			break;
	        		case 250065:	
	        			assertTrue(re.getInfoRecord().getField(Options.Info_From).equals("1"));
	        			assertTrue(re.getInfoRecord().getField("SOMATIC") != null);
	        			break;
	        		case 260065:	
	        			re.getInfo().equals(Options.Info_From +"=2");
	        			break;
	        		default:
	        			throw new IllegalArgumentException("unexpected record: " + re.toString());
        		}
        	}
       }
 	}
	

	public static void createVcf(String file) throws IOException{
		final List<String> data = new ArrayList<String>();
		data.add(VcfHeaderUtils.STANDARD_FILE_VERSION + "=VCFv4.0");
		data.add(VcfHeaderUtils.STANDARD_FILE_DATE  + "=20150819");
		data.add(VcfHeaderUtils.STANDARD_UUID_LINE + "=abcd_12345678_xzy_999666333");
		data.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXTERN-MELA-20140526-102\tEXTERN-MELA-20140526-07");       
		data.add("chr1\t248845\trs2000390\tC\tT\t2370.26\tPASS\t.\tGT:AD:NNS\t1/1:0,66:52\t1/1:1,74:61");
		data.add("chr1\t248945\trs11485825\tG\tA\t1469.90\tPASS\t.\tGT:AD:NNS\t1/1:0,38:36\t1/1:0,60:47");
		data.add("chr1\t249065\t.\tG\tA\t309.32\tCOVN12\tSOMATIC\tGT:AD:NNS\t.:.:0\t0/1:3,10:8");
		  
		try(BufferedWriter out = new BufferedWriter(new FileWriter(file));) {          
			for (final String line : data)   out.write(line +"\n");                  
		}  
	}
}
