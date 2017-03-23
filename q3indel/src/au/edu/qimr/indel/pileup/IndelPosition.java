/**
 * © Copyright The University of Queensland 2010-2014.  This code is released under the terms outlined in the included LICENSE file.
 */
package au.edu.qimr.indel.pileup;

import java.util.ArrayList;
import java.util.List;

import org.qcmg.common.model.ChrRangePosition;
import org.qcmg.common.util.Constants;
import org.qcmg.common.util.IndelUtils;
import org.qcmg.common.util.IndelUtils.SVTYPE;
import org.qcmg.common.vcf.VcfFormatFieldRecord;
import org.qcmg.common.vcf.VcfRecord;
import org.qcmg.common.vcf.VcfUtils;
import org.qcmg.common.vcf.header.VcfHeaderUtils;
public class IndelPosition {
	
	private final List<VcfRecord> vcfs ; 	
	private final ChrRangePosition position; 
	private final int indelStart, indelEnd; 
	
	private final SVTYPE mutationType;	
	private IndelPileup tumourPileup;
	private IndelPileup normalPileup;

	/**
	 * retrive information from a vcf record
	 * @param re: a vcf record
	 * @param type
	 * @param indelFileType
	 * @param isGermline
	 */
	public IndelPosition(VcfRecord re, SVTYPE type){
		vcfs = new ArrayList<VcfRecord>(); 
		
		VcfRecord vcf = re; 				
		vcfs.add(re);		
		
		this.mutationType = type;
		String fullChromosome = IndelUtils.getFullChromosome(vcf.getChromosome());
	 	
		
		if (isInsertion()) {
			this.indelStart = vcf.getPosition();
			this.indelEnd = vcf.getPosition()+1;
		} else if (isDeletion()){
			this.indelStart = vcf.getPosition() + 1;
			this.indelEnd = vcf.getPosition() + vcf.getRef().length() - 1;
		} else {
			this.indelStart = -1;
			this.indelEnd = -1; 
		}
		
		int start = vcf.getPosition();
		int end = vcf.getChrPosition().getEndPosition();
		position = new ChrRangePosition(fullChromosome, start, end);
		
	}

	//next job: check all vcfs are same type, start and end
	public IndelPosition(List<VcfRecord> res, SVTYPE type) throws Exception{
		this(res.get(0), type);		
		//append all vcfs
		vcfs.clear();
		vcfs.addAll(res);

	}
	
	public IndelPosition(VcfRecord re){
		this(re, IndelUtils.getVariantType(re.getRef(), re.getAlt() ));		
	}
	
	public IndelPosition(List<VcfRecord> res ) throws Exception{
		this(res, IndelUtils.getVariantType(res.get(0).getRef(), res.get(0).getAlt() ) );
				
	}
	
	public VcfRecord getIndelVcf(int index){		
		return vcfs.get(index); 
	}
	
	public List<VcfRecord> getIndelVcfs( ){		
		return vcfs; 
	}
	
	public void addVcf(VcfRecord vcf) throws Exception{
		if(!IndelUtils.getVariantType(vcf.getRef(), vcf.getAlt() ).equals(mutationType) )
			throw new Exception();
				
		vcfs.add(vcf);
	}
	
	public SVTYPE getIndelType(){
		return mutationType; 
	}

	/**
	 * 
	 * @return ChrRangePosition from vcfRecord but convert chr to full chromosome
	 * eg return new ChrRangePosition("chr1", 99, 99 ) for vcf "1 99 . A ATT"
	 *   return new ChrRangePosition("chr1", 99, 101 ) for vcf "1 99 . ATT A"
	 * 
	 */
	public ChrRangePosition getChrRangePosition(){		
		return position; 
	}
	
	/**
	 * 
	 * @return Indel start.
	 *  eg. return 99 for INS "chr 99 . A ATT"
	 *      return 100 for DEL "chr 99 . ATT A"
	 */
	public int getStart() {	
		return indelStart;
	}
	
	/**
	 * 
	 * @return Indel end.
	 *  eg. return 100 for INS "chr 99 . A ATT"
	 *      return 101 for DEL "chr 99 . ATT A"
	 */
	public int getEnd() {		 
		return indelEnd; 
	}
	
	public String getFullChromosome() {		
		return position.getChromosome();	 
    }
	
	public boolean isInsertion() {		
		return  mutationType.equals(SVTYPE.INS); 
	}

	public boolean isDeletion() {
		return  mutationType.equals(SVTYPE.DEL); 		 
	}		
	
	public boolean isComplex() {
		return mutationType.equals(SVTYPE.CTX);
	}
	
	public String getMotif(int index) {	
		
		if(mutationType.equals(SVTYPE.INS))
			return vcfs.get(index).getAlt().substring(1);			 
		else if(mutationType.equals(SVTYPE.DEL))
			return vcfs.get(index).getRef().substring(1);
		
		return null; 
	}
	
	
	public List<String> getMotifs( ) {	
		List<String> motifs = new ArrayList<String>();
		for(int i = 0; i < vcfs.size(); i ++)
			motifs.add(i, getMotif(i));
		
		return motifs; 
	}
	
	
	@Override
    public boolean equals(final Object o) {
	       
        if (!(o instanceof IndelPosition)) return false;
        
        final IndelPosition other = (IndelPosition) o;
             	
    	if (! this.mutationType .equals(other.mutationType))  
    		return false;
    	
    	if( ! this.position.equals(other.position))
    		return false;
    	
     
        return true; 
    }
	    
    @Override
    public int hashCode() { 
				
		return position.hashCode();
   }
    
	@Override
	public String toString() {
		return this.position.toString() + ":" + mutationType + ": variants number: " + vcfs.size();
	}

	public void setPileup(boolean isTumour, IndelPileup pileup){
		if(isTumour)
			 this.tumourPileup = pileup;
		else
			this.normalPileup = pileup;
	}
	
	
	public VcfRecord getPileupedVcf(int index, final int gematic_nns, final float gematic_soi){
		VcfRecord re = vcfs.get(index);		
		
		//not interested in these indels since over coverage		
		if(tumourPileup != null  && tumourPileup.getTotalCount() > 1000){
			re.setFilter( IndelUtils.FILTER_HCOVT);
			return re; 
		}
	 	
		if(normalPileup != null  && normalPileup.getTotalCount() > 1000){
			re.setFilter( IndelUtils.FILTER_HCOVN);
			return re; 
		}		
		
		boolean somatic = (re.getFilter().equals(ReadIndels.FILTER_SOMATIC))? true : false;
 		if(normalPileup != null && somatic){
 			if( normalPileup.getsuportReadCount(index) > gematic_nns )
				somatic = false;
			else if(normalPileup.getInformativeCount() > 0){
				int scount =   normalPileup.getsuportReadCount(index);
				int icount =   normalPileup.getInformativeCount();
				if((100 * scount / icount) >= (gematic_soi * 100)) 
					somatic = false; 
			}		
		}
		
		if( somatic) 
			re.appendInfo(VcfHeaderUtils.INFO_SOMATIC);

		//set default filter as PASS
		re.setFilter(VcfHeaderUtils.FILTER_PASS);
		String td = ".", nd = ".";	
		
		if(tumourPileup != null){ 		
			if(tumourPileup.getTotalCount() > 0)
			 td = String.format("%d,%d,%d,%d[%d,%d],%d[%d],%d,%d,%d", tumourPileup.getstrongsupportNovelStartReadCount(index),tumourPileup.getTotalCount(),tumourPileup.getInformativeCount(), 
					tumourPileup.getstrongSuportReadCount(index),tumourPileup.getforwardsuportReadCount(index),tumourPileup.getbackwardsuportReadCount(index),tumourPileup.getsuportReadCount(index),
					tumourPileup.getsupportNovelStartReadCount(index), tumourPileup.getparticalReadCount(index),tumourPileup.getNearbyIndelCount(),tumourPileup.getNearybySoftclipCount());
			
			if(!somatic && tumourPileup.getTotalCount() < 8)
				VcfUtils.updateFilter(re,  IndelUtils.FILTER_COVT);		 
			if(tumourPileup.getparticalReadCount(index) >= 3 &&
					(100 * tumourPileup.getparticalReadCount(index) / tumourPileup.getTotalCount()) > 10)
				VcfUtils.updateFilter(re,IndelUtils.FILTER_TPART);
			if(somatic && tumourPileup.getsuportReadCount(index) >=3 && tumourPileup.hasStrandBias(index, 0.1))
				VcfUtils.updateFilter(re,IndelUtils.FILTER_TBIAS);			
		}	
		
		//String nd = "ND=0:0:0:0:0:0:0";
		if(normalPileup != null){
			if(normalPileup.getTotalCount() > 0)
				nd = String.format("%d,%d,%d,%d[%d,%d],%d[%d],%d,%d,%d", normalPileup.getstrongsupportNovelStartReadCount(index),normalPileup.getTotalCount(),normalPileup.getInformativeCount(), 
					normalPileup.getstrongSuportReadCount(index),normalPileup.getforwardsuportReadCount(index),normalPileup.getbackwardsuportReadCount(index),normalPileup.getsuportReadCount(index),
					normalPileup.getsupportNovelStartReadCount(index),  normalPileup.getparticalReadCount(index),normalPileup.getNearbyIndelCount(),normalPileup.getNearybySoftclipCount());
			 			
			if(somatic && normalPileup.getTotalCount() < 12)
				VcfUtils.updateFilter(re,IndelUtils.FILTER_COVN12);
			if(!somatic && normalPileup.getTotalCount() < 8)
				VcfUtils.updateFilter(re,IndelUtils.FILTER_COVN8);			
			if(somatic && normalPileup.getsupportNovelStartReadCount(index) > 0)
				VcfUtils.updateFilter(re,IndelUtils.FILTER_MIN);
			if(normalPileup.getparticalReadCount(index) >= 3 &&
					(100 * normalPileup.getparticalReadCount(index) / normalPileup.getTotalCount()) > 5)
				VcfUtils.updateFilter(re,IndelUtils.FILTER_NPART);
			if( !somatic && normalPileup.getsuportReadCount(index) >=3 && normalPileup.hasStrandBias(index, 0.05))
				VcfUtils.updateFilter(re,IndelUtils.FILTER_NBIAS);	
		}
					 
		//check GT	
		//control always on first column and then test, must have control, test	
		List<String> genotypeField =  re.getFormatFields();	
		String normal_gt = (genotypeField.size() > 1)? genotypeField.get(1)+ ":" + nd : nd;
		String tumour_gt =  (genotypeField.size() > 2)? genotypeField.get(2) + ":" + td: td;	
		
		if(normal_gt.startsWith(Constants.MISSING_DATA_STRING + Constants.COLON_STRING))
			normal_gt =  ( normalPileup != null && normalPileup.getTotalCount() > 2 ) ? normal_gt.replaceFirst(".:", "0/0:") : normal_gt.replaceFirst(".:", "./.:");   			
		if( tumour_gt.startsWith(Constants.MISSING_DATA_STRING + Constants.COLON_STRING) )			 			
			tumour_gt = (tumourPileup != null && tumourPileup.getTotalCount() > 2 ) ? tumour_gt.replaceFirst(".:", "0/0:") : tumour_gt.replaceFirst(".:", "./.:");   
			 							
		List<String> field = new ArrayList<String>();
		field.add(0,  (genotypeField.size() > 0)? genotypeField.get(0) + ":" + IndelUtils.FORMAT_ACINDEL : IndelUtils.FORMAT_ACINDEL );
		field.add(1, normal_gt );
		field.add(2, tumour_gt );		
				
		re.setFormatFields( field ); 					 
		IndelPileup pileup = (somatic)? tumourPileup: normalPileup; 				
		
		if(pileup !=null  && pileup.getstrongsupportNovelStartReadCount(index) < 4)
			VcfUtils.updateFilter(re,IndelUtils.FILTER_NNS);
		
		float nn = (pileup == null || pileup.getTotalCount() == 0)? 0 : (float) pileup.getNearbyIndelCount() / pileup.getTotalCount();
		float ss = (pileup == null || pileup.getInformativeCount() == 0)? 0 : (float) pileup.getstrongSuportReadCount(index) / pileup.getInformativeCount();		
		re.appendInfo( String.format( IndelUtils.INFO_NIOC + "=" + ((nn == 0 )? "0" : String.format("%.3f", nn))) );		
		re.appendInfo( String.format( IndelUtils.INFO_SSOI + "=" + ((ss == 0 )? "0" : String.format("%.3f", ss))) );	
				
//		re.appendInfo("SVTYPE=" + this.mutationType.name());
//		re.appendInfo("END=" + indelEnd);		
		re.appendInfo( IndelUtils.INFO_SVTYPE + this.mutationType.name() );
//		re.appendId( IndelUtils.INFO_END + indelEnd );		
		re.appendInfo( VcfHeaderUtils.INFO_MERGE_IN + "=1" );
					
		return re; 	
	}

}
