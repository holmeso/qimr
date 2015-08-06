package au.edu.qimr.clinvar.model;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.util.SequenceUtil;

import org.junit.Test;

public class ProbeTest {
	
	@Test
	public void getSequenceProbe7() {
		//public Probe(int id, String dlsoSeq, String dlsoSeqRC, String ulsoSeq, String ulsoSeqRC, int p1Start, int p1End, int p2Start, int p2End, String subseq, int ssStart, int ssEnd, String chr) {
		String subseq = "GGTTCACGCCATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGACTACAGGCACCTGCCACCAGGCCCAGCTACTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTGCTGACCTCCTGATCCACCCGCTTCGGCCTCCCAAAGTGCTGGGATTACGGGCGTGCGCCACTGCACCTGGCGACCTTGTAATTCTTTATGTTTTCCAGTTTTCTTCCTTTGTCTCCCTCTTCCCTCAGCTTCAGCTTTGTAAATGCTTTTGAGTCTTTGAGGGGAATAGTTAAATGAATTGCTTAGTGTTCTCTTTATGGAAGTAGCCATTAAGTTTTTTGTTGTTTGTTTTTATGATAATGATGACTTCCTATTACTTCCATTCTGACGCTAGACATGGATTCCTCTAGGAATATACTGGCATTTTGGCTGTCTACCCAGCTGTTGCTTCCATTTGCCTATTCCTTGAGTGTAAGGCAATTAATAACTTACACTTGTCTTTATGTTCCAGCCTGAAAGAATAGACCCAAGCGCATCACGACAAGGATATGATGTCCGCTCTGATGTCTGGAGTTTGGGGATCACATTGGTATGTTTATGCTGATTCAACCTTGCCACAGTAGCGTAACAATAAGAAATTTAGAAGTGAAAGAAAACTTAATCAGACTTCCCCGTTCGTTAAGAACTATAATCACAGACACTATGGTTTTAAGTTGCTGAAAAAAAAAAAGTATGTATTTATTTACTTTAAAAATCAAATCAGACTTGATTATTTCCCTTGAAGTTATGTGAAGTGTCTAGAGCCCTTTAATGTTTTAACTGGATCTCTTGGCCACAGATAATAGGAGTCAACTATTATCTGCAACTGCACTAAGAATGGGAACAGGACAAGCCAGCTTACCTGCAGTCAATTCATTGATGAAGGCACATAGGCTTCTCCCATCTCAGGCAAGGTGTAGTATGTCCAGTAATACTGTATATAGTGTGTATATTAAGTATACGTATTTTTATATGTGGAAAATAGTTTCTACTTGCAATATTCTATTCACTGAATAGTAAAATCACAAGTGAATGAGCTGAAAGTAATTCAGTTTTTGAGGGACAATATATTAGT";
		Probe p = new Probe(7, "AATTGCCTTACACTCAAGGAATAGG", "CCTATTCCTTGAGTGTAAGGCAATT", "GTTCTTAACGAACGGGGAAGTCTGATT", "AATCAGACTTCCCCGTTCGTTAAGAAC",
				12028556, 12028580, 12028758, 12028784, subseq, 12028086, 12029213, "chr17", false,"name1");
		
		String probeRef = p.getReferenceSequence();
		System.out.println("p7 ref: " + probeRef);
		assertEquals(p.getExpectedFragmentLength(), probeRef.length());
		
		String bin = "CCTATTCCTTGAGTGTAAGGCAATTAATAACTTACACTTGTCTTTATGTTCCAGCCTGAAAGAATAGACCCAAGCGCATCACGACAAGGATATGATGTCCGCTCTGATGTCTGGAGTTTGGGGATCACATTGGTATGTTTATGCTGATTCAACCTTGCCACAGTAGCGTAACAATAAGAAATTTAGAAGTGAAAGAAAACTTAATCAGACTTCCCCGTTCGTTAAGAAC";
		// need to RC bin as we are on the reverse strand
		
		assertEquals(probeRef, bin);
//		assertEquals(probeRef, SequenceUtil.reverseComplement(bin));
	}
	@Test
	public void getSequenceProbe12() {
		//public Probe(int id, String dlsoSeq, String dlsoSeqRC, String ulsoSeq, String ulsoSeqRC, int p1Start, int p1End, int p2Start, int p2End, String subseq, int ssStart, int ssEnd, String chr) {
		String subseq = "ATACCTGGGCATTGAGATTTTTAAAATCTCCCCCAGGTGATTGTAAGGTACAGCAAACCTTGGGAATAGTTGTACTAACTGAAGATTGAGGCACTGGAATAACTAATGTTGGGTAAGATAGGGTAAGAACTTTTTCGCAATAAGGTTTGCAGCTCTTGGGATTCTCACTTGATAATGCTGTTTGGATAAACCATCATCCAAGTTTGAATAAGATCTTATATTTAGAAATGAAAACTGATCTGCAGAATTAGTTTAATCAATAAATGGGAATTTTTAAAAATCTTCTTGGCAATTTGTTTTTATGAATCCCGTTGAAGCTGTGTCTATTGACTACCTAAATCTAGATCATTCTTAATCCAGTGTGTAAAAAACCAGGATGACACAAATGAAAAACTTCAAAAACCTGGAGGTCAGACTATTTTAGTAATTTAGAAAACATTTTTCCCACACATTAATCAGTACTAAAAGAAAAAAGTTAAAACCTATTTAAAATGTGGAAAAATTGCTTCCCAATATTTTAACAGAGAGAGACTGAGAACACACAGCATTGAGTCATCAGGAAAACTGAAGATCTCCCCTGAACAACACTGGGATTTCACTGCAGAGGACTTGAAAGACCTTGGAGAAATTGGACGAGGAGCTTATGGTTCTGTCAACAAAATGGTCCACAAACCAAGTGGGCAAATAATGGCAGTTAAAGTAGGTGATGCCATGATTATTTTTGGTACTTTAATCCATTAGGTGAAATTTCATGGTGCAGTAATACCACTGTTGTTGTGTTCCTACTTTTGTGGTAAATGTGGGTGTTTAAAAAATTGTTTCTCCAACTCCTTTGAGGGTGTTCTGTGTAGAGGTTTCTTATTGGTGGCACATTTTCTATCTCTTTGGAAACATGAGTGTATGAAGTGTGCATGTTGATTGCATTTTTGGCATGATGCATTAACATGTTTTAAACTTCAGGCCCTTCTACTGCCAAGGTGAGTTCAGGCTGGGCGGCTGCACCCCTGGGAGCAGGGCAGTGCTGCACTGAGCCAGGCGGGAGCTGGAAGAAGACGCAGCACACTGGGTTGGGCAAGGTGCGGGGCTAGGGCTAACAGCAGTCTTACTGAAGGTTTCCTGGAAACCACGCACATGCTGTTGCCACTAACCTCAACCTTACTCGGTCCTGACCGGCTCGGCTTCTGTTTGTTTATTTCATCTCTACTCAGTACTGCCCTGTTCCCT";
		Probe p = new Probe(12, "TGAAATTTCATGGTGCAGTAATACCACTG", "CAGTGGTATTACTGCACCATGAAATTTCA", "GCATTGAGTCATCAGGAAAACTGAAGA", "TCTTCAGTTTTCCTGATGACTCAATGC",
				11984692, 11984718, 11984890, 11984918, subseq, 11984148, 11985372, "chr17", true,"name2");
		
		String probeRef = p.getReferenceSequence();
		System.out.println("p12 ref: " + probeRef);
		assertEquals(p.getExpectedFragmentLength(), probeRef.length());
		
		String bin = "CAGTGGTATTACTGCACCATGAAATTTCACCTAATGGATTAAAGTACCAAAAATAATCATGGCATCACCTACTTTAACTGCCATTATTTGCCCACTTGGTTTGTGGACCATTTTGTTGACAGAACCATAAGCTCCTCGTCCAATTTCTCCAAGGTCTTTCAAGTCCTCTGCAGTGAAATCCCAGTGTTGTTCAGGGGAGATCTTCAGTTTTCCTGATGACTCAATGC";
		bin = SequenceUtil.reverseComplement(bin);
		System.out.println("p12 bin: " + bin);
		
		// is bin in subseq
		assertEquals(true, subseq.contains(bin));
		assertEquals(probeRef, bin);
	}
	
	@Test
	public void getPositinOfOffset() {
		String subseq = "AAAAAAAAAABBBBBBBBBBCCCCCCCCCCDDDDDDDDDDEEEEEEEEEEFFFFFFFFFFGGGGGGGGGGHHHHHHHHHH";
		Probe p = new Probe(7, "", "", "", "", 11, 20, 71, 70, subseq, 1, 80, "chrOllie", false,"name1");
		
		assertEquals(1, p.getSubReferencePosition(subseq));
		assertEquals(11, p.getSubReferencePosition("BBBBBBBBBB"));
		assertEquals(21, p.getSubReferencePosition("CCCCCCCCCC"));
		assertEquals(-1, p.getSubReferencePosition("XYZ"));
	}

}
