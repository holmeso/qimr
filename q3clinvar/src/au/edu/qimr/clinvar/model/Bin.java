package au.edu.qimr.clinvar.model;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.string.StringUtils;

public class Bin implements Comparable<Bin> {
	
	private final static int MAX_NUMBER_OF_DIFFERENCES = 16;		// increase this if we want to allow reads containing more than 16 bases diff from sequence
	private final static QLogger logger = QLoggerFactory.getLogger(Bin.class);
	
	private final String sequence;
	private final int length;
	private final int[] counts;
	private final Map<PosBase, AtomicInteger> diffs = new HashMap<>();
	
	public Bin(String sequence, int exactMatches) {
		this.sequence = sequence;
		this.length = this.sequence.length();
		counts = new int[MAX_NUMBER_OF_DIFFERENCES];	
		counts[0] = exactMatches;
	}
	
	public void addSequence(String sequenceToAdd) {
		// need position and base at which this differs from sequence
		
		if (StringUtils.isNullOrEmpty(sequenceToAdd)) {
			throw new IllegalArgumentException("Null or empty sequence passed to Bin.addSequence(): " + sequenceToAdd);
		}
		
		if (sequenceToAdd.length() != length) {
			throw new IllegalArgumentException("Sequence passed to Bin.addSequence() is a different length to bins sequence: " + sequence + ", sequenceToAdd: " + sequenceToAdd);
		}
		
		// walk the 2 strings, noting the position and base of any differences
		
		int totalDiffs = 0;
		for (int i = 0 ; i < length ; i++) {
			char a = sequence.charAt(i);
			char b = sequenceToAdd.charAt(i);
			
			if (a != b) {
				totalDiffs++;
				updateDiffs(b, i);
			}
		}
		// update our count array
		if (totalDiffs < MAX_NUMBER_OF_DIFFERENCES) {
			counts[totalDiffs]++;
		} else {
			logger.warn("More than " + MAX_NUMBER_OF_DIFFERENCES + " differences found between sequence: " + sequence + " and sequenceToAdd: " + sequenceToAdd + ", no of diffs: " + totalDiffs);
		}
	}
	
	public String getSequence() {
		return sequence;
	}
	
	public int getLength() {
		return length;
	}
	
	public String getDifferences() {
		if ( ! diffs.isEmpty()) {
			StringBuilder sb =new StringBuilder();
			for (Entry<PosBase, AtomicInteger> entry : diffs.entrySet()) {
				sb.append(entry.getKey().getPos()).append(":").append(entry.getKey().getBase()).append(":").append(entry.getValue().get()).append(",");
			}
			return sb.toString();
		} else {
			return "";
		}
	}
	
	public int getReadCount() {
		int totalCount = 0;
		for (int i : counts) {
			if (i > 0) {
				totalCount += i;
			}
		}
		return totalCount;
	}

	private void updateDiffs(char b, int i) {
		PosBase pb = new PosBase(i, b);
		AtomicInteger ai = diffs.get(pb);
		if (null == ai) {
			ai = new AtomicInteger(1);
			diffs.put(pb, ai);
		} else {
			ai.incrementAndGet();
		}
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((sequence == null) ? 0 : sequence.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Bin other = (Bin) obj;
		if (sequence == null) {
			if (other.sequence != null)
				return false;
		} else if (!sequence.equals(other.sequence))
			return false;
		return true;
	}

	@Override
	public int compareTo(Bin b) {
		return  b.getReadCount() - getReadCount();
	}

}
