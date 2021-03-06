package au.edu.qimr.panel.model;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.ChrPointPosition;
import org.qcmg.common.model.ChrPosition;
import org.qcmg.common.model.ChrRangePosition;
import org.qcmg.common.util.ChrPositionUtils;
import org.qcmg.common.util.Constants;

public class PositionChrPositionMap {
	
	private static final QLogger logger = QLoggerFactory.getLogger(PositionChrPositionMap.class);
	private final Map<ChrPosition, LongRange> chrPosToPositionRange = new HashMap<>();
	
	
	public PositionChrPositionMap() {}
	
	public void loadMap(List<String> chrPosStartPos) {
		for (String s : chrPosStartPos) {
			if (s.startsWith("##chr") || s.startsWith("##GL000")) {
				String [] params = s.substring(2).split(Constants.COLON_STRING);
				int chrLength = Integer.parseInt(params[1]);
				long startOffset = Long.parseLong(params[2]);
				ChrPosition cp = new ChrRangePosition(params[0], 1, chrLength);
				LongRange range = new LongRange(startOffset, startOffset + chrLength -1);
				chrPosToPositionRange.put(cp,  range);
			}
		}
	}
	
	/**
	 * Return the start position of the ChrPosition object as a long primative
	 * @param cp
	 * @return
	 */
	public long getLongStartPositionFromChrPosition(ChrPosition cp) {
		/*
		 * loop through the keys in the map
		 */
		long startPosition = -1;
		for (Entry<ChrPosition, LongRange> entry : chrPosToPositionRange.entrySet()) {
			if (ChrPositionUtils.doChrPositionsOverlap(entry.getKey(), cp)) {
				startPosition = entry.getValue().from + cp.getStartPosition() - 1;
				break;
			}
		}
		return startPosition;
	}
	
	public ChrPosition getChrPositionFromLongPosition(long position) {
		/*
		 * Need to loop through our map values, and check each one to see if the position falls within the range.
		 * Should only every have 1 range that encompases a position....
		 */
		
		for (Entry<ChrPosition, LongRange> entry : chrPosToPositionRange.entrySet()) {
			LongRange lr = entry.getValue();
			if (lr.isPositionWithinRange(position)) {
				/*
				 * Calculate position within contig
				 */
				long positionWithinContig = (position - lr.getFrom()) + 1;
				if (positionWithinContig > Integer.MAX_VALUE) {
					// oops
					logger.warn("positionWithinContig can't be cast to int without overflow!!!");
				}
				
				return ChrPointPosition.valueOf(entry.getKey().getChromosome(), (int) positionWithinContig);
			}
		}
		logger.warn("Could not find ChrPosition for postion: " + position);
		return null;
	}
	
	public ChrPosition getBufferedChrPositionFromLongPosition(long position, int length, int buffer) {
		/*
		 * Need to loop through our map values, and check each one to see if the position falls within the range.
		 * Should only every have 1 range that encompases a position....
		 */
		
		for (Entry<ChrPosition, LongRange> entry : chrPosToPositionRange.entrySet()) {
			LongRange lr = entry.getValue();
			if (lr.isPositionWithinRange(position)) {
				/*
				 * Calculate position within contig
				 */
				long positionWithinContig = (position - lr.getFrom()) + 1;
				if (positionWithinContig > Integer.MAX_VALUE) {
					// oops
					logger.warn("positionWithinContig can't be cast to int without overflow!!!");
				}
				
				return new ChrRangePosition(entry.getKey().getChromosome(), Math.max(1, (int) positionWithinContig - buffer), (int) positionWithinContig + length + buffer);
			}
		}
		logger.warn("Could not find ChrPosition for postion: " + position);
		return null;
	}
	
	private static class LongRange {
		private final long from;
		private final long to;
		
		public LongRange(long from, long to) {
			this.from = from;
			this.to = to;
		}
		public boolean isPositionWithinRange(long position) {
			return from <= position && position <= to;
		}
		public long getFrom() {
			return from;
		}
	}

}
