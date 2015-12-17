package au.edu.qimr.clinvar.model;

import gnu.trove.list.array.TIntArrayList;

import org.qcmg.common.log.QLogger;
import org.qcmg.common.log.QLoggerFactory;
import org.qcmg.common.model.ChrPosition;

public class Fragment implements Comparable<Fragment> {
	
	private final static int MAX_NUMBER_OF_DIFFERENCES = 16;		// increase this if we want to allow reads containing more than 16 bases diff from sequence
	private final static QLogger logger = QLoggerFactory.getLogger(Fragment.class);
	
	private final int id;
	private final String fragment;
	private final int fragmentLength;
	private final  ChrPosition bestTiledLocation;
	private int fsCount;
	private int rsCount;
	private ChrPosition actualLocation;
	private String [] smithWatermanDiffs;
	private final TIntArrayList overlapDistribution;
	
	public Fragment(int id,String sequence, int fsCount, int rsCount, ChrPosition bestTiledLocation, TIntArrayList overlapDist) {
		this.id = id;
		this.fragment = sequence;
		this.fragmentLength = this.fragment.length();
		this.fsCount = fsCount;
		this.rsCount = rsCount;
		this.bestTiledLocation = bestTiledLocation;
		this.overlapDistribution = overlapDist;
	}
	
	public void setForwardStrandCount(int fsCount) {
		this.fsCount = fsCount;
	}
	public void setReverseStrandCount(int rsCount) {
		this.rsCount = rsCount;
	}
	
	public void setSWDiffs(String [] diffs) {
		this.smithWatermanDiffs = diffs;
	}
	public String [] getSmithWatermanDiffs() {
		return smithWatermanDiffs;
	}
	
	public String getSequence() {
		return fragment;
	}
	
	public int getLength() {
		return fragmentLength;
	}
	
	public int getRecordCount() {
		return rsCount + fsCount;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((fragment == null) ? 0 : fragment.hashCode());
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
		Fragment other = (Fragment) obj;
		if (fragment == null) {
			if (other.fragment != null)
				return false;
		} else if (!fragment.equals(other.fragment))
			return false;
		return true;
	}

	@Override
	public int compareTo(Fragment b) {
		return  b.getRecordCount() - getRecordCount();
	}

	public int getId() {
		return id;
	}

	public void setActualPosition(ChrPosition actualCP) {
		this.actualLocation = actualCP;
	}
	public ChrPosition  getActualPosition() {
		return actualLocation;
	}

	public ChrPosition getBestTiledLocation() {
		return bestTiledLocation;
	}

	public int getFsCount() {
		return fsCount;
	}

	public int getRsCount() {
		return rsCount;
	}
	
	public TIntArrayList getOverlapDistribution() {
		return overlapDistribution;
	}

	public void addOverlapDistribution(TIntArrayList overlapDistribution2) {
		overlapDistribution.addAll(overlapDistribution2);
	}

}