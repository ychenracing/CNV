package utils;

public class Region implements Comparable<Region> {
    private String chr;
    private long   start;
    private long   end;

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public long getStart() {
        return start;
    }

    public void setStart(long start) {
        this.start = start;
    }

    public long getEnd() {
        return end;
    }

    public void setEnd(long end) {
        this.end = end;
    }

    public Region(String chr, long start, long end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    public Region(String chr, String start, String end) {
        this.chr = chr;
        this.start = Long.parseLong(start);
        this.end = Long.parseLong(end);
    }

    public long getLength() {
        return this.getEnd() - this.getStart();
    }

    /**
     * judge if this region overlaps with other region
     * @param other
     * @return 
     */
    public boolean intersact(Region other) {
        if (!this.chr.equals(other.chr))
            return false;
        if (this.start <= other.start && this.end >= other.end
            || this.start <= other.end && this.end >= other.end
            || this.start >= other.start && this.end <= other.end) {
            return true;
        }
        return false;
    }

    /**
     * return the overlap length between this region and other region.
     * @param other
     * @return overlapLength
     */
    public long overlap(Region other) {
        if (!this.intersact(other))
            return 0;
        long overlapLength;
        if (other.getStart() >= this.getStart()) {
            if (other.getEnd() >= this.getEnd()) {
                overlapLength = this.getEnd() - other.getStart();
            } else {
                overlapLength = other.getEnd() - other.getStart();
            }
        } else {
            overlapLength = other.getEnd() - this.getStart();
        }
        return overlapLength;
    }

    @Override
    public int hashCode() {
        return (int) (this.chr.hashCode() + this.start + this.end);
    }

    @Override
    public boolean equals(Object other) {
        Region region = (Region) other;
        return this.chr.equals(region.chr) && this.start == region.start && this.end == region.end;
    }

    @Override
    public int compareTo(Region other) {
        if (!this.chr.equals(other.chr))
            return this.chr.compareTo(other.chr) < 0 ? -1 : 1;
        if (this.start < other.start)
            return -1;
        else if (this.start > other.start)
            return 1;
        else {
            if (this.end < other.end)
                return -1;
            else if (this.end > other.end)
                return 1;
            else
                return 0;
        }
    }

    @Override
    public String toString() {
        return chr + "\t" + start + "\t" + end;
    }

}
