package utils;

public class Region implements Comparable<Region> {
    private String chr;
    private int    start;
    private int    end;

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public Region(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    public Region(String chr, String start, String end) {
        this.chr = chr;
        this.start = Integer.parseInt(start);
        this.end = Integer.parseInt(end);
    }

    @Override
    public String toString() {
        return chr + "\t" + start + "\t" + end;
    }

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

    @Override
    public int hashCode() {
        return this.chr.hashCode() + this.start + this.end;
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

}
