package utils;

public class Region implements Comparable<Region> {
    private String chr;
    private long   start;
    private long   end;
    private Type   type;

    public enum Type {
        GAIN, LOSS;
    }

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
        this.type = Type.GAIN;
    }

    public Region(String chr, long start, long end, Type type) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.type = type;
    }

    public Region(String chr, String start, String end) {
        this(chr, Long.parseLong(start), Long.parseLong(end));
    }

    public Region(String chr, String start, String end, String type) {
        this(chr, Long.parseLong(start), Long.parseLong(end), Type.valueOf(type));
    }

    public long getLength() {
        return this.getEnd() - this.getStart();
    }

    /**
     * judge if number located at the position between start and end.
     * @param number
     * @param start
     * @param end
     * @return
     */
    private boolean isBetween(long number, long start, long end) {
        if (number >= start && number <= end) {
            return true;
        }
        return false;
    }

    /**
     * judge if this region overlaps with other region
     * @param other
     * @return 
     */
    public boolean isOverlapped(Region other) {
        if (!this.chr.equals(other.chr))
            return false;
        if (isBetween(this.getStart(), other.getStart(), other.getEnd())
            || isBetween(this.getEnd(), other.getStart(), other.getEnd())
            || isBetween(other.getStart(), this.getStart(), this.getEnd())
            || isBetween(other.getEnd(), this.getStart(), this.getEnd())) {
            return true;
        }
        return false;
    }

    /**
     * judge if this region overlaps with other region. If type of this region does not meet other region's type, return false.
     * @param other
     * @return 
     */
    public boolean isOverlappedWithType(Region other) {
        if (!this.chr.equals(other.chr) || this.type != other.type)
            return false;
        if (isBetween(this.getStart(), other.getStart(), other.getEnd())
            || isBetween(this.getEnd(), other.getStart(), other.getEnd())
            || isBetween(other.getStart(), this.getStart(), this.getEnd())
            || isBetween(other.getEnd(), this.getStart(), this.getEnd())) {
            return true;
        }
        return false;
    }

    /**
     * return the overlap length between this region and other region.
     * @param other
     * @return overlapLength
     */
    public long getOverlapLength(Region other) {
        if (!this.isOverlapped(other))
            return 0;
        long overlapLength;
        if (other.getStart() >= this.getStart()) {
            if (other.getEnd() >= this.getEnd()) {
                overlapLength = this.getEnd() - other.getStart();
            } else {
                overlapLength = other.getEnd() - other.getStart();
            }
        } else {
            if (other.getEnd() >= this.getEnd()) {
                overlapLength = this.getEnd() - this.getStart();
            } else {
                overlapLength = other.getEnd() - this.getStart();
            }
        }
        return overlapLength;
    }

    /**
     * return the overlap length between this region and other region. If type of this region does not meet other region's type, return 0.
     * @param other
     * @return overlapLength
     */
    public long getOverlapLengthWithType(Region other) {
        if (!this.isOverlappedWithType(other))
            return 0;
        long overlapLength;
        if (other.getStart() >= this.getStart()) {
            if (other.getEnd() >= this.getEnd()) {
                overlapLength = this.getEnd() - other.getStart();
            } else {
                overlapLength = other.getEnd() - other.getStart();
            }
        } else {
            if (other.getEnd() >= this.getEnd()) {
                overlapLength = this.getEnd() - this.getStart();
            } else {
                overlapLength = other.getEnd() - this.getStart();
            }
        }
        return overlapLength;
    }

    /**
     * return sum of the base length of two regions overlapped.
     * @param other
     * @return sum of the base length, 0 if the two regions not overlapped.
     */
    public long getOverlapBaseLength(Region other) {
        if (!this.isOverlapped(other))
            return 0;
        long sumBaseLength;
        if (other.getStart() >= this.getStart()) {
            if (other.getEnd() >= this.getEnd()) {
                sumBaseLength = other.getEnd() - this.getStart();
            } else {
                sumBaseLength = this.getEnd() - this.getStart();
            }
        } else {
            if (other.getEnd() >= this.getEnd()) {
                sumBaseLength = other.getEnd() - other.getStart();
            } else {
                sumBaseLength = this.getEnd() - other.getStart();
            }
        }
        return sumBaseLength;
    }

    /**
     * return sum of the base length of two regions overlapped. If type of this region does not meet other region's type, return 0.
     * @param other
     * @return sum of the base length, 0 if the two regions not overlapped.
     */
    public long getOverlapBaseLengthWithType(Region other) {
        if (!this.isOverlapped(other))
            return 0;
        long sumBaseLength;
        if (other.getStart() >= this.getStart()) {
            if (other.getEnd() >= this.getEnd()) {
                sumBaseLength = other.getEnd() - this.getStart();
            } else {
                sumBaseLength = this.getEnd() - this.getStart();
            }
        } else {
            if (other.getEnd() >= this.getEnd()) {
                sumBaseLength = other.getEnd() - other.getStart();
            } else {
                sumBaseLength = this.getEnd() - other.getStart();
            }
        }
        return sumBaseLength;
    }

    @Override
    public int hashCode() {
        return (int) (this.chr.hashCode() + this.start + this.end + this.type.hashCode());
    }

    @Override
    public boolean equals(Object other) {
        Region region = (Region) other;
        return this.chr.equals(region.chr) && this.start == region.start && this.end == region.end
               && this.type == region.type;
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
            else {
                if (this.type == other.type)
                    return 0;
                else if (this.type == Type.GAIN)
                    return 1;
                else
                    return -1;
            }
        }
    }

    @Override
    public String toString() {
        return type.toString() + "\t" + chr + "\t" + start + "\t" + end;
    }

    public Type getType() {
        return type;
    }

    public void setType(String type) {
        this.type = Type.valueOf(type);
    }

}
