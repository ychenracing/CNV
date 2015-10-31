package utils;

public class Pair<T, R> {
    T t;
    R r;

    public Pair() {
    }

    public Pair(T t, R r) {
        this.t = t;
        this.r = r;
    }

    public T getFirst() {
        return t;
    }

    public void setFirst(T t) {
        this.t = t;
    }

    public R getSecond() {
        return r;
    }

    public void setSecond(R r) {
        this.r = r;
    }

    @Override
    public String toString() {
        return this.t.toString() + ":" + this.r.toString();
    }
}
