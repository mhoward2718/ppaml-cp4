package com.kdtree;

/**
 * Created by eecs on 10/24/2014.
 */
// a data point is a collection of values
public class DataPoint {
    double [] values;

    public DataPoint() {}
    public void setValues(double[] values) {
        this.values = values;
    }
    public int getDim() {
        return this.values.length;
    }
    public double[] getValues() {
        return values;
    }
}
