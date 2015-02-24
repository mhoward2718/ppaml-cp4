package com.kdtree;

import java.util.Random;

/**
 * Created by eecs on 10/24/2014.
 */
public class KDTree {

    public static double BuildTree(DataSet data, int minLeaf, int NA, int NB)
    {
        int NDataA = data.getDataA().size();
        int NDataB = data.getDataB().size();

        DataSet leftData = new DataSet();
        DataSet rightData = new DataSet();
        DataSet [] leftrightData = new DataSet[2];

        if((NDataA < 2*minLeaf) || (NDataB < 2*minLeaf))
        {
           // not enough points to split, so compute the variation distance at this leaf
           //double val = 0.5*(Math.abs(((double)NDataA/NA) - ((double)NDataB/NB)));
           //return val;
           return 0.5*(Math.abs(((double)NDataA/NA) - ((double)NDataB/NB)));
        }
        else
        {
            // select a dimension at random
            int minDim = 0;
            int maxDim = data.getDataA().get(0).getDim() - 1;
            Random rnd = new Random();
            int randomDim = rnd.nextInt((maxDim - minDim) + 1) + minDim;
            //randomDim = 4;
            //int randomDim = rnd.nextInt((maxDim - minDim) + 1) + minDim;
            // Compute the median of the union of the two data sets along this dimension
            double med = Utils.median(data.getDataA(), data.getDataB(), randomDim);
            //split the data according to the median
            leftrightData = Utils.divideElements(data.getDataA(), data.getDataB(), randomDim, med);
            leftData = leftrightData[0];
            rightData = leftrightData[1];

            if((leftData.dataA.size() == 0 && leftData.dataB.size() == 0) ||(rightData.dataA.size() == 0 && rightData.dataB.size() == 0))
            {
                return 0.5*(Math.abs(((double)NDataA/NA) - ((double)NDataB/NB)));
            }

            //double leftval = BuildTree(leftData, minLeaf, NA, NB);
            //double rightval = BuildTree(rightData, minLeaf, NA, NB);
            //double totalval = leftval + rightval;
            //return totalval;
            // combine the results of the left and right recursive calls
            return BuildTree(leftData, minLeaf, NA, NB) + BuildTree(rightData, minLeaf, NA, NB);
        }
    }
}