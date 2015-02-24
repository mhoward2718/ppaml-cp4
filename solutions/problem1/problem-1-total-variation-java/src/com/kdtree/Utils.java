package com.kdtree;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by eecs on 10/25/2014.
 */
public class Utils {
    public static double median(ArrayList<DataPoint> A, ArrayList<DataPoint> B, int dim)
    {
        int totalSize = A.size()+B.size();
        int middle = totalSize/2;
        double [] array = new double[totalSize];

        int i,j;

        for(i=0;i<A.size();i++)
        {
            array[i] = A.get(i).getValues()[dim];
        }
        for(j=0;j<B.size();j++,i++)
        {
            array[i] = B.get(j).getValues()[dim];
        }

        Arrays.sort(array);
        if(totalSize %2 == 1)
            return array[middle];
        else
            return (array[middle-1] + array[middle])/2.0;
    }

    public static DataSet [] divideElements(ArrayList<DataPoint> dataA, ArrayList<DataPoint> dataB, int dim, double med)
    {
      ArrayList<DataPoint> leftDataA = new ArrayList<DataPoint>();
      ArrayList<DataPoint> leftDataB = new ArrayList<DataPoint>();
      ArrayList<DataPoint> rightDataA = new ArrayList<DataPoint>();
      ArrayList<DataPoint> rightDataB = new ArrayList<DataPoint>();

      DataSet [] leftrightData = new DataSet[2];

      for(int i=0;i<dataA.size();i++)
      {
         if(dataA.get(i).getValues()[dim] <= med)
         {
            leftDataA.add(dataA.get(i));
         }
         else
         {
            rightDataA.add(dataA.get(i));
         }
      }
      for(int i=0;i<dataB.size();i++)
        {
            if(dataB.get(i).getValues()[dim] <= med)
            {
                leftDataB.add(dataB.get(i));
            }
            else
            {
                rightDataB.add(dataB.get(i));
            }
        }
      leftrightData[0] = new DataSet(leftDataA,leftDataB);
      leftrightData[1] = new DataSet(rightDataA,rightDataB);

      return leftrightData;
    }
}
