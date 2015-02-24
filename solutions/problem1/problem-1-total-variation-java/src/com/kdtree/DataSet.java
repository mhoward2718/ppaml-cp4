package com.kdtree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by eecs on 10/24/2014.
 */
public class DataSet {
    ArrayList<DataPoint> dataA;
    ArrayList<DataPoint> dataB;

    public DataSet()
    {
        dataA = new ArrayList<DataPoint>();
        dataB = new ArrayList<DataPoint>();
    }

    public DataSet(ArrayList<DataPoint> A, ArrayList<DataPoint> B)
    {
        dataA = new ArrayList<DataPoint>(A);
        dataB = new ArrayList<DataPoint>(B);
    }

    public ArrayList<DataPoint> getDataA() {
        return dataA;
    }
    public ArrayList<DataPoint> getDataB() { return dataB; }

    // read a data file
    // data file should contain comma separated numbers
    // # columns = # dimensions, # rows = # data points
    // # columns should be equal in both data files
    private ArrayList<DataPoint> ReadFile(String filename)
    {
       String cvsSplitBy = ",";
       String line;

       // a data point is an array of values
       DataPoint dataP;

       // a data file has several of these data points
       ArrayList<DataPoint> data = new ArrayList<DataPoint>();
       File f = new File(filename);
       if(!f.exists())
       {
           System.out.println("Missing File: " + filename);
           System.out.println("EXIT!");
           System.exit(1);
       }
       try {
           BufferedReader br = new BufferedReader(new FileReader(filename));
           while ((line = br.readLine()) != null)
           {
               String[] cols = line.split(cvsSplitBy);
               double [] data_vals = new double[cols.length];

               for(int i=0;i<cols.length;i++)
                   data_vals[i] = Double.parseDouble(cols[i]);

               dataP = new DataPoint();

               dataP.setValues(data_vals);
               data.add(dataP);
           }
       }catch(IOException e)
       {
           System.out.println("Error reading file: " + filename + ". EXITING!");
           System.exit(1);
       }
       return data;
    }

    // data file should contain comma separated numbers
    // # columns = # dimensions, # rows = # data points
    // # columns should be equal in both data files
    public void ReadData(String fileA, String fileB)
    {
       dataA = ReadFile(fileA);
       dataB = ReadFile(fileB);

       // if the data dimensions do not match, then exit
       if(dataA.get(0).getDim() != dataB.get(0).getDim())
       {
           System.out.println("Data dimensions does not match! EXITING!");
           System.exit(1);
       }
    }
}
