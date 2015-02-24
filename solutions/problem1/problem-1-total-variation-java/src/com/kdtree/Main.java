package com.kdtree;

public class Main {

    //static String dataFileA = "problem-1-solution-samples.csv";
    //static String dataFileB = "mu2data_1M_50D.csv";
    // data file A and B should contain comma separated numbers
    // # columns = # dimensions, # rows = # data points
    // # columns should be equal in both data files

    static int ITERATIONS =  1000;
    static int minLeaf = 1;

    static DataSet data;

    public static void main(String[] args)
	{
       if(args.length <2)
       {
         System.out.println("Missing Ground Truth or Input File Argument!");
         System.out.println("java -jar totalvar.jar <ground-truth-file-path> <input-file-path>");
         System.out.println("EXIT.");
         System.exit(1);
       }
       if(args.length >= 3)
           try {
               minLeaf = Integer.parseInt(args[2]);
           }
           catch(Exception e)
           {
               System.out.println("minLeaf argument Error. Exit.");
               System.exit(1);
           }
       computeScore(args[0], args[1]);
    }
	
	public static void computeScore(String groundTruth, String input)
	{
	    double [] score = new double[ITERATIONS];
        double sumScore = 0;
        double avgScore = 0;

        data = new DataSet();
        data.ReadData(groundTruth, input);

        // NA = number of points in data set A
        // NB = number of points in data set B
        int NA = data.getDataA().size();
        int NB = data.getDataB().size();

		System.out.println("Beginning score computation...");
		
        for(int i=0;i<ITERATIONS;i++) {
            System.out.println("Starting Iteration: " + i);
            score[i] = KDTree.BuildTree(data, minLeaf, NA, NB);
            sumScore = sumScore + score[i];
            System.out.println("Iteration: " + i + " Score: " + score[i]);
        }
        avgScore = sumScore/ITERATIONS;
        System.out.println("Average score over " + ITERATIONS + " iterations: " + avgScore);
	}
	
}
