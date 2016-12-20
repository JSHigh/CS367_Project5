import java.util.Random;

///////////////////////////////////////////////////////////////////////////////
//
//Title:            TestSort
//Files:            TestSort.java, ComparisonSort.java, SortObject.java, Questions.txt
//Semester:         CS 367 Fall 2016
//
//Author:           Justin High
//Email:            jshigh@wisc.edu
//CS Login:         high
//Lecturer's Name:  Charles Fischer
//Lab Section:      004
//
///////////////////////////////////////////////////////////////////////////////
//
//Pair Partner:     Aaron Gordner
//Email:            agord@wisc.edu
//CS Login:         gordner
//Lecturer's Name:  Charles Fischer
//Lab Section:      004
//////////////////////////////////////////////////////////////////////////////
public class TestSort {

    /**
     * Main method to run the ComparisonSort class.
     * 
     * @param args  a two-value array: first the number of items in the input
     *              array, then the random number seed (integer)to use in
     *              generating values
     */
    public static void main(String[] args) {

        if (args.length != 2) {
            System.err.println("Expected 2 but got " + args.length);
            System.err.println("Arguments expected:");
            System.err.println("  # items in input array");
            System.err.println("  random # seed");
            System.exit(1);
        }
        int arrSize = Integer.parseInt(args[0]);
        int seed = Integer.parseInt(args[1]);

        System.out.println("Parameters used:");
        System.out.println("  # items in input array: " + arrSize);
        System.out.println("  random # seed: " + seed);

        // Create the input array of unsorted objects.
        SortObject[] arr = new SortObject[arrSize];

        // It is important to give the seed so you can reproduce results.
        Random random = new Random(seed);
        for (int k = 0; k < arrSize; k++)
            arr[k] = new SortObject(random.nextInt());

        // Run all the sorts on the array of random integers.
        ComparisonSort.runAllSorts(arr);
    }
}