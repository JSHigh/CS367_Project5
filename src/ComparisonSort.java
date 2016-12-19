import javax.xml.transform.Templates;
import java.util.Calendar;


/**
 * This class implements six different comparison sorts (and an optional
 * seventh sort for extra credit):
 * <ul>
 * <li>selection sort</li>
 * <li>insertion sort</li>
 * <li>merge sort</li>
 * <li>quick sort</li>
 * <li>heap sort</li>
 * <li>selection2 sort</li>
 * <li>(extra credit) insertion2 sort</li>
 * </ul>
 * It also has a method that runs all the sorts on the same input array and
 * prints out statistics.
 */

public class ComparisonSort {

    /**
     * Sorts the given array using the selection sort algorithm. You may use
     * either the algorithm discussed in the on-line reading or the algorithm
     * discussed in lecture (which does fewer data moves than the one from the
     * on-line reading). Note: after this method finishes the array is in sorted
     * order.
     * 
     * @param <E>  the type of values to be sorted
     * @param A    the array to sort
     */
    public static <E extends Comparable<E>> void selectionSort(E[] A) {
        int j, k, minIndex;
        E min;
        int N = A.length;

        for (k = 0; k < N; k++) {
            min = A[k];
            minIndex = k;
            for (j = k+1; j < N; j++) {
                if (A[j].compareTo(min) < 0) {
                    min = A[j];
                    minIndex = j;
                }
            }
            A[minIndex] = A[k];
            A[k] = min;
        }
    }

    /**
     * Sorts the given array using the insertion sort algorithm. Note: after
     * this method finishes the array is in sorted order.
     * 
     * @param <E>  the type of values to be sorted
     * @param A    the array to sort
     */
    public static <E extends Comparable<E>> void insertionSort(E[] A) {
    	int k, j;
        E tmp;
        int N = A.length;
          
        for (k = 1; k < N; k++) {
            tmp = A[k];
            j = k - 1;
            while ((j >= 0) && (A[j].compareTo(tmp) > 0)) {
                A[j+1] = A[j]; // move one value over one place to the right
                j--;
            }
        A[j+1] = tmp;    // insert kth value in correct place relative to previous values
        }
    }

    /**
     * Sorts the given array using the merge sort algorithm. Note: after this
     * method finishes the array is in sorted order.
     * 
     * @param <E>  the type of values to be sorted
     * @param A    the array to sort
     */
    public static <E extends Comparable<E>> void mergeSort(E[] A) {
    	mergeAux(A, 0, A.length - 1); // call the aux. function to do all the work
    }
    
    /**
     * Helper function to recursively perform a merge sort
     * @param A - array to be sorted
     * @param low - starting index to sort
     * @param high - ending index to sort
     */
    @SuppressWarnings("unchecked")
	private static <E extends Comparable<E>> void mergeAux(E[] A, int low, int high) {
        // base case
        if (low == high) return;
     
        // recursive case
        // Step 1: Find the middle of the array (conceptually, divide it in half)
        int mid = (low + high) / 2;
         
        // Steps 2 and 3: Sort the 2 halves of A
        mergeAux(A, low, mid);
        mergeAux(A, mid+1, high);
     
        // Step 4: Merge sorted halves into an auxiliary array
        E[] tmp = (E[])(new Comparable[high-low+1]);
        int left = low;    // index into left half
        int right = mid+1; // index into right half
        int pos = 0;       // index into tmp
         
        while ((left <= mid) && (right <= high)) {
        // choose the smaller of the two values "pointed to" by left, right
        // copy that value into tmp[pos]
        // increment either left or right as appropriate
        // increment pos
        if (A[left].compareTo(A[right]) <= 0) {
          	tmp[pos] = A[left];
           	left++;
        }
        else {
        	tmp[pos] = A[right];
            right++;
        }
        pos++;
        }
        
        // when one of the two sorted halves has "run out" of values, but
        // there are still some in the other half, copy all the remaining 
        // values to tmp
        // Note: only 1 of the next 2 loops will actually execute
        while (left <= mid) {
        	tmp[pos] = A[left];  
            left++;
            pos++;
        }
        while (right <= high) { 
        	tmp[pos] = A[right]; 
            right++;
            pos++;
        }
     
        // all values are in tmp; copy them back into A
        System.arraycopy(tmp, 0, A, low, tmp.length);
    }

    /**
     * Sorts the given array using the quick sort algorithm, using the median of
     * the first, last, and middle values in each segment of the array as the
     * pivot value. Note: after this method finishes the array is in sorted
     * order.
     * 
     * @param <E>  the type of values to be sorted
     * @param A   the array to sort
     */
    public static <E extends Comparable<E>> void quickSort(E[] A) {
    	quickAux(A, 0, A.length-1);
    }
    
    /**
     * recursive helper function to sort an array using quick sort
     * @param A - array to be sorted
     * @param low - starting value to sort
     * @param high - ending value to sort
     */
    private static <E extends Comparable<E>> void quickAux(E[] A, int low, int high) {
    	if (high - low < 2) {
    		//base case of two items or fewer, sort manually
    		if (high - low == 1) {
    			if (A[low].compareTo(A[high]) > 0) {
    				swap(A, low, high);
    			}
    		}
    	}
    	
    	else {
    		int right = partition(A, low, high);
            quickAux(A, low, right);
            quickAux(A, right+2, high);
    	}
    }
    
    /**
     * partitions array into two halves around a calculated pivot point
     * @param A - array to partition
     * @param low - start of data
     * @param high - end of data
     * @return position where swapped pivot was moved
     */
    private static <E extends Comparable<E>> int partition(E[] A, int low, int high) {
    	E pivot = medianOfThree(A, low, high); // this does step 1
    	int left = low+1;
    	int right = high-2;
        while ( left <= right ) {
            while (A[left].compareTo(pivot) < 0) left++;
            while (A[right].compareTo(pivot) > 0) right--;
            if (left <= right) {
                swap(A, left, right);
                left++;
                right--;
            }
        }

        swap(A, right+1, high-1); // step 4
        return right;
    }
    
    /**
     * Swaps two values in an array
     * @param a - array of data
     * @param i - index 1 to swap with index 2
     * @param j - index 2 to swap with index 1
     */
	private static <E extends Comparable<E>> void swap(E[] a, int i, int j) {
		E tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
	}
	
	/**
	 * calculates the median of three values low, high, and (low+high)/2
	 * @param a - array of values
	 * @param low - lowest value in array
	 * @param high - highest value in array
	 * @return value in the middle
	 */
	private static <E extends Comparable<E>> E medianOfThree(E[] a, int low, int high) {
		int avg = (low + high) / 2;
		//compare average (middle) and low
		if (a[low].compareTo(a[avg]) > 0) {
			swap(a, low, avg);
		}
		
		//compare low and high
		if (a[low].compareTo(a[high]) > 0) {
			swap(a, low, high);
		}
		
		//compare average and high
		if (a[avg].compareTo(a[high]) > 0) {
			swap(a, avg, high);
		}
		
		return a[high -1];
	}

	/**
     * Sorts the given array using the heap sort algorithm outlined below. Note:
     * after this method finishes the array is in sorted order.
     * <p>
     * The heap sort algorithm is:
     * </p>
     * 
     * <pre>
     * for each i from 1 to the end of the array
     *     insert A[i] into the heap (contained in A[0]...A[i-1])
     *     
     * for each i from the end of the array up to 1
     *     remove the max element from the heap and put it in A[i]
     * </pre>
     * 
     * @param <E>  the type of values to be sorted
     * @param A    the array to sort
     */
    @SuppressWarnings("unchecked")
	public static <E extends Comparable<E>> void heapSort(E[] A) {
    	int length = A.length;
    	int size = 0;
    	E[] heap = (E[])(new Comparable[length + 1]);
    	
    	//for each i from 1 to end of the array, insert A[i] into heap
    	for (int i = 0; i < length; i++) {
    		heap[size++] = A[i];
    		
    		// Heapify by swapping the value up
    		int child = size;
    		while (heap[child / 2] != null && heap[child / 2].compareTo(heap[child]) < 0) {
    			// Swap the value up because the parent is less
    			E temp = heap[child / 2];
    			heap[child / 2] = heap[child];
    			heap[child] = temp;

    			// Do we need to swap again?
    			child = child / 2;
    		}
   		}

   		// for each i from the end of the array up to 1
   		// remove the max element from the heap and put it in A[i]
   		for (int i = length - 1; i >= 0; i--) {
   			// Save the root as the value to put at the end of the array
   			A[i] = heap[1];
			// Set the last child as the root
			heap[1] = heap[size];
    		// Set the old last child as null
    		heap[size--] = null;

    		// Heapify by swapping down
    		int parent = 1;
    		while (parent * 2 + 1 < heap.length && ((heap[parent * 2] != null && heap[parent * 2].compareTo(heap[parent]) > 0) || (heap[parent * 2 + 1] != null && heap[parent * 2 + 1].compareTo(heap[parent]) > 0))) {
    			// Swap the parent with the child if the children are bigger
    			E temp = heap[parent];
    			
    			// If both children are bigger, pick the biggest and swap
    			if (heap[parent * 2] != null && heap[parent * 2 + 1] != null) {
    				if (heap[parent * 2].compareTo(heap[parent * 2 + 1]) > 0) {
    					// The left is bigger, swap with the parent
    					heap[parent] = heap[parent * 2];
    					heap[parent * 2] = temp;
    					parent *= 2;
    				} 
    				
    				else {
    					// The right is bigger, swap with the parent
    					heap[parent] = heap[parent * 2 + 1];
    					heap[parent * 2 + 1] = temp;
    					parent = parent * 2 + 1;
    				}
  				} 
    			
    			else if (heap[parent * 2] != null) {
    				// Only the left child is bigger swap with the parent
    				heap[parent] = heap[parent * 2];
    				heap[parent * 2] = temp;
    				parent *= 2;
  				} 
    			
    			else {
    				// Only the right child is bigger, swap with the parent
    				heap[parent] = heap[parent * 2 + 1];
    				heap[parent * 2 + 1] = temp;
    				parent = parent * 2 + 1;
   				}
   			}
    	}
    }

    /**
     * Sorts the given array using the selection2 sort algorithm outlined
     * below. Note: after this method finishes the array is in sorted order.
     * <p>
     * The selection2 sort is a bi-directional selection sort that sorts
     * the array from the two ends towards the center. The selection2 sort
     * algorithm is:
     * </p>
     * 
     * <pre>
     * begin = 0, end = A.length-1
     * 
     * // At the beginning of every iteration of this loop, we know that the 
     * // elements in A are in their final sorted positions from A[0] to A[begin-1]
     * // and from A[end+1] to the end of A.  That means that A[begin] to A[end] are
     * // still to be sorted.
     * do
     *     use the MinMax algorithm (described below) to find the minimum and maximum 
     *     values between A[begin] and A[end]
     *     
     *     swap the maximum value and A[end]
     *     swap the minimum value and A[begin]
     *     
     *     ++begin, --end
     * until the middle of the array is reached
     * </pre>
     * <p>
     * The MinMax algorithm allows you to find the minimum and maximum of N
     * elements in 3N/2 comparisons (instead of 2N comparisons). The way to do
     * this is to keep the current min and max; then
     * </p>
     * <ul>
     * <li>take two more elements and compare them against each other</li>
     * <li>compare the current max and the larger of the two elements</li>
     * <li>compare the current min and the smaller of the two elements</li>
     * </ul>
     * 
     * @param <E>  the type of values to be sorted
     * @param A    the array to sort
     */
    public static <E extends Comparable<E>> void selection2Sort(E[] A) {
        // TODO: implement this sorting algorithm
    	int begin = 0;
    	int end = A.length - 1;
    }

    
    /**
     * <b>Extra Credit:</b> Sorts the given array using the insertion2 sort 
     * algorithm outlined below.  Note: after this method finishes the array 
     * is in sorted order.
     * <p>
     * The insertion2 sort is a bi-directional insertion sort that sorts the 
     * array from the center out towards the ends.  The insertion2 sort 
     * algorithm is:
     * </p>
     * <pre>
     * precondition: A has an even length
     * left = element immediately to the left of the center of A
     * right = element immediately to the right of the center of A
     * if A[left] > A[right]
     *     swap A[left] and A[right]
     * left--, right++ 
     *  
     * // At the beginning of every iteration of this loop, we know that the elements
     * // in A from A[left+1] to A[right-1] are in relative sorted order.
     * do
     *     if (A[left] > A[right])
     *         swap A[left] and A[right]
     *  
     *     starting with with A[right] and moving to the left, use insertion sort 
     *     algorithm to insert the element at A[right] into the correct location 
     *     between A[left+1] and A[right-1]
     *     
     *     starting with A[left] and moving to the right, use the insertion sort 
     *     algorithm to insert the element at A[left] into the correct location 
     *     between A[left+1] and A[right-1]
     *  
     *     left--, right++
     * until left has gone off the left edge of A and right has gone off the right 
     *       edge of A
     * </pre>
     * <p>
     * This sorting algorithm described above only works on arrays of even 
     * length.  If the array passed in as a parameter is not even, the method 
     * throws an IllegalArgumentException
     * </p>
     *
     * @param  A the array to sort
     * @throws IllegalArgumentException if the length or A is not even
     */    
    public static <E extends Comparable<E>> void insertion2Sort(E[] A) { 
        // TODO: implement this sorting algorithm 
    }

    /**
     * Internal helper for printing rows of the output table.
     * 
     * @param sort          name of the sorting algorithm
     * @param compares      number of comparisons performed during sort
     * @param moves         number of data moves performed during sort
     * @param milliseconds  time taken to sort, in milliseconds
     */
    private static void printStatistics(String sort, int compares, int moves,
                                        long milliseconds) {
        System.out.format("%-23s%,15d%,15d%,15d\n", sort, compares, moves, 
                          milliseconds);
    }

    /**
     * Sorts the given array using the six (seven with the extra credit)
     * different sorting algorithms and prints out statistics. The sorts 
     * performed are:
     * <ul>
     * <li>selection sort</li>
     * <li>insertion sort</li>
     * <li>merge sort</li>
     * <li>quick sort</li>
     * <li>heap sort</li>
     * <li>selection2 sort</li>
     * <li>(extra credit) insertion2 sort</li>
     * </ul>
     * <p>
     * The statistics displayed for each sort are: number of comparisons, 
     * number of data moves, and time (in milliseconds).
     * </p>
     * <p>
     * Note: each sort is given the same array (i.e., in the original order) 
     * and the input array A is not changed by this method.
     * </p>
     * 
     * @param A  the array to sort
     */
    static public void runAllSorts(SortObject[] A) {
        System.out.format("%-23s%15s%15s%15s\n", "algorithm", "data compares", 
                          "data moves", "milliseconds");
        System.out.format("%-23s%15s%15s%15s\n", "---------", "-------------", 
                          "----------", "------------");

        // TODO: run each sort and print statistics about what it did
        long startTime = 0;
        long endTime = 0;
        long ellapsedTime = 0;
        String sortName = "";
        int compares = 0;
        int moves = 0;
        
        // Sort 1: selection sort
        sortName = "Selection Sort";
        SortObject[] selSort = A;
        startTime = System.currentTimeMillis();
        ComparisonSort.selectionSort(selSort);
        endTime = System.currentTimeMillis();
        ellapsedTime = endTime - startTime;
       	compares += SortObject.getCompares();
        SortObject.resetCompares(); // rest the object's compares before the next sort

        ComparisonSort.printStatistics(sortName, compares, moves, ellapsedTime);
        
        // Sort 2: 
    }
}