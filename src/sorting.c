#include<stdio.h>

static int DEBUG=0;

//#define GRAPH_PATH 1
#include<SparseBLAS.h>

/* Function to sort an array using insertion sort*/
void insertionSort(COOE *arr, int n, Comparing comp, Ordering *order)
{
  COOE key;
  int i,  j;
  if (DEBUG) {
    printf("%d insert>: ",n );
    for (int i=0; i<n; i++) printf("(%d,%d,%d) ",
				   arr[i].m,arr[i].n,(int)arr[i].value);
    printf("\n");
  }
  
  for (i = 1; i < n; i++)
    {
      key = arr[i];
      j = i-1;

      /* Move elements of arr[0..i-1], that are
	 greater than key, to one position ahead
	 of their current position */
      while (j >= 0 && comp(arr+j, &key, order)>0)
	{
	  arr[j+1] = arr[j];
	  j = j-1;
	}
      arr[j+1] = key;
    }

    if (DEBUG) {
      printf("%d insert<: ",n);
      for (int i=0; i<n; i++) printf("(%d,%d,%d) ",
				     arr[i].m,arr[i].n,(int)arr[i].value);
      printf("\n");
    }


}


// A utility function to swap two elements
void swap(COOE* a, COOE* b)
{
  COOE t = *a;
  *a = *b;
  *b = t;
}

/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
int partition (COOE *arr, int low, int high, Comparing comp, Ordering *order)
{
  COOE pivot = arr[high]; // pivot
  int i = (low - 1); // Index of smaller element

  for (int j = low; j <= high-1; j++)
    {
      // If current element is smaller than or
      // equal to pivot
      if (comp(arr +j,&pivot,order)<=0)
	{
	  i++; // increment index of smaller element
	  swap(arr +i, arr +j);
	}
    }
  swap(arr +i + 1, arr +high);
  return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSort(COOE *arr, int low, int high, Comparing comp, Ordering *order)
{
  if (low < high )
    {
      if (DEBUG) {
	printf("%d %d quick >: ", low, high);
	for (int i=low; i<=high; i++) printf("(%d,%d,%d) ",
					     arr[i].m,arr[i].n,(int)arr[i].value);
	printf("\n");
      }
      
      
      if (high-low <7) 
	insertionSort(arr+low,high-low+1, comp, order)	;
      else {
      
	/* pi is partitioning index, arr[p] is now
	   at right place */
	int pi = partition(arr, low, high,comp, order);

	if (DEBUG) {
	  printf("%d pivot  quick >: ", pi);
	  printf("(%d,%d,%d) ",arr[pi].m,arr[pi].n,(int)arr[pi].value);
	  printf("\n");
	}

	// Separately sort elements before
	// partition and after partition
	quickSort(arr, low, pi - 1,comp,order);
	quickSort(arr, pi + 1, high,comp,order);
      }

      if (DEBUG) {
	printf("< %d %d quick: ", low, high);
	for (int i=low; i<=high; i++) printf("(%d,%d,%d) ",arr[i].m,arr[i].n,(int)arr[i].value);
	printf("\n");
      }
      
    }
  
}

/***
 *  Using coordinate elements in sparse matrices, the order is row
 *  major then column major (row, col) ~ (row*N + col). To transpose
 *  the matrix we need to change the order so that (col, row) ~ (row +
 *  col*M).
 *  
 *  here we can see two example how the transposing is done by
 *  sorting.
 */     


int roworder(COOE *a, COOE *b, Ordering *order) {

  return (a->m*order->N + a->n)-(b->m*order->N+b->n);
  
}

int colorder(COOE *a, COOE *b, Ordering *order) {

  return (a->m+ order->M *a->n)-(b->m + order->M*b->n);
  
}

void columnsort(COO *M) {
  Ordering order = {M->M, M->N, 1 } ;
  quickSort(M->data,0, M->length-1, colorder, &order);
}

void rowsort(COO *M) {
  Ordering order = {M->M, M->N, 1 } ;
  quickSort(M->data,0, M->length-1, roworder, &order);
}




/* Function to sort an array using insertion sort*/
void insertionSort_b(COOB *arr, int n, ComparingB comp, Ordering *order)
{
  COOB key;
  int i,  j;
  if (DEBUG) {
    printf("%d insert>: ",n );
    for (int i=0; i<n; i++) printf("(%d,%d,%d) ",
				   arr[i].m,arr[i].n,(int)arr[i].value[0]);
    printf("\n");
  }
  
  for (i = 1; i < n; i++)
    {
      key = arr[i];
      j = i-1;

      /* Move elements of arr[0..i-1], that are
	 greater than key, to one position ahead
	 of their current position */
      while (j >= 0 && comp(arr+j, &key, order)>0)
	{
	  arr[j+1] = arr[j];
	  j = j-1;
	}
      arr[j+1] = key;
    }

    if (DEBUG) {
      printf("%d insert<: ",n);
      for (int i=0; i<n; i++) printf("(%d,%d,%d) ",
				     arr[i].m,arr[i].n,(int)arr[i].value[0]);
      printf("\n");
    }


}


// A utility function to swap two elements
void swap_b(COOB* a, COOB* b)
{
  COOB t = *a;
  *a = *b;
  *b = t;
}

/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
int partition_b (COOB *arr, int low, int high, ComparingB comp, Ordering *order)
{
  COOB pivot = arr[high]; // pivot
  int i = (low - 1); // Index of smaller element

  for (int j = low; j <= high-1; j++)
    {
      // If current element is smaller than or
      // equal to pivot
      if (comp(arr +j,&pivot,order)<=0)
	{
	  i++; // increment index of smaller element
	  swap_b(arr +i, arr +j);
	}
    }
  swap_b(arr +i + 1, arr +high);
  return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSort_b(COOB *arr, int low, int high, ComparingB comp, Ordering *order)
{
  if (low < high )
    {
      if (DEBUG) {
	printf("%d %d quick >: ", low, high);
	for (int i=low; i<=high; i++) printf("(%d,%d,%d) ",
					     arr[i].m,arr[i].n,(int)arr[i].value[0]);
	printf("\n");
      }
      
      
      if (high-low <7) 
	insertionSort_b(arr+low,high-low+1, comp, order)	;
      else {
      
	/* pi is partitioning index, arr[p] is now
	   at right place */
	int pi = partition_b(arr, low, high,comp, order);

	if (DEBUG) {
	  printf("%d pivot  quick >: ", pi);
	  printf("(%d,%d,%d) ",arr[pi].m,arr[pi].n,(int)arr[pi].value[0]);
	  printf("\n");
	}

	// Separately sort elements before
	// partition and after partition
	quickSort_b(arr, low, pi - 1,comp,order);
	quickSort_b(arr, pi + 1, high,comp,order);
      }

      if (DEBUG) {
	printf("< %d %d quick: ", low, high);
	for (int i=low; i<=high; i++) printf("(%d,%d,%d) ",arr[i].m,arr[i].n,(int)arr[i].value[0]);
	printf("\n");
      }
      
    }
  
}

/***
 *  Using coordinate elements in sparse matrices, the order is row
 *  major then column major (row, col) ~ (row*N + col). To transpose
 *  the matrix we need to change the order so that (col, row) ~ (row +
 *  col*M).
 *  
 *  here we can see two example how the transposing is done by
 *  sorting.
 */     


int roworder_b(COOB *a, COOB *b, Ordering *order) {

  return (a->m*order->N + a->n)-(b->m*order->N+b->n);
  
}

int colorder_b(COOB *a, COOB *b, Ordering *order) {

  return (a->m+ order->M *a->n)-(b->m + order->M*b->n);
  
}

void columnsort_b(COOMB *M) {
  Ordering order = {M->M, M->N, 1 } ;
  quickSort_b(M->data,0, M->length-1, colorder_b, &order);
}

void rowsort_b(COOMB *M) {
  Ordering order = {M->M, M->N, 1 } ;
  quickSort_b(M->data,0, M->length-1, roworder_b, &order);
}
