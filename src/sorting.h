


#ifndef SORT_DEFINITION
#define SORT_DEFINITION 1

void quickSort(COOE *arr, int low, int high, Comparing comp, Ordering *order);
extern void columnsort(COO *M);
extern void rowsort(COO *M); 
extern int roworder(COOE *a, COOE *b, Ordering *order);
extern int colorder(COOE *a, COOE *b, Ordering *order);
extern void columnsort_b(COOMB *M);
extern void rowsort_b(COOMB *M); 
extern int roworder_b(COOB *a, COOB *b, Ordering *order);
extern int colorder_b(COOB *a, COOB *b, Ordering *order);
#endif
