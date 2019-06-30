#ifndef SORT_DEFINITION
#define SORT_DEFINITION 1

extern void columnsort(COO *M);
extern void rowsort(COO *M); 
extern int roworder(COOE *a, COOE *b, Ordering *order);
extern int colorder(COOE *a, COOE *b, Ordering *order);
#endif
