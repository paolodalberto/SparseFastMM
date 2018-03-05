#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

//Matrix generation function
void gen_M(int n, int sparse_factor);

//write Matrix to file in Regular Format
void write_M(int **M, int n, int sparse_factor);

//write Matrix to file in Compressed Row Format
void write_CRF_M(int *val, int *col, int *rowPtr, int n, int k, int sparse_factor);

//read Matrix from regular format file
void read_M(int **res, int n, int sparse_factor);

//read Matrix from CRF format file
void read_CRF_M(int *val, int *col, int *rowPtr, int n, int k, int sparse_factor);

//read number of elements in sparse matrix
void read_K(int *k, int n, int sparse_factor);