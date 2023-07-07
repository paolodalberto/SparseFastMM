typedef COO  (*MatrixComputation)(COO A, COO B);
typedef struct operands_addition TAddOperands;

struct operands_addition { 
  int  pi;
  MatrixComputation m;  // C = A*B 
  COO   *c;
  COO   a;
  COO   b;
} ;



#ifndef PARALLEL_PARSE_COO
#define PARALLEL_PARSE_COO 1


extern void MatrixComputations(TAddOperands *args, int len);
extern COO merge_alt( COO C, COO T);

extern COO matmul_coo_par(COO C,COO A,COO B,
			  int Ps /* number of threads */
			  );
extern COO matmul_coo_par_basic(
				int *CX, int *CY, Mat *CV,
				long unsigned int LC, int MC, int NC,
				int *AX, int *AY, Mat *AV,
				long unsigned int LA, int MA, int NA,
				int *BX, int *BY, Mat *BV,
				long unsigned int LB, int MB, int NB, int Ps /* number of threads */
				);

#endif
