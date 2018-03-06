#ifndef PARALLEL_PARSE_COO
#define PARALLEL_PARSE_COO 1

extern COO matmul_coo_par(COO C,COO A,COO B,
			  int Ps /* number of threads */
			  );

#endif
