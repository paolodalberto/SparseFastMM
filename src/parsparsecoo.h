#ifndef PARALLEL_PARSE_COO
#define PARALLEL_PARSE_COO 1

extern COO matmul_coo_par(COO C,COO A,COO B,
			  int Ps /* number of threads */
			  );
extern COO matmul_coo_par_basic(
				int *C,
				long unsigned LC,
				int MC,
				int NC,
				int *A,
				long unsigned LA,
				int MA,
				int NA,
				int *B,
				long unsigned LB,
				int MB,
				int NB,
				int Ps /* number of threads */
				);

#endif
