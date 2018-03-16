#include "matGen.h"
#include "vecGen.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#define THREAD_NUM 8

void *tMul(void *arg) {
	/* code */
	int tid = *(int *)arg;
}

int main(int argc, char **argv)
{
	int **A,*B, *C, *D;
	int i, j, k, n, actualVal=0, first, same=1, sparse_factor;
	int *val;
	int *col, *rowPtr;

	pthread_t *threads;
	clock_t t;

	printf("Give matrix size (N):\n");
	scanf("%d",&n);
	printf("N = %d\n",n);

	printf("Give sparse percentage:\n");
	scanf("%d",&sparse_factor);
	printf("Sparse percentage: %d\n",sparse_factor);

	srand(time(0));

//	for(n=20000;n<=25000;n+=5000){
		A=(int **)malloc(n*sizeof(int *));
		for(i=0; i<n; i++)
			A[i]=(int *)malloc(n*sizeof(int*));
		B=(int *)malloc(n*sizeof(int*));
		C=(int *)malloc(n*sizeof(int*));
		D=(int *)malloc(n*sizeof(int*));
		rowPtr=(int *)malloc(n*sizeof(int*));

		threads = (pthread_t *) malloc(num_threads * sizeof(pthread_t));

		//gen_V(n);
		read_V(B,n);

// 		for(sparse_factor=95;sparse_factor<100;sparse_factor++){

			//gen_M(n, sparse_factor);
			read_M(A,n,sparse_factor);
			read_K(&k,n,sparse_factor);
			val=(int*)malloc(k*sizeof(int*));
			col=(int*)malloc(k*sizeof(int*));
			read_CRF_M(val,col,rowPtr,n,k,sparse_factor);

			printf("\nComputing naive multiplication...\n");
			t=clock();

			for(i=0;i<n;i++){
				for(j=0;j<n;j++){
					C[i]+=A[i][j]*B[j];
				}
			}

			t=clock()-t;
			double time_taken = ((double)t)/CLOCKS_PER_SEC;

			printf("\nComputing compressed row multiplication...\n");

			t=clock();

			for(i=0;i<n;i++){
				k=i+1;
				if(rowPtr[i]>=0){
					while((rowPtr[k]<0)&&(k<n))
						k++;
					for(j=rowPtr[i]; (k==n) ? (j<actualVal) : (j<rowPtr[k]) ;j++){
						D[i]+=val[j]*B[col[j]];
					}
				}
			}
// function to multiply each val with B
			for(i=0;i<n/THREAD_NUM);i++){
			//	i*thre
			}
			for ( i = 0; i < THREAD_NUM; ++i ) {
    		int *tid;
    		tid = (int *) malloc( sizeof(int) );
    		*tid = i;
    		pthread_create( &threads[i], NULL, worker, (void *)tid );
			}

			t=clock()-t;
			double time_taken_sp = ((double)t)/CLOCKS_PER_SEC;

			printf("\tTime for naive multiply: %f s\n", time_taken);
			printf("\tTime for sparse multiply computation: %f s\n", time_taken_sp);
			free(val);
			free(col);
//		}

		free(rowPtr);
		free(B);
		free(C);
		free(D);
		for(i=0; i<n; i++)
			free(A[i]);
		free(A);

//	}

	return 0;
}
