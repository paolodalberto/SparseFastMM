#include "matGen.h"
#include "vecGen.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#define THREAD_NUM 8

int **A,*B, *C, *D, *E;
int i, j, k, n, valNum, first, same, sparse_factor;
int *val;
int *col, *rowPtr;

void *tMul(void *arg);

int main(int argc, char **argv)
{

	pthread_t *threads;
	clock_t t;

	valNum=0;
	same=1;

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

		threads = (pthread_t *) malloc(THREAD_NUM * sizeof(pthread_t));

		//gen_V(n);
		read_V(B,n);

// 		for(sparse_factor=95;sparse_factor<100;sparse_factor++){

			//gen_M(n, sparse_factor);
			read_M(A,n,sparse_factor);
			read_K(&valNum,n,sparse_factor);
			val=(int*)malloc(valNum*sizeof(int*));
			col=(int*)malloc(valNum*sizeof(int*));
			E=(int *)malloc(valNum*sizeof(int*));
			read_CRF_M(val,col,rowPtr,n,valNum,sparse_factor);

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
					for(j=rowPtr[i]; (k==n) ? (j<valNum) : (j<rowPtr[k]) ;j++){
						D[i]+=val[j]*B[col[j]];
					}
				}
			}

			t=clock()-t;
			double time_taken_sp = ((double)t)/CLOCKS_PER_SEC;

// function to multiply each val with B
			t=clock();

			for(i=0;i<n/THREAD_NUM;i++){
				int *tid;
		    tid = (int *) malloc( sizeof(int) );
		    *tid = i;
		    pthread_create( &threads[i], NULL, tMul, (void *)tid );
			}

			for ( i = 0; i < THREAD_NUM; ++i ) {
    		pthread_join( threads[i], NULL );
			}

			for(i=0;i<n;i++){
				k=i+1;
				if(rowPtr[i]>=0){
					while((rowPtr[k]<0)&&(k<n))
						k++;
					for(j=rowPtr[i]; (k==n) ? (j<valNum) : (j<rowPtr[k]) ;j++){
						D[i]+=E[j];
					}
				}
			}

			t=clock()-t;
			double time_taken_s_par = ((double)t)/CLOCKS_PER_SEC;

			printf("\tTime for naive multiply: %f s\n", time_taken);
			printf("\tTime for sparse multiply computation: %f s\n", time_taken_sp);
			printf("\tTime for parallel sparse multiply computation: %f s\n", time_taken_s_par);
			free(val);
			free(col);
			free(E);
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

void *tMul(void *arg){
	int i,tid;
	if(i*THREAD_NUM+tid<valNum){
			E[i*THREAD_NUM+tid]=val[i*THREAD_NUM+tid]*B[col[i*THREAD_NUM+tid]];
	}
}
