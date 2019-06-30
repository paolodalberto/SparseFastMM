#include "matGen.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void gen_M(int n, int sparse_factor){

	int **A,*B, *C, *D;
	int i, j, k, actualVal=0, first, same=1;
	int *val;
	int *col, *rowPtr;
	clock_t t;

//Allocating memory for matrix & vectors

	A=(int **)malloc(n*sizeof(int *));
	for(i=0; i<n; i++)
		A[i]=(int *)malloc(n*sizeof(int));

	rowPtr=(int *)malloc(n*sizeof(int));

	srand(time(0));

//Generating random matrix and sparsifying it

	printf("Generating matrix...\n");

	for(i=0;i<n;i++){
		rowPtr[i]=-1;

		for(j=0;j<n;j++){
			A[i][j]=rand()%1000;
			//printf("A[%d][%d] = %d || ", i, j, A[i][j]);
		}
		//printf("\n----------------------\n");
	}

	printf("\nSparsifying matrix...\n");

//	printf("%f\n%f\n",(double)actualVal/(n*n),(double)sparse_factor/100);

	while(((double)actualVal/(n*n))<((double)sparse_factor/100)){
		i=rand()%n;
		j=rand()%n;
		printf("%d ",actualVal);
		if(A[i][j]!=0){
			A[i][j]=0;
			actualVal++;
		}
	}

	printf("Non-sparse values: %d/%d",actualVal,n*n);

//Generating column, value and row pointer vectors

	printf("\nGenerating sparse vectors...\n");
	t=clock();

	val=(int *)malloc(actualVal*sizeof(int));
	col=(int *)malloc(actualVal*sizeof(int));
	k=0;
	for(i=0;i<n;i++){
		first=1;
		for(j=0;j<n;j++){
			if(A[i][j]!=0){
				val[k]=A[i][j];
				col[k]=j;
				if(first!=0){
					rowPtr[i]=k;
					first=0;
				}
				k++;
			}
		}
	}

	t=clock()-t;
	double time_taken_sparse_vectors = ((double)t)/CLOCKS_PER_SEC;
	printf("\tTime for sparse vectors generation: %f s\n", time_taken_sparse_vectors);

	write_M(A,n,sparse_factor);
	write_CRF_M(val, col, rowPtr, n, k, sparse_factor);

	for(i=0; i<n; i++)
		free(A[i]);
	free(A);
}

void write_M(int **M, int n, int sparse_factor){

	int i,j,x;
	int bufSize=40;
	FILE *f;
	char fname[bufSize];
	x=snprintf(fname, bufSize, "../src/data/Mat%dx%d_%d.txt",n,n,sparse_factor);
	printf("x=%d\n",x);
	printf("%s\n",fname);
	f=fopen(fname,"w+");
	if(f!=NULL){
		for(i=0; i<n; i++){
			for(j=0; j<n; j++){
				fprintf(f,"%d ", M[i][j]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}else{
		printf("Error opening file!\n");
	}
}

void write_CRF_M(int *val, int *col, int *rowPtr, int n, int k, int sparse_factor){

	int i,j,x;
	int bufSize=40;
	FILE *f;
	char fname[bufSize];
	x=snprintf(fname, bufSize, "../src/data/CRF%dx%d_%d.txt",n,n,sparse_factor);
	printf("x=%d\n",x);
	printf("%s\n",fname);
	f=fopen(fname,"w+");
	if(f!=NULL){
		fprintf(f, "%d\n",k);
		for(i=0; i<k; i++){
		fprintf(f, "%d %d\n",val[i], col[i]);
		}
		fprintf(f,"\n");
		for(i=0; i<n; i++){
			fprintf(f, "%d ",rowPtr[i]);
		}
		fclose(f);
	}else{
		printf("Error opening file!\n");
	}
}

void read_M(int **res, int n, int sparse_factor){

	int i,j, x;
	int bufSize=40;
	FILE *f;
	char fname[bufSize];
	x=snprintf(fname, bufSize, "../src/data/Mat%dx%d_%d.txt",n,n,sparse_factor);
	printf("x=%d\n",x);
	printf("%s\n",fname);
	f=fopen(fname,"r+");

	//res=(int **)malloc(n*sizeof(int *));
	//for(i=0; i<n; i++)
	//	res[i]=(int *)malloc(n*sizeof(int*));

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			fscanf(f, "%d ",&res[i][j]);
		}
	}

	fclose(f);
}

void read_CRF_M(int *val, int *col, int *rowPtr, int n, int k, int sparse_factor){

	int i,j,x;
	int bufSize=40;
	FILE *f;
	char fname[bufSize];
	x=snprintf(fname, bufSize, "../src/data/CRF%dx%d_%d.txt",n,n,sparse_factor);
	printf("x=%d\n",x);
	printf("%s\n",fname);
	f=fopen(fname,"r+");
	if(f!=NULL){
		for(i=0; i<k; i++){
			fscanf(f, "%d",&val[i]);
			fscanf(f, "%d",&col[i]);
			//printf("val & col: %d %d\n",val[i],col[i]);
		}

		for(i=0; i<n; i++){
			fscanf(f, "%d",&rowPtr[i]);
		}
		fclose(f);
	}else{
		printf("Error opening file!\n");
	}
}

void read_K(int *k, int n, int sparse_factor){
	int i,j,x;
	int bufSize=40;
	FILE *f;
	char fname[bufSize];
	x=snprintf(fname, bufSize, "../src/data/CRF%dx%d_%d.txt",n,n,sparse_factor);
	printf("x=%d\n",x);
	printf("%s\n",fname);
	f=fopen(fname,"r+");
	if(f!=NULL){
		fscanf(f, "%d\n",k);
		//printf("k is %d\n",*k);
		fclose(f);
	}else{
		printf("Error opening file!\n");
	}
}
