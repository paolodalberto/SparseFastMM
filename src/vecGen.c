#include "vecGen.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void gen_V(int n){

	int j, *B;
	//Generating dense vector

	printf("Generating dense vector... \n");

	B=(int *)malloc(n*sizeof(int*));
	for(j=0;j<n;j++){
		B[j]=rand()%1000;
	}

	write_V(B, n);
}

void write_V(int *V, int n){

	int i,x;
	int bufSize=40;
	FILE *f;
	char fname[bufSize];
	x=snprintf(fname, bufSize, "../src/data/Vec%d.txt",n);
	printf("x=%d\n",x);
	printf("%s\n",fname);
	f=fopen(fname,"w+");
	if(f!=NULL){
		for(i=0; i<n; i++){
			fprintf(f,"%d ", V[i]);
		}
		fclose(f);
	}else{
		printf("Error opening file!\n");
	}
}

void read_V(int *V, int n){

	int i,x;
	int bufSize=40;
	FILE *f;
	char fname[bufSize];
	x=snprintf(fname, bufSize, "../src/data/Vec%d.txt",n);
	printf("x=%d\n",x);
	printf("%s\n",fname);
	f=fopen(fname,"r+");
	if(f!=NULL){
		for(i=0; i<n; i++){
			fscanf(f,"%d ", &V[i]);
		}
		fclose(f);
	}else{
		printf("Error opening file!\n");
	}
}
