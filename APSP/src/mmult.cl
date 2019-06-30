/**********
Copyright (c) 2018, Xilinx, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********/

//This function represents an OpenCL kernel. The kernel will be call from
//host application. The pointers in kernel parameters with the global 
//keyword represents cl_mem objects on the FPGA DDR memory. Array partitioning
//and loop unrolling is done to achieve better performance.

#define MAX_SIZE 4
#define BSIZE 2

#define add(a,b) (((a)<(b))?(a):(b))
#define e_a INT_MAX
#define mul(a,b) ((a)+(b))
#define e_m 0

kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void rKleene_mmultA( __global int* inA,  //Read-only input matrix A
  		      __global int* inB,  //Read-only input matrix B
  		      __global int* inC,  //Read-only input matrix C
  		      __global int* inD,  //Read-only input matrix D
            	      __global int* outB,  //Output matrix B
            	      __global int* outC,  //Output matrix C
            	      __global int* outD,  //Output matrix D
            	      int dim             //One dimension of the matrix
          	    )
{
    //Local memory to store input matrices
    __local int local_inA[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));
    __local int local_inB[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_inC[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_inD[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));
 
    __local int B[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int C[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int D[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));

    __local int local_outB[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_outC[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_outD[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));  


    //Burst reads on input matrices from DDR memory
    //Burst read for matrix local_in1 and local_in2
    read_in1A: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inA[i][j] = inA[iter];
    }
    read_in2A: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inB[i][j] = inB[iter];
    }
    read_in3A: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inC[i][j] = inC[iter];
    }
    read_in4A: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inD[i][j] = inD[iter];
    }

    //Based on the functionality the number of iterations 
    //to be executed for "loop_3" must be "dim" size. 
    //But for the pipeline to happen in the "loop_2" the
    //"loop_3" must be unrolled, to unroll the size cannot be dynamic.
    //It gives better throughput with usage of additional resources. 
    loop_1A_B: for(int i = 0; i < dim; i++){
        __attribute__((xcl_pipeline_loop))
        loop_2A_B: for(int j = 0; j < dim; j++){
		B[i][j] = 0;
            __attribute__((opencl_unroll_hint))
            loop_3A_B: for(int k = 0; k < MAX_SIZE; k++){
                B[i][j] += local_inA[i][k] * local_inB[k][j];
            }
	    local_outB[i][j] = B[i][j] + local_inB[i][j];
        }
    }    

    //Based on the functionality the number of iterations 
    //to be executed for "loop_3" must be "dim" size. 
    //But for the pipeline to happen in the "loop_2" the
    //"loop_3" must be unrolled, to unroll the size cannot be dynamic.
    //It gives better throughput with usage of additional resources. 
    loop_1A_C: for(int i = 0; i < dim; i++){
        __attribute__((xcl_pipeline_loop))
        loop_2A_C: for(int j = 0; j < dim; j++){
		C[i][j] = 0;
            __attribute__((opencl_unroll_hint))
            loop_3A_C: for(int k = 0; k < MAX_SIZE; k++){
                C[i][j] += local_inC[i][k] * local_inA[k][j];
            }
	    local_outC[i][j] = C[i][j] + local_inC[i][j];
        }
    }

    //Based on the functionality the number of iterations 
    //to be executed for "loop_3" must be "dim" size. 
    //But for the pipeline to happen in the "loop_2" the
    //"loop_3" must be unrolled, to unroll the size cannot be dynamic.
    //It gives better throughput with usage of additional resources. 
    loop_1A_D: for(int i = 0; i < dim; i++){
        __attribute__((xcl_pipeline_loop))
        loop_2A_D: for(int j = 0; j < dim; j++){
		D[i][j] = 0;
            __attribute__((opencl_unroll_hint))
            loop_3A_D: for(int k = 0; k < MAX_SIZE; k++){
                D[i][j] += local_outC[i][k] * local_outB[k][j];
            }
	    local_outD[i][j] = D[i][j] + local_inD[i][j];
        }
    }      

    //Burst write from local to DDR memory
    write_out1A: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        outB[iter] = local_outB[i][j];
    }
    write_out2A: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        outC[iter] = local_outC[i][j];
    }
    write_out3A: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        outD[iter] = local_outD[i][j];
    }
}

kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void rKleene_mmultB( __global int* inA,  //Read-only input matrix A
  		      __global int* inB,  //Read-only input matrix B
  		      __global int* inC,  //Read-only input matrix C
  		      __global int* inD,  //Read-only input matrix D
            	      __global int* outB,  //Output matrix B
            	      __global int* outC,  //Output matrix C
            	      __global int* outA,  //Output matrix D
            	      int dim             //One dimension of the matrix
          	    )
{
    //Local memory to store input matrices
    __local int local_inA[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));
    __local int local_inB[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_inC[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_inD[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));
 
    __local int B[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int C[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int A[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));

    __local int local_outB[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_outC[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 
    __local int local_outA[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2)));  


    //Burst reads on input matrices from DDR memory
    //Burst read for matrix local_in1 and local_in2
    read_in1B: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inA[i][j] = inA[iter];
    }
    read_in2B: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inB[i][j] = inB[iter];
    }
    read_in3B: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inC[i][j] = inC[iter];
    }
    read_in4B: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_inD[i][j] = inD[iter];
    }

    //Based on the functionality the number of iterations 
    //to be executed for "loop_3" must be "dim" size. 
    //But for the pipeline to happen in the "loop_2" the
    //"loop_3" must be unrolled, to unroll the size cannot be dynamic.
    //It gives better throughput with usage of additional resources. 
    loop_1B_B: for(int i = 0; i < dim; i++){
        __attribute__((xcl_pipeline_loop))
        loop_2B_B: for(int j = 0; j < dim; j++){
		B[i][j] = 0;
            __attribute__((opencl_unroll_hint))
            loop_3B_B: for(int k = 0; k < MAX_SIZE; k++){
                B[i][j] += local_inB[i][k] * local_inD[k][j];
            }
	    local_outB[i][j] = B[i][j] + local_inB[i][j];
        }
    }    

    //Based on the functionality the number of iterations 
    //to be executed for "loop_3" must be "dim" size. 
    //But for the pipeline to happen in the "loop_2" the
    //"loop_3" must be unrolled, to unroll the size cannot be dynamic.
    //It gives better throughput with usage of additional resources. 
    loop_1B_C: for(int i = 0; i < dim; i++){
        __attribute__((xcl_pipeline_loop))
        loop_2B_C: for(int j = 0; j < dim; j++){
		C[i][j] = 0;
            __attribute__((opencl_unroll_hint))
            loop_3B_C: for(int k = 0; k < MAX_SIZE; k++){
                C[i][j] += local_inD[i][k] * local_inC[k][j];
            }
	    local_outC[i][j] = C[i][j] + local_inC[i][j];
        }
    }

    //Based on the functionality the number of iterations 
    //to be executed for "loop_3" must be "dim" size. 
    //But for the pipeline to happen in the "loop_2" the
    //"loop_3" must be unrolled, to unroll the size cannot be dynamic.
    //It gives better throughput with usage of additional resources. 
    loop_1B_A: for(int i = 0; i < dim; i++){
        __attribute__((xcl_pipeline_loop))
        loop_2B_A: for(int j = 0; j < dim; j++){
		A[i][j] = 0;
            __attribute__((opencl_unroll_hint))
            loop_3B_A: for(int k = 0; k < MAX_SIZE; k++){
                A[i][j] += local_outB[i][k] * local_outC[k][j];
            }
	    local_outA[i][j] = A[i][j] + local_inA[i][j];
        }
    }      

    //Burst write from local to DDR memory
    write_out1B: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        outB[iter] = local_outB[i][j];
    }
    write_out2B: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        outC[iter] = local_outC[i][j];
    }
    write_out3B: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        outA[iter] = local_outA[i][j];
    }
}

kernel __attribute__((reqd_work_group_size(1, 1, 1)))
void mmult( __global int* in1,  //Read-only input matrix1
         __global int* out,  //Output matrix
         int dim             //One dimension of the matrix
       )
{
    //Local memory to store input matrix
    __local int local_in1[MAX_SIZE][MAX_SIZE] __attribute__((xcl_array_partition(complete, 2))); 


    //Burst reads on input matrices from DDR memory
    //Burst read for matrix local_in1 and local_in2
    read_in1: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        local_in1[i][j] = in1[iter];
    }

    //Based on the functionality the number of iterations 
    //to be executed for "loop_3" must be "dim" size. 
    //But for the pipeline to happen in the "loop_2" the
    //"loop_3" must be unrolled, to unroll the size cannot be dynamic.
    //It gives better throughput with usage of additional resources. 
    loop_1: for(int i = 0; i < dim; i++){
        __attribute__((xcl_pipeline_loop))
        loop_2: for(int j = 0; j < dim; j++){
            __attribute__((opencl_unroll_hint))
            loop_3: for(int k = 0; k < MAX_SIZE; k++){
		if(local_in1[j][k] > local_in1[j][i] + local_in1[i][k])
	            local_in1[j][k] = local_in1[j][i] + local_in1[i][k];
            }
        }
    }    

    //Burst write from local to DDR memory
    write_out: for(int iter = 0, i = 0, j = 0; iter < dim * dim; iter++, j++){
        if(j == dim){ j = 0; i++; }
        out[iter] = local_in1[i][j];
    }
}
