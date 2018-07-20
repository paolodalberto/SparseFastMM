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

/*****
* Changes from Lab 3:

To achieve better performance on the matrix multiplication the input and output 
arrays are partitioned and the innermost loop is unrolled in the kernel.
*****/

//OpenCL utility layer include
#include "xcl2.hpp"
#include "stdlib.h"
#include "limits.h"
#include <vector>
#include <fstream>

//Max Array Size
#define MAX_SIZE 19

//Array Size to access 
#define DATA_SIZE 19 

//Block Size
#define BSIZE 2

#define add(a,b) (((a)<(b))?(a):(b))
#define e_a 10000
#define mul(a,b) ((a)+(b))
#define e_m 0

uint64_t get_duration_ns (const cl::Event &event) {
    uint64_t nstimestart, nstimeend;
    event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_START,&nstimestart);
    event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_END,&nstimeend);
    return(nstimeend-nstimestart);
}

//CPU implementation of first part of Kleene Matrix Multiplication
//The inputs are of the size (DATA_SIZE x DATA_SIZE)
void rKleene_mmult_cpu1 (
    int *inA,   //Input Matrix 1
    int *inB,   //Input Matrix 2
    int *inC,   //Input Matrix 3
    int *inD,   //Input Matrix 4
    int *outB,   //Output Matrix 1
    int *outC,   //Output Matrix 2
    int *outD,   //Output Matrix 3
    int s,     //Short dimension of matrix
    int l     //Long dimension of matrix
)
{

    // Display the numbers produced:
    std::cout << "The mmult1 A results are: ";
    for (int ct = 0; ct < s*s; ct++){
	if(ct % s == 0)    std::cout << std::endl;
        std::cout <<  inA[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult1 B results are: ";
    for (int ct = 0; ct < s*l; ct++){
	if(ct % l == 0)    std::cout << std::endl;
        std::cout <<  inB[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult1 C results are: ";
    for (int ct = 0; ct < l*s; ct++){
	if(ct % s == 0)    std::cout << std::endl;
        std::cout <<  inC[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult1 D results are: ";
    for (int ct = 0; ct < l*l; ct++){
	if(ct % l == 0)    std::cout << std::endl;
        std::cout <<  inD[ct] << " ";
    }
    std::cout << std::endl;
    //Perform Matrix multiply B += A x B
    for(int i = 0; i < s; i++) {
        for(int j = 0; j < l; j++) {
	    outB[i * l + j] = e_a;
            for(int k = 0; k < s; k++) {
                outB[i * l + j] = add(outB[i * l + j], mul(inA[i * s + k], inB[k * l + j]));
            }
	    outB[i * l + j] = add(inB[i * l + j], outB[i * l + j]);
        }
    }

    //Perform Matrix multiply C += C x A
    for(int i = 0; i < l; i++) {
        for(int j = 0; j < s; j++) {
	    outC[i * s + j] = e_a;
            for(int k = 0; k < s; k++) {
                outC[i * s + j] = add(outC[i * s + j], mul(inC[i * s + k], inA[k * s + j]));
            }
	    outC[i * s + j] = add(inC[i * s + j], outC[i * s + j]);
        }
    }

    //Perform Matrix multiply D += C x B
    for(int i = 0; i < l; i++) {
        for(int j = 0; j < l; j++) {
	    outD[i * l + j] = e_a;
            for(int k = 0; k < s; k++) {
                outD[i * l + j] = add(outD[i * l + j], mul(outC[i * s + k], outB[k * l + j]));
            }
	    outD[i * l + j] = add(inD[i * l + j], outD[i * l + j]);
	}
    }
    // Display the numbers produced:
    std::cout << "The mmult1 outB results are: ";
    for (int ct = 0; ct < s*l; ct++){
	if(ct % l == 0)    std::cout << std::endl;
        std::cout <<  outB[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult1 outC results are: ";
    for (int ct = 0; ct < l*s; ct++){
	if(ct % s == 0)    std::cout << std::endl;
        std::cout <<  outC[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult1 outD results are: ";
    for (int ct = 0; ct < l*l; ct++){
	if(ct % l == 0)    std::cout << std::endl;
        std::cout <<  outD[ct] << " ";
    }
    std::cout << std::endl;
}  

//CPU implementation of second part of Kleene Matrix Multiplication
//The inputs are of the size (DATA_SIZE x DATA_SIZE)
void rKleene_mmult_cpu2 (
    int *inA,   //Input Matrix 1
    int *inB,   //Input Matrix 2
    int *inC,   //Input Matrix 3
    int *inD,   //Input Matrix 4
    int *outB,   //Output Matrix 1
    int *outC,   //Output Matrix 2
    int *outA,   //Output Matrix 3
    int s,     //Short dimension of matrix
    int l     //Long dimension of matrix
)
{


    // Display the numbers produced:
    std::cout << "The mmult2 A results are: ";
    for (int ct = 0; ct < s*s; ct++){
	if(ct % s == 0)    std::cout << std::endl;
        std::cout <<  inA[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult2 B results are: ";
    for (int ct = 0; ct < s*l; ct++){
	if(ct % l == 0)    std::cout << std::endl;
        std::cout <<  inB[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult2 C results are: ";
    for (int ct = 0; ct < l*s; ct++){
	if(ct % s == 0)    std::cout << std::endl;
        std::cout <<  inC[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult2 D results are: ";
    for (int ct = 0; ct < l*l; ct++){
	if(ct % l == 0)    std::cout << std::endl;
        std::cout <<  inD[ct] << " ";
    }
    std::cout << std::endl;
    //Perform Matrix multiply B += B x D
    for(int i = 0; i < s; i++) {
        for(int j = 0; j < l; j++) {
	    outB[i * l + j] = e_a;
            for(int k = 0; k < l; k++) {
                outB[i * l + j] = add(outB[i * l + j], mul(inB[i * l + k], inD[k * l + j]));
            }
	    outB[i * l + j] = add(inB[i * l + j], outB[i * l + j]);
	}
    }

    //Perform Matrix multiply C += D x C
    for(int i = 0; i <l; i++) {
        for(int j = 0; j < s; j++) {
	    outC[i * s + j] = e_a;
            for(int k = 0; k < l; k++) {
                outC[i * s + j] = add(outC[i * s + j], mul(inD[i * l + k], inC[k * s + j]));
            }
	    outC[i * s + j] = add(inC[i * s + j], outC[i * s + j]);
	}
    }

    //Perform Matrix multiply A += B x C
    for(int i = 0; i < s; i++) {
        for(int j = 0; j < s; j++) {
	    outA[i * s + j] = e_a;
            for(int k = 0; k < l; k++) {
                outA[i * s + j] = add(outA[i * s + j], mul(outB[i * l + k], outC[k * s + j]));
            }
	    outA[i * s + j] = add(inA[i * s + j], outA[i * s + j]);
	}
    }
    // Display the numbers produced:
    std::cout << "The mmult2 outB results are: ";
    for (int ct = 0; ct < s*l; ct++){
	if(ct % l == 0)    std::cout << std::endl;
        std::cout <<  outB[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult2 outC results are: ";
    for (int ct = 0; ct < l*s; ct++){
	if(ct % s == 0)    std::cout << std::endl;
        std::cout <<  outC[ct] << " ";
    }
    std::cout << std::endl;
    // Display the numbers produced:
    std::cout << "The mmult2 outA results are: ";
    for (int ct = 0; ct < s*s; ct++){
	if(ct % s == 0)    std::cout << std::endl;
        std::cout <<  outA[ct] << " ";
    }
    std::cout << std::endl;
}  

//CPU implementation of Floyd-Warshall
//The inputs are of the size (DATA_SIZE x DATA_SIZE)
void FW_cpu (
    int *in,   //Input Matrix 
    int *out,   //Output Matrix
    int dim     //One dimension of matrix
)
{
    //Initialize output matrix to input matrix
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
	    out[i * dim + j] = in[i * dim + j];
	}
    }
    //Perform Floyd-Warshall
    for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
            for(int k = 0; k < dim; k++) {
		if(out[j * dim + i] + out[i * dim + k] < out[j * dim + k])
		    out[j * dim + k] = out[j * dim + i] + out[i * dim + k];
            }
        }
    }
    // Display the numbers produced:
    std::cout << "The FW results are: ";
    for (int ct = 0; ct < dim*dim; ct++){
	if(ct % dim == 0)    std::cout << std::endl;
        std::cout <<  out[ct] << " ";
    }
    std::cout << std::endl;
} 

void RKleene_cpu (
    int *in,   //Input Matrix
    int *out,   //Output Matrix
    int dim     //One dimension of matrix
)
{

    //check if block is small enough for FW, else do R-Kleene.
    if(dim <= BSIZE){
	FW_cpu(in, out, dim);
    }else{
	int s = dim/2, l = dim - dim/2;
	int Aind = 0, Bind = 0, Cind = 0, Dind = 0;
	int A[(s)*(s)], B[(s)*(l)], C[(l)*(s)], D[(l)*(l)],
        Atmp[(s)*(s)], Btmp[(s)*(l)], Ctmp[(l)*(s)], Dtmp[(l)*(l)];
	//Initialize A, B, C & D
    	for(int i = 0; i < dim; i++) {
            for(int j = 0; j < dim; j++) {
		if(i < dim/2 && j < dim/2){
		    A[Aind] = in[i * dim + j];
		    Aind++;
		}
		if(i < dim/2 && j >= dim/2){
		    B[Bind] = in[i * dim + j];
		    Bind++;
		}
		if(i >= dim/2 && j < dim/2){
		    C[Cind] = in[i * dim + j];
		    Cind++;
		}
		if(i >= dim/2 && j >= dim/2){
		    D[Dind] = in[i * dim + j];
		    Dind++;
		}
	    }
    	}

	//Perform R-Kleene computations
	RKleene_cpu(A, Atmp, s);
	rKleene_mmult_cpu1 (Atmp, B, C, D, Btmp, Ctmp, Dtmp, s, l);
	RKleene_cpu(Dtmp, D, l);
	rKleene_mmult_cpu2 (Atmp, Btmp, Ctmp, D, B, C, A, s, l);

	//write to out
	Aind = 0, Bind = 0, Cind = 0, Dind = 0;
    	for(int i = 0; i < dim; i++) {
            for(int j = 0; j < dim; j++) {
		if(i < dim/2 && j < dim/2){
		    out[i * dim + j] = A[Aind];
		    Aind++;
		}
		if(i < dim/2 && j >= dim/2){
		    out[i * dim + j] = B[Bind];
		    Bind++;
		}
		if(i >= dim/2 && j < dim/2){
		    out[i * dim + j] = C[Cind];
		    Cind++;
		}
		if(i >= dim/2 && j >= dim/2){
		    out[i * dim + j] = D[Dind];
		    Dind++;
		}
	    }
    	}
    }
}

//Functionality to setup OpenCL context and trigger the Kernel
uint64_t RKleene_fpga (
    std::vector<int,aligned_allocator<int>>& source_in1,   //Input Matrix 1
    std::vector<int,aligned_allocator<int>>& source_fpga_results,    //Output Matrix
    int dim                                         //One dimension of matrix
)
{
    int size = dim, halfsize = size/2;    
    size_t matrix_size_bytes = sizeof(int) * size * size;

    cl::Event event;
    uint64_t kernel_duration = 0;

    //The get_xil_devices will return vector of Xilinx Devices 
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];

    //Creating Context and Command Queue for selected Device
    cl::Context context(device);
    cl::CommandQueue q1(context, device, CL_QUEUE_PROFILING_ENABLE);
    cl::CommandQueue q2(context, device, CL_QUEUE_PROFILING_ENABLE);
    std::string device_name = device.getInfo<CL_DEVICE_NAME>(); 

    //import_binary() command will find the OpenCL binary file created using the 
    //xocc compiler load into OpenCL Binary and return as Binaries
    //OpenCL and it can contain many functions which can be executed on the
    //device.
    std::string binaryFile = xcl::find_binary_file(device_name,"mmult");
    cl::Program::Binaries bins = xcl::import_binary_file(binaryFile);
    devices.resize(1);
    cl::Program program(context, devices, bins);

    //This call will extract a kernel out of the program we loaded in the
    //previous line. A kernel is an OpenCL function that is executed on the
    //FPGA. This function is defined in the src/mmult.cl file.
    cl::Kernel kernel1(program,"rKleene_mmultA");
    //cl::Kernel kernel2(program,"rKleene_mmultB");
    cl::Kernel kernel3(program,"mmult");

    //check if block is small enough for FW, else do R-Kleene.
    if(source_in1.size() <= BSIZE * BSIZE){
	/*********************************************
			    FW
	*********************************************/
	//These commands will allocate memory on the FPGA. The cl::Buffer
	//objects can be used to reference the memory locations on the device.
	//The cl::Buffer object cannot be referenced directly and must be passed
	//to other OpenCL functions.
	cl::Buffer buffer_in1(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,source_in1.data());    
	cl::Buffer buffer_output(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, 
            matrix_size_bytes,source_fpga_results.data());

    	//These commands will load the source_in1 and source_in2 vectors from the host
    	//application into the buffer_in1 and buffer_in2 cl::Buffer objects. The data
    	//will be be transferred from system memory over PCIe to the FPGA on-board
    	//DDR memory.
    	q1.enqueueMigrateMemObjects({buffer_in1},0/* 0 means from host*/);

    	//Set the kernel arguments
    	int narg = 0;
    	kernel3.setArg(narg++, buffer_in1);
    	kernel3.setArg(narg++, buffer_output);
    	kernel3.setArg(narg++, size);

    	//Launch the kernel
    	q1.enqueueTask(kernel3, NULL, &event);

    	//The result of the previous kernel execution will need to be retrieved in
    	//order to view the results. This call will write the data from the
    	//buffer_output cl_mem object to the source_fpga_results vector
    	q1.enqueueMigrateMemObjects({buffer_output},CL_MIGRATE_MEM_OBJECT_HOST);
    	q1.finish();

	kernel_duration += get_duration_ns(event);
    }else{
	int Aind = 0, Bind = 0, Cind = 0, Dind = 0;
	std::vector<int,aligned_allocator<int>> A, B, C, D,
         Atmp, Btmp, Ctmp, Dtmp;
	//Initialize A, B, C & D
    	for(int i = 0; i < dim; i++) {
            for(int j = 0; j < dim; j++) {
		if(i < dim/2 && j < dim/2){
		    A[Aind] = source_in1[i * dim + j];
		    Aind++;
		}
		if(i < dim/2 && j >= dim/2){
		    B[Bind] = source_in1[i * dim + j];
		    Bind++;
		}
		if(i >= dim/2 && j < dim/2){
		    C[Cind] = source_in1[i * dim + j];
		    Cind++;
		}
		if(i >= dim/2 && j >= dim/2){
		    D[Dind] = source_in1[i * dim + j];
		    Dind++;
		}
	    }
    	}
	//Perform R-Kleene computations
	//kernel_duration += RKleene_fpga(A, Atmp, dim/2);

	/*********************************************
			r-kleene_mmultA
	*********************************************/
	//These commands will allocate memory on the FPGA. The cl::Buffer
	//objects can be used to reference the memory locations on the device.
	//The cl::Buffer object cannot be referenced directly and must be passed
	//to other OpenCL functions.
	cl::Buffer buffer_inA(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,Atmp.data());    
	cl::Buffer buffer_inB(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,B.data()); 
	cl::Buffer buffer_inC(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,C.data()); 
	cl::Buffer buffer_inD(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,D.data()); 
	cl::Buffer buffer_outB(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, 
            matrix_size_bytes,Btmp.data());
	cl::Buffer buffer_outC(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, 
            matrix_size_bytes,Ctmp.data());
	cl::Buffer buffer_outD(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, 
            matrix_size_bytes,Dtmp.data());

    	//These commands will load the source_in1 and source_in2 vectors from the host
    	//application into the buffer_in1 and buffer_in2 cl::Buffer objects. The data
    	//will be be transferred from system memory over PCIe to the FPGA on-board
    	//DDR memory.
    	q2.enqueueMigrateMemObjects({buffer_inA, buffer_inB, buffer_inC, buffer_inD},0/* 0 means from host*/);

    	//Set the kernel arguments
    	int narg = 0;
    	kernel1.setArg(narg++, buffer_inA);
    	kernel1.setArg(narg++, buffer_inB);
    	kernel1.setArg(narg++, buffer_inC);
    	kernel1.setArg(narg++, buffer_inD);
    	kernel1.setArg(narg++, buffer_outB);
    	kernel1.setArg(narg++, buffer_outC);
    	kernel1.setArg(narg++, buffer_outD);
    	kernel1.setArg(narg++, halfsize);

    	//Launch the kernel
    	q2.enqueueTask(kernel1, NULL, &event);

    	//The result of the previous kernel execution will need to be retrieved in
    	//order to view the results. This call will write the data from the
    	//buffer_output cl_mem object to the source_fpga_results vector
    	q2.enqueueMigrateMemObjects({buffer_outB, buffer_outC, buffer_outD},CL_MIGRATE_MEM_OBJECT_HOST);
    	q2.finish();//?

	kernel_duration += get_duration_ns(event);

	//kernel_duration += RKleene_fpga(Dtmp, D, dim/2);

	/*********************************************
			r-kleene_mmultB
	*********************************************//*
	//These commands will allocate memory on the FPGA. The cl::Buffer
	//objects can be used to reference the memory locations on the device.
	//The cl::Buffer object cannot be referenced directly and must be passed
	//to other OpenCL functions.
	cl::Buffer buffer2_inA(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,Atmp.data());    
	cl::Buffer buffer2_inB(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,Btmp.data()); 
	cl::Buffer buffer2_inC(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,Ctmp.data()); 
	cl::Buffer buffer2_inD(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, 
            matrix_size_bytes,D.data()); 
	cl::Buffer buffer2_outB(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, 
            matrix_size_bytes,B.data());
	cl::Buffer buffer2_outC(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, 
            matrix_size_bytes,C.data());
	cl::Buffer buffer2_outA(context,CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, 
            matrix_size_bytes,A.data());

    	//These commands will load the source_in1 and source_in2 vectors from the host
    	//application into the buffer_in1 and buffer_in2 cl::Buffer objects. The data
    	//will be be transferred from system memory over PCIe to the FPGA on-board
    	//DDR memory.
    	q2.enqueueMigrateMemObjects({buffer_inA, buffer_inB, buffer_inC, buffer_inD},0*//* 0 means from host*//*);

    	//Set the kernel arguments
    	int narg2 = 0;
    	kernel2.setArg(narg2++, buffer2_inA);
    	kernel2.setArg(narg2++, buffer2_inB);
    	kernel2.setArg(narg2++, buffer2_inC);
    	kernel2.setArg(narg2++, buffer2_inD);
    	kernel2.setArg(narg2++, buffer2_outB);
    	kernel2.setArg(narg2++, buffer2_outC);
    	kernel2.setArg(narg2++, buffer2_outA);
    	kernel2.setArg(narg2++, halfsize);

    	//Launch the kernel
    	q2.enqueueTask(kernel2, NULL, &event);

    	//The result of the previous kernel execution will need to be retrieved in
    	//order to view the results. This call will write the data from the
    	//buffer_output cl_mem object to the source_fpga_results vector
    	q2.enqueueMigrateMemObjects({buffer2_outB, buffer2_outC, buffer2_outA},CL_MIGRATE_MEM_OBJECT_HOST);
    	q2.finish();

	kernel_duration += get_duration_ns(event);*/
    }
    return kernel_duration;
}

int main(int argc, char** argv)
{
    if (DATA_SIZE > MAX_SIZE) {
        std::cout << "Size is bigger than internal buffer size,"
        << " please use a size smaller than " << MAX_SIZE << "!" << std::endl;
        return EXIT_FAILURE;
    }

    //Allocate Memory in Host Memory
    int size = DATA_SIZE, count=0;    
    size_t matrix_size_bytes = sizeof(int) * size * size;

    //When creating a buffer with user pointer, under the hood user ptr is
    //used if and only if it is properly aligned (page aligned). When not 
    //aligned, runtime has no choice but to create its own host side buffer
    //that backs user ptr. This in turn implies that all operations that move
    //data to/from device incur an extra memcpy to move data to/from runtime's
    //own host buffer from/to user pointer. So it is recommended to use this 
    //allocator if user wish to Create Buffer/Memory Object to align user buffer
    //to the page boundary. It will ensure that user buffer will be used when 
    //user create Buffer/Mem Object.
    std::vector<int,aligned_allocator<int>> source_in1(matrix_size_bytes);
    std::vector<int,aligned_allocator<int>> source_fpga_results(matrix_size_bytes);
    std::vector<int,aligned_allocator<int>> source_cpu_results(matrix_size_bytes);

    std::ifstream inputFile("input.txt");        // Input file stream object

    // Check if exists and then open the file.
    /*if (inputFile.is_open()) {
        // Push items into a vector
        int current_number = 0;
        while (inputFile >> current_number && count < size*size){
	    if(current_number == 0)
		source_in1[count] = e_a;
	    else
		source_in1[count] = current_number;
	    count++;
        }

        // Close the file.
        inputFile.close();
    }else {
        std::cout << "Error! Cannot open file.\n";
        return 0;
    }*/

    //Create the test data and Software Result 
    for(int i = 0 ; i < DATA_SIZE * DATA_SIZE ; i++){
        source_in1[i] = rand() % size;
        source_cpu_results[i] = 0;
        source_fpga_results[i] = 0;
    }

        // Display the numbers read:
        std::cout << "The numbers are: ";
        for (int ct = 0; ct < size*size; ct++){
	    if(ct % size == 0)    std::cout << std::endl;
            std::cout << source_in1[ct] << " ";
        }
        std::cout << std::endl;

    uint64_t kernel_duration = 0;  

    FW_cpu(source_in1.data(), source_cpu_results.data(), size);

    //Compute CPU Results
    RKleene_cpu(source_in1.data(), source_cpu_results.data(), size);

        // Display the numbers produced:
        std::cout << "The results are: ";
        for (int ct = 0; ct < size*size; ct++){
	    if(ct % size == 0)    std::cout << std::endl;
            std::cout <<  source_cpu_results[ct] << " ";
        }

        std::cout << std::endl;

    //Compute FPGA Results
    kernel_duration = RKleene_fpga(source_in1, source_fpga_results, size);

    //Compare the results of the FPGA to CPU 
    bool match = true;
    for (int i = 0 ; i < size * size; i++){
        if (source_fpga_results[i] != source_cpu_results[i]){
            std::cout << "Error: Result mismatch" << std::endl;
            std::cout << "i = " << i << " CPU result = " << source_cpu_results[i]
                << " FPGA result = " << source_fpga_results[i] << std::endl;
            match = false;
            break;
        }
    }

    std::cout << "TEST " << (match ? "PASSED" : "FAILED") << std::endl; 

    std::cout << "Wall Clock Time (Kernel execution): " << kernel_duration << std::endl;
    std::cout << "Note: Wall Clock Time is meaningful for real hardware execution only,"  
            << "not for emulation." << std::endl; 

    return (match ? EXIT_SUCCESS :  EXIT_FAILURE);
}
