/* ************************************************************************
 * Copyright (c) 2018-2020 Advanced Micro Devices, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ************************************************************************ */

#include <cstdlib>
#include <iostream>
#include <rocalution.hpp>

using namespace rocalution;


void regular(LocalMatrix<double> &mat,
	     int times) {

  // rocALUTION objects
  LocalVector<double> x;
  LocalVector<double> rhs;

  // Allocate vectors
  x.Allocate(  "x"  , mat.GetM());
  rhs.Allocate("rhs", mat.GetN());

  // Convert matrix to ELL format
  // mat.ConvertToCOO();

  double time = rocalution_time();
  mat.MoveToAccelerator();
  x.MoveToAccelerator();
  rhs.MoveToAccelerator();
  time = rocalution_time() - time;
  std::cout << "Comm time:" << time << std::endl;

  // Print matrix info
  //mat.Info();
  // Print vector info
  //x.Info();
  //rhs.Info();
  
  // Set rhs to 1
  rhs.Ones();
    
  time = rocalution_time();
  // x = mat * rhs
  for (int i=0; i<times; i++)
    mat.Apply(rhs, &x);
  time = rocalution_time() - time;
  double gflops = (((double)2)*mat.GetNnz()*times/ (time/1000000))/1000000000;
    
  // Print dot product <x, rhs>
  std::cout << "Compute time:" << time <<  \
    " nnz:" << mat.GetNnz()  << std::endl;
  std::cout << "GFLOPS:" << gflops << std::endl;

  x.Clear();
  rhs.Clear();
    

}



void split(LocalMatrix<double> &A,
	   LocalMatrix<double> &mat,
	   LocalMatrix<double> &mat_d,
	   LocalMatrix<double> &mat_h,
	   int times,
	   int a, int b, int c) {

  // rocALUTION objects
  LocalVector<double> x1;
  LocalVector<double> rhs1;
  // rocALUTION objects
  LocalVector<double> x2;
  LocalVector<double> rhs2;
  // rocALUTION objects
  LocalVector<double> x3;
  LocalVector<double> rhs3;

  // Allocate vectors
  x1.Allocate(  "x1"  , mat.GetM());
  rhs1.Allocate("rhs1", mat.GetN());
  if (b) {
    // Allocate vectors
    x2.Allocate(  "x2"  , mat_d.GetM());
    rhs2.Allocate("rhs2", mat_d.GetN());
  }
  // Allocate vectors
  if (c) {
    x3.Allocate(  "x3"  , mat_h.GetM());
    rhs3.Allocate("rhs3", mat_h.GetN());
  }
  // Convert matrix to ELL format
  // mat.ConvertToCOO();

  double time = rocalution_time();
  mat.MoveToAccelerator();
  x1.MoveToAccelerator();
  rhs1.MoveToAccelerator();
  if (b) {
    mat_d.MoveToAccelerator();
    x2.MoveToAccelerator();
    rhs2.MoveToAccelerator();
    rhs2.Ones();
  }
  if (c) {
    mat_h.MoveToAccelerator();
    x3.MoveToAccelerator();
    rhs3.MoveToAccelerator();
    rhs3.Ones();
  }
  time = rocalution_time() - time;
  std::cout << "Comm time:" << time << std::endl;
  
  
  // Set rhs to 1
  rhs1.Ones();
  
    
  time = rocalution_time();
  // x = mat * rhs
  for (int i=0; i<times; i++) {
    mat.Apply(rhs1, &x1);
    if (b) mat_d.Apply(rhs2, &x2);
    if (c) mat_h.Apply(rhs3, &x3);
  }
  time = rocalution_time() - time;
  double gflops = (((double)2)*A.GetNnz()*times/ (time/1000000))/1000000000;
    
  // Print dot product <x, rhs>
  std::cout << "Compute time:" << time <<  \
    " nnz:" << A.GetNnz() <<  std::endl;
  std::cout << "GFLOPS:" << gflops << std::endl;


  x1.Clear();  x2.Clear();  x3.Clear();
  rhs1.Clear();rhs2.Clear();rhs3.Clear();
}



int main(int argc, char* argv[])
{
    // Check command line parameters
    if(argc == 1)
    {
      std::cerr << argv[0] << " device trails <matrix> <matrix_s> <matrix_d> <matrix_h> " << std::endl;
        exit(1);
    }

    set_device_rocalution(atoi(argv[1]));
    int times = atoi(argv[2]);

    // Initialize rocALUTION
    init_rocalution();
    
    

    // Print rocALUTION info
    info_rocalution();



    LocalMatrix<double> mat;
    LocalMatrix<double> mat_s;
    LocalMatrix<double> mat_d;
    LocalMatrix<double> mat_h;
    int a=1;
    int b=0;
    int c=0;
    
    // Read matrix from MTX file
    mat.ReadFileMTX(  std::string(argv[3]));
    mat_s.ReadFileMTX(std::string(argv[4])); 
    try {
      std::string a(argv[5]);
      std::cout << a << std::endl;
      if (a != "NONE")  {
	b = 1;
	mat_d.ReadFileMTX(a);
      } 
    } catch (int e) {
      printf("D is wrong \n");
    }
    try {
      std::string a(argv[6]);
      std::cout << a << std::endl;
      if (a != "NONE")  {
	c = 1;
	mat_h.ReadFileMTX(a);
      }
    } catch (int e) {
      printf("H is wrong \n");
    }

    printf("COO \n");
    mat.ConvertToCOO();
    regular(mat,times);
    mat.ConvertToCSR();
    printf("CSR \n");
    regular(mat,times);
    if ((long)mat.GetM()*mat.GetN() < ((long)8)*1024*1024*1024) {
      printf("DENSE \n");
      mat.ConvertToDENSE();
      regular(mat,times);
    }
    printf("split # %d \n", times);
    mat_s.Info();
    mat_d.Info();
    mat_h.Info();
    if (0 & b) mat_d.ConvertToDENSE();
    if (0 & c) mat_h.ConvertToDENSE();
    mat.ConvertToCOO();
    mat_s.ConvertToCOO();
    printf("COO \n");split(mat,mat_s,mat_d,mat_h,times,a,b,c);
    mat.ConvertToCSR();
    mat_s.ConvertToCSR();
    printf("CSR \n");split(mat,mat_s,mat_d,mat_h,times,a,b,c);

    
    // Stop rocALUTION platform
    stop_rocalution();

    return 0;
}
