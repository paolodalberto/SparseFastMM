#include <Python.h>
#include "arrayobject.h"
#include <stdlib.h>
#define ALGEBRA_PATH 1

static int DEBUG = 0;
static int DEBUG2=0;
#include <SparseBLAS.h>
#include <parsparsecoo.h>


#define GETTIME



#ifdef GETTIME
#include <sys/time.h>
struct timeval _t1,_t2;
double duration;

#define START_CLOCK   gettimeofday(&_t1,NULL ); 
#define END_CLOCK   gettimeofday(&_t2,NULL);   duration = (_t2.tv_sec-_t1.tv_sec)+ (double)(_t2.tv_usec-_t1.tv_usec)/1000000;    printf("----------> get time %e sec<------\n",duration); 
#endif /*  GETTIME */


#ifdef CLOCK
#include <time.h>
clock_t _t1,_t2;
double duration;

#define START_CLOCK   _t1 = clock(); 
#define END_CLOCK     _t2 = clock(); duration =  ((double)(_t2-_t1))/CLOCKS_PER_SEC; \
  printf("clock time %e s \n", duration); 
#endif


PyObject *pcoomul(PyObject *self, PyObject *args); 

/**
 * An array of functions made available through the module.
 */
static PyMethodDef methods[] = {
    { "pcoomul", pcoomul, METH_VARARGS, "Sparse COO matrix R = C + A*B." },
    { NULL, NULL, 0, NULL }
};

/**
 * Initialization function for the module.
 */
void initsparsecoo (void)
{
    Py_InitModule ("sparsecoo", methods);
    import_array();
}



PyObject *pcoomul(PyObject *self, PyObject *args) {
  
  PyArrayObject *CX,*CY,*CV; int CM, CN;
  PyArrayObject *AX,*AY,*AV; int AM, AN;
  PyArrayObject *BX,*BY,*BV; int BM, BN;
  PyArrayObject *c,*r, *v; 
  int P;
  long unsigned int dim[1];
  COO R;
  PyObject *result;
  long unsigned int LC,LA,LB;
  COOE *_C;
  COOE *_A;
  COOE *_B;
  
  printf("Parsing basic iO!O!O!ii\n");
  
  if (0 &&! PyArg_ParseTuple (args, "iO!O!O!ii",
			  &P,&PyArray_Type,&CX,&PyArray_Type,&CY,&PyArray_Type,&CV,&CM,&CN
			  )) {
    printf("Argh\n");
    return NULL;
  }
  printf("Parsing basic iO!O!O!iiO!O!O!iiO!O!O!ii\n");
  if ( ! PyArg_ParseTuple (args, "iO!O!O!iiO!O!O!iiO!O!O!ii",
			  &P,
			  &PyArray_Type,&CX,&PyArray_Type,&CY,&PyArray_Type,&CV,&CM,&CN,
			  &PyArray_Type,&AX,&PyArray_Type,&AY,&PyArray_Type,&AV,&AM,&AN,
			  &PyArray_Type,&BX,&PyArray_Type,&BY,&PyArray_Type,&BV,&BM,&BN
			  )) {
    printf("Argh\n");
    return NULL;
  }


  printf("Calling C basic\n");
  printf("C X %d Y %d V %d\n",CX->dimensions[0],CY->dimensions[0],CV->dimensions[0]);
  printf("C CM %d CN %d\n",CM,CN);
  printf("A X %d Y %d V %d\n",AX->dimensions[0],AY->dimensions[0],AV->dimensions[0]);
  printf("A AM %d AN %d\n",AM,AN);
  printf("B X %d Y %d V %d\n",BX->dimensions[0],BY->dimensions[0],BV->dimensions[0]);
  printf("B BM %d BN %d\n",BM,BN);
  LC = (long unsigned int) CX->dimensions[0];
  LA = (long unsigned int) AX->dimensions[0];
  LB = (long unsigned int) BX->dimensions[0];
#ifdef GRAPH_PATH
  printf("C  %d  %d %d \n",((int*)CX->data)[3],((int*)CY->data)[3],((Mat *)CV->data)[3]);
#else
  printf("C  %d  %d %d %e \n",((int*)CX->data)[3],((int*)CY->data)[3],sizeof(Mat),((Mat *)CV->data)[3]);
#endif
  printf("Clocking\n");

  _C = from_three_to_one((int *)CX->data,(int *)CY->data,(Mat *)CV->data,LC);
  _A = from_three_to_one((int *)AX->data,(int *)AY->data,(Mat *)AV->data,LA);;
  _B = from_three_to_one((int *)BX->data,(int *)BY->data,(Mat *)BV->data,LB);; 
  
  
  START_CLOCK;
  {
    COO C = { _C, LC, CM, CN};
    COO A = { _A, LA, AM, AN};
    COO B = { _B, LB, BM, BN};

    if (!validate(C)) {
      printf("Problems with C\n");
    } 
    if (!validate(A)) {
      printf("Problems with A\n");
    }
    if (!validateT(B)) {
      printf("Problems with B\n");
    }
    
    R = matmul_coo_par(C,A,B,P);
    if (!validate(R)) {
      printf("Problems with R\n");
    } 
    
  }
  END_CLOCK;

  if (DEBUG2) printf("Done matmul_coo_par\n");
  free(_C);
  free(_A);
  free(_B);
  
  dim[0]  = R.length;
  printf("C  Computed  \n");
      
  /* Make a new float vector of same dimension */
  //PyArray_ZEROS(int nd, npy_intp* dims, int type_num, int fortran)¶
  
  c=(PyArrayObject *) PyArray_ZEROS(1,dim,NPY_INT,0);
  printf("c intermediate %d %d \n", PyArray_Size(c),c->data[0]);
  r=(PyArrayObject *) PyArray_ZEROS(1,dim,NPY_INT,0);
#ifdef GRAPH_PATH
  v=(PyArrayObject *) PyArray_ZEROS(1,dim,NPY_INT,0);
#else
  v=(PyArrayObject *) PyArray_ZEROS(1,dim,NPY_DOUBLE,0);
#endif

  printf("setting the elements\n");
  for (long unsigned int i=0; i<R.length; i++) {
    if (R.data[i].n<0 || R.data[i].m<0) {
      printf("Negative %lu %d %d\n", i,R.data[i],R.data[i]); 

    }
    ((int*)c->data)[i] = R.data[i].n;
    ((int*)r->data)[i] = R.data[i].m;
    ((Mat*)v->data)[i] = R.data[i].value;
  }

  free(R.data);
  printf("setting the output\n");
  result = PyTuple_New (3);
  printf("setting 0 \n");
  PyTuple_SetItem (result, 0, r);
  printf("setting 1 \n");
  PyTuple_SetItem (result, 1, c);
  PyTuple_SetItem (result, 2, v);
  
  printf("return  \n");
  return result;
  
}
