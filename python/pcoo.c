#include <Python.h>
#include "arrayobject.h"

#define GRAPH_PATH 1

static int DEBUG = 0;
static int DEBUG2=1;
#include <SparseBLAS.h>
#include <parsparsecoo.h>


PyObject *pcoomul(PyObject *self, PyObject *args); 

/**
 * An array of functions made available through the module.
 */
static PyMethodDef anomalymethods[] = {
    { "pcoomul", pcoomul, METH_VARARGS, "Sparse COO matrix R = C + A*B." },
    { NULL, NULL, 0, NULL }
};

/**
 * Initialization function for the module.
 */
void initpcoo (void)
{
    Py_InitModule ("sparsecoo", anomalymethods);
}



PyObject *pcoomul(PyObject *self, PyObject *args) {
  
  PyArrayObject *C; int CM, CN;
  PyArrayObject *A; int AM, AN;
  PyArrayObject *B; int BM, BN;
  PyArrayObject *c; 
  int P;
  int dim[2];
  COO R;

  if (! PyArg_ParseTuple (args, "iO!iiO!iiO!ii",
			  &P,
			  &PyArray_Type,&C,&CM,&CN,
			  &PyArray_Type,&A,&AM,&AN,
			  &PyArray_Type,&B,&BM,&BN
			  ))
    return NULL;

      
  R  = matmul_coo_par_basic(
			    (int*)C->data,
			    C->dimensions[1],
			    CM,CN,
			    (int*)A->data,
			    A->dimensions[1],
			    AM,AN,
			    (int*)B->data,
			    B->dimensions[1],
			    BM,BN,
			    P);
  
  dim[0] = 3; dim[1] = R.length;
      
      
  /* Make a new float vector of same dimension */
  c=(PyArrayObject *) PyArray_FromDims(2,dim,NPY_INT);
  
  memcpy( (void*) c->data, (void*) R.data,
	  dim[0]*dim[1]*sizeof(int)
	  );
  free(R.data);
  
  return PyArray_Return(c);
  
}
