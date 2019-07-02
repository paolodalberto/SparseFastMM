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
}



PyObject *pcoomul(PyObject *self, PyObject *args) {
  
  PyArrayObject *CX,*CY,*CV; int CM, CN;
  PyArrayObject *AX,*AY,*AV; int AM, AN;
  PyArrayObject *BX,*BY,*BV; int BM, BN;
  PyArrayObject *c,*r, *v; 
  int P;
  int dim[2];
  COO R;
  PyObject *result;

  if (! PyArg_ParseTuple (args, "iO!O!O!iiO!O!O!iiO!O!O!ii",
			  &P,
			  &PyArray_Type,&CX,&PyArray_Type,&CY,&PyArray_Type,&CV,
			  &CM,&CN,
			  &PyArray_Type,&AX,&PyArray_Type,&AY,&PyArray_Type,&AV,
			  &AM,&AN,
			  &PyArray_Type,&BX,&PyArray_Type,&BY,&PyArray_Type,&BV,
			  &BM,&BN
			  ))
    return NULL;

      
  R  = matmul_coo_par_basic(
			    (int*)CX->data,(int*)CY->data,(Mat*)CV->data,
			    CX->dimensions[1],
			    CM,CN,
			    (int*)AX->data,(int*)AY->data,(Mat*)AV->data,
			    AX->dimensions[1],
			    AM,AN,
			    (int*)BX->data,(int*)BY->data,(Mat*)BV->data,
			    BX->dimensions[1],
			    BM,BN,
			    P);
  
  dim[0] = 3; dim[1] = R.length;
      
      
  /* Make a new float vector of same dimension */
  c=(PyArrayObject *) PyArray_FromDims(1,dim[1],NPY_INT);
  r=(PyArrayObject *) PyArray_FromDims(1,dim[1],NPY_INT);
#ifdef GRAPH_PATH
  v=(PyArrayObject *) PyArray_FromDims(1,dim[1],NPY_INT);
#else
  v=(PyArrayObject *) PyArray_FromDims(1,dim[1],NPY_FLOAT);
#endif

  for (long unsigned int i; i<R.length; i++) {
    c->data[i] = R.data[i].n;
    r->data[i] = R.data[i].m;
    v->data[i] = R.data[i].value;
  }
    
  free(R.data);
  result = PyTuple_New (6);
  PyTuple_SetItem (result, 0, r);
  PyTuple_SetItem (result, 1, c);
  PyTuple_SetItem (result, 2, v);
  PyTuple_SetItem (result, 3, R.length);
  PyTuple_SetItem (result, 4, R.M);
  PyTuple_SetItem (result, 5, R.N);

  return result;
  
}
