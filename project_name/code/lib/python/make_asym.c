#include<Python.h>
#include "numpy/arrayobject.h"
#include <stdio.h>
#include <stdlib.h>

//#define IND(a,i) *((double *)(a->data+i*a->strides[0]))
#define IND(a,i,j) *((double *)(a->data+i*a->strides[0]+j*a->strides[1]))

static PyObject *make_asym(PyObject *self, PyObject *args);

static PyObject *make_asym(PyObject *self, PyObject *args)
{
  PyArrayObject *data, *dis, *weights;
  int w_truth;

  if(!PyArg_ParseTuple(args,"OOOi",&data,&dis,&weights,&w_truth)){
    return NULL;
  }

  if (!PyArray_SAMESHAPE(data,dis)){
    PyErr_Format(PyExc_ValueError,
		 "Shape of data array not equal to that of distance");
    return NULL;
  }

  int dis_i,dis_j,dis_ii,dis_jj,dis_dim0,dis_dim1;
  double r,weight_sum,sum_weight_mean,mu,core,asym,var;
  asym = 0;

  dis_dim0 = dis->dimensions[0];
  dis_dim1 = dis->dimensions[1];


  for(dis_i=0;dis_i<dis_dim0;dis_i++){
    for(dis_j=0;dis_j<dis_dim1;dis_j++){
      r = IND(dis,dis_i,dis_j);
      weight_sum = 0;
      sum_weight_mean = 0;
      core = 0;
      mu =0;
      var=0;

      for(dis_ii=0;dis_ii<dis_dim0;dis_ii++){
	for(dis_jj=0;dis_jj<dis_dim1;dis_jj++){

	  if(IND(dis,dis_ii,dis_jj)==r){
	    weight_sum += IND(weights,dis_ii,dis_jj);
	    sum_weight_mean += IND(weights,dis_ii,dis_jj)*
	      IND(data,dis_ii,dis_jj);
	  }
	}
      }
      mu = sum_weight_mean/weight_sum;
      for(dis_ii=0;dis_ii<dis_dim0;dis_ii++){
	for(dis_jj=0;dis_jj<dis_dim1;dis_jj++){

	  if(IND(dis,dis_ii,dis_jj)==r){
	    core += IND(weights,dis_ii,dis_jj)*(IND(data,dis_ii,dis_jj)-mu)
	      *(IND(data,dis_ii,dis_jj)-mu);
	      }
	}
      }
      var = core/weight_sum;
      if(w_truth == 1){
	var = var*var;
      }
      asym += var;
	}
  }
  return Py_BuildValue("d",asym);
}

static PyMethodDef make_asym_methods[] = {
  {"make_asym",make_asym,METH_VARARGS},{NULL}};

void initmake_asym()
{
  Py_InitModule("make_asym",make_asym_methods);
  import_array();
}
