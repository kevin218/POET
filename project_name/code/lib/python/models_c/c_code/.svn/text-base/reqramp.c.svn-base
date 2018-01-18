#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *reqramp(PyObject *self, PyObject *args, PyObject *keywds);

static PyObject *reqramp(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *etc;
  PyArrayObject *x,*y, *rampparams;
  double goal,m,x0,a,b,c,x1;
  int i;
  npy_intp dims[1];

  //  etc = PyList_New(0);

  static char *kwlist[] = {"rampparams","x","etc",NULL};

  if(!PyArg_ParseTupleAndKeywords(args,keywds,"OO|O",kwlist,&rampparams,&x,&etc))
    {
      return NULL;
    }

  goal = IND(rampparams,0);
  m    = IND(rampparams,1);
  x0   = IND(rampparams,2);
  a    = IND(rampparams,3);
  b    = IND(rampparams,4);
  c    = IND(rampparams,5);
  x1   = IND(rampparams,6);

  dims[0] = x->dimensions[0];

  y = (PyArrayObject *) PyArray_SimpleNew(1,dims,PyArray_DOUBLE);
  #pragma omp parallel for
  for(i=0;i<dims[0];i++)
    {
      IND(y,i) = goal*(1-exp(-1*m*(IND(x,i)-x0)))+a*pow((IND(x,i)-x1),2)	\
						       +b*(IND(x,i)-x1)+c;
    }
  return PyArray_Return(y);
}

static PyMethodDef reqramp_methods[] = {
  {"reqramp",(PyCFunction)reqramp,METH_VARARGS|METH_KEYWORDS},{NULL}};

void initreqramp(void)
{
  Py_InitModule("reqramp",reqramp_methods);
  import_array();
}
