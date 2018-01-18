#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *felramp(PyObject *self, PyObject *args, PyObject *keywds);

static PyObject *felramp(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *etc;
  PyArrayObject *t,*y,*rampparams;
  double goal,m,t0,a,t1;
  int i;
  npy_intp dims[1];

  //etc = PyList_New(0);

  static char *kwlist[] = {"rampparams","t","etc",NULL};

  if(!PyArg_ParseTupleAndKeywords(args,keywds,"OO|O"		\
				  ,kwlist,&rampparams,&t,&etc))
    {
      return NULL;
    }

  goal = IND(rampparams,0);
  m    = IND(rampparams,1);
  t0   = IND(rampparams,2);
  a    = IND(rampparams,3);
  t1   = IND(rampparams,4);

  dims[0] = t->dimensions[0];

  y = (PyArrayObject *) PyArray_SimpleNew(1,dims,PyArray_DOUBLE);

  for(i=0;i<dims[0];i++)
    {
      IND(y,i) = goal*(1+exp(-1*m*(IND(t,i)-t0)))+a*(IND(t,i)-t1);
    }

  return PyArray_Return(y);
}


static PyMethodDef felramp_methods[] = {
  {"felramp",(PyCFunction)felramp,METH_VARARGS|METH_KEYWORDS},{NULL}};

void initfelramp(void)
{
  Py_InitModule("felramp",felramp_methods);
  import_array();
}
