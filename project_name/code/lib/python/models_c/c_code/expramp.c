#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *expramp(PyObject *self, PyObject *args, PyObject *keywds);

static PyObject *expramp(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *etc;
  PyArrayObject *t,*y, *rampparams;
  double goal,m,a;
  int i;
  npy_intp dims[1];

  //  etc = PyList_New(0);

  static char *kwlist[] = {"rampparams","t","etc",NULL};

  if(!PyArg_ParseTupleAndKeywords(args,keywds,"OO|O"		\
				  ,kwlist,&rampparams,&t,&etc))
    {
      return NULL;
    }

  goal = IND(rampparams,0);
  m    = IND(rampparams,1);
  a    = IND(rampparams,2);

  dims[0] = t->dimensions[0];

  y = (PyArrayObject *) PyArray_SimpleNew(1,dims,PyArray_DOUBLE);



  #pragma omp parallel for
  for(i=0;i<dims[0];i++)
    {
      IND(y,i) = goal-a*exp(-m*IND(t,i));
    }

  return PyArray_Return(y);
}


static PyMethodDef expramp_methods[] = {
  {"expramp",(PyCFunction)expramp,METH_VARARGS|METH_KEYWORDS},{NULL}};

void initexpramp(void)
{
  Py_InitModule("expramp",expramp_methods);
  import_array();
}
