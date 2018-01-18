#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *chisq(PyObject *self, PyObject *args);

static PyObject *chisq(PyObject *self, PyObject *args)
{
  PyArrayObject *mod, *data, *errors;
  int i,lim;
  double chi;

  if(!PyArg_ParseTuple(args,"OOO", &mod, &data, &errors))
    {
      return NULL;
    }
  
  lim = mod->dimensions[0];
  chi = 0;

  //CANNOT PERFORM PARALLEL CALCULATION HERE
  for(i=0;i<lim;i++)
    {
      chi += pow((IND(mod,i)-IND(data,i))/IND(errors,i),2);
    }
  return Py_BuildValue("d",chi);
}

static char chisq_doc[]="\
   This function creates the chi squared statistic given input parameters.\n\
\n\
   Parameters\n\
   ----------\n\
   mod:   1D NPY ARRAY - contains the model to be tested\n\
   data:  1D NPY ARRAY - contains the actual measurements\n\
   errors 1D NPY ARRAY - errors made on the meaurements (not weights)\n\
\n\
   Returns\n\
   -------\n\
   Float - the chi squared value given the model and weights\n\
\n\
   Revisions\n\
   ---------\n\
   2011-01-08    Nate Lust, UCF\n\
                 natelust at linux dot com\n\
                 Initial version, as c extension\n\
";

static PyMethodDef chisq_methods[] = {
  {"chisq", chisq,METH_VARARGS,chisq_doc},{NULL}};

void initchisq(void)
{
  Py_InitModule("chisq", chisq_methods);
  import_array();
}

