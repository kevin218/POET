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


static PyMethodDef module_methods[] = {
  {"expramp",(PyCFunction)expramp,METH_VARARGS|METH_KEYWORDS},{NULL}};

static char module_docstring[] =
"This module is used to calcuate the expramp";

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_expramp(void)
#else
initexpramp(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
PyObject *module;
static struct PyModuleDef moduledef = {
PyModuleDef_HEAD_INIT,
"expramp",             /* m_name */
module_docstring,    /* m_doc */
-1,                  /* m_size */
module_methods,      /* m_methods */
NULL,                /* m_reload */
NULL,                /* m_traverse */
NULL,                /* m_clear */
NULL,                /* m_free */
};
#endif

#if PY_MAJOR_VERSION >= 3
module = PyModule_Create(&moduledef);
if (!module)
return NULL;
/* Load `numpy` functionality. */
import_array();
return module;
#else
PyObject *m = Py_InitModule3("expramp", module_methods, module_docstring);
if (m == NULL)
return;
/* Load `numpy` functionality. */
import_array();
#endif
}
