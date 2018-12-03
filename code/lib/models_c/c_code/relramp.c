#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *relramp(PyObject *self, PyObject *args, PyObject *keywds);

static PyObject *relramp(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *etc;
  PyArrayObject *x,*y, *rampparams;
  double goal,m,x0,a,b,x1;
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
  x1   = IND(rampparams,5);

  dims[0] = x->dimensions[0];

  y = (PyArrayObject *) PyArray_SimpleNew(1,dims,PyArray_DOUBLE);
  #pragma omp parallel for
  for(i=0;i<dims[0];i++)
    {
      IND(y,i) = goal*(1-exp(-1*m*(IND(x,i)-x0)))+a*(IND(x,i)-x1)+b;
    }
  return PyArray_Return(y);
}

static PyMethodDef module_methods[] = {
  {"relramp",(PyCFunction)relramp,METH_VARARGS|METH_KEYWORDS,module_docstring},{NULL}};

static char module_docstring[]="\
    This function NEEDS A DOC_STRING.\n\
  \n\
    Parameters\n\
    ----------\n\
  \n\
    Returns\n\
    -------\n\
  \n\
    Revisions\n\
    ---------\n\
    2010-07-30    Kevin Stevenson, UCF  \n\
                  kevin218@knights.ucf.edu\n\
                  Original version\n\n\
    2010-12-24    Nate Lust, UCF\n\
                  natelust at linux dot com\n\
                  Converted to C\n\n\
    2018-11-22    Jonathan Fraine, SSI\n\
                  jfraine at spacescience.org\n\
                  Updated c extensions to python3, with support for python2.7\n\n\
";

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
    PyInit_relramp(void)
#else
    initrelramp(void)
#endif
{
    #if PY_MAJOR_VERSION >= 3
        PyObject *module;
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "relramp",             /* m_name */
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
        PyObject *m = Py_InitModule3("relramp", module_methods, module_docstring);
        if (m == NULL)
            return;
        /* Load `numpy` functionality. */
        import_array();
    #endif
}
