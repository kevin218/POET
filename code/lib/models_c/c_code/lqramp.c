#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *lqramp(PyObject *self, PyObject *args, PyObject *keywds);

static PyObject *lqramp(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *etc;
  PyArrayObject *x,*y, *rampparams;
  double x0,a,b,c,d,x1,temp;
  int i;
  npy_intp dims[1];

  //etc = PyList_New(0);

  static char *kwlist[] = {"rampparams","x","etc",NULL};

  if(!PyArg_ParseTupleAndKeywords(args,keywds,"OO|O",kwlist,&rampparams,&x,&etc))
    {
      return NULL;
    }

  x0 = IND(rampparams,0);
  a  = IND(rampparams,1);
  b  = IND(rampparams,2);
  c  = IND(rampparams,3);
  d  = IND(rampparams,4);
  x1 = IND(rampparams,5);

  dims[0] = x->dimensions[0];

  y = (PyArrayObject *) PyArray_SimpleNew(1,dims,PyArray_DOUBLE);

  for(i=0;i<dims[0];i++)
    {
      if(IND(x,i)<=x0)
        IND(y,i) = 1;
      else
        IND(y,i) = a*log(IND(x,i)-x0)+b*pow((IND(x,i)-x1),2)+c*(IND(x,i)-x1)+d;
    }
  return PyArray_Return(y);
}

static char module_docstring[]="\
  NAME:\n\
	LQRAMP\n\
\n\
 PURPOSE:\n\
	This function creates a model that fits a log + quadratic ramped eclipse\n\
\n\
 CATEGORY:\n\
	Astronomy.\n\
\n\
 CALLING SEQUENCE:\n\
\n\
	Result = LQRAMP([a,b,c,d],x)\n\
\n\
 INPUTS:\n\
    x0: phase offset\n\
	a:	log(x) term\n\
	b:	quadratic term\n\
	c:	linear term\n\
    d:  constant term\n\
	x:	Array of time/phase points\n\
\n\
 OUTPUTS:\n\
	This function returns an array of y values...\n\
\n\
 EXAMPLE:\n\
\n\
\n\
 MODIFICATION HISTORY:\n\
 	Written by:	Kevin Stevenson, UCF  	2009-11-28\n\
			kevin218@knights.ucf.edu\n\

	2018-11-22  Jonathan Fraine, SSI\n\
                jfraine at spacescience.org\n\
                Updated c extensions to python3, with support for python2.7\n\
";

static PyMethodDef module_methods[] = {
  {"lqramp",(PyCFunction)lqramp,METH_VARARGS|METH_KEYWORDS,module_docstring},{NULL}};

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
    PyInit_lqramp(void)
#else
    initlqramp(void)
#endif
{
    #if PY_MAJOR_VERSION >= 3
        PyObject *module;
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "lqramp",             /* m_name */
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
        PyObject *m = Py_InitModule3("lqramp", module_methods, module_docstring);
        if (m == NULL)
            return;
        /* Load `numpy` functionality. */
        import_array();
    #endif
}
