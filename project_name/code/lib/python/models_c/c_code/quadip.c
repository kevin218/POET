#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))
//#define IND_arr(a,i) (PyArrayObject *)(a->data+i*a->strides[0])
#define IND2(a,i,j) *((double *)(a->data+i*a->strides[0]+j*a->strides[1]))

static PyObject *quadip(PyObject *self, PyObject *args, PyObject *keywds);

static PyObject *quadip(PyObject *self, PyObject *args, PyObject *keywds)
{
  PyObject *etc;
  PyArrayObject *out, *ipparams, *position;
  double a,b,c,d,e,f;
  int i;
  npy_intp dims[1];

  static char *kwlist[] = {"ipparams","position","etc",NULL};

  if(!PyArg_ParseTupleAndKeywords(args,keywds,"OO|O",kwlist,&ipparams,&position,&etc))
    {
      return NULL;
    }

  a = IND(ipparams,0);
  b = IND(ipparams,1);
  c = IND(ipparams,2);
  d = IND(ipparams,3);
  e = IND(ipparams,4);
  f = IND(ipparams,5);

  dims[0] = PyArray_DIM(position, 1);

  out = (PyArrayObject *) PyArray_SimpleNew(1,dims,PyArray_DOUBLE);
  
  #pragma omp parallel for
  for(i=0;i<dims[0];i++)
    {
      IND(out,i) = a*pow(IND2(position,0,i),2)+b*pow(IND2(position,1,i),2)+  \
                   c*IND2(position,0,i)*IND2(position,1,i)+d*IND2(position,0,i)+e*IND2(position,1,i)+f;
    }
  return PyArray_Return(out);
}

static char quadip_doc[]="\
  This function fits the intra-pixel sensitivity effect using a 2D quadratic.\n\
\n\
  Parameters\n\
  ----------\n\
    a: quadratic coefficient in y\n\
    b: quadratic coefficient in x\n\
    c: coefficient for cross-term\n\
    d: linear coefficient in y\n\
    e: linear coefficient in x\n\
    f: constant\n\
\n\
  Returns\n\
  -------\n\
    returns the flux values for the intra-pixel model\n\
\n\
  Revisions\n\
  ---------\n\
  2008-07-05	Kevin Stevenson, UCF  \n\
			kevin218@knights.ucf.edu\n\
		Original version\n\
  2011-01-05    Nate Lust, UCF\n\
                natelust at linux dot com\n\
                Converted to c extention function\n\
\n\
";

static PyMethodDef quadip_methods[] = {
  {"quadip",(PyCFunction)quadip,METH_VARARGS|METH_KEYWORDS,quadip_doc},{NULL}};

void initquadip(void)
{
  Py_InitModule("quadip",quadip_methods);
  import_array();
}
