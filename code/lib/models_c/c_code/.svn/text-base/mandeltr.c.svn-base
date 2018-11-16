#include<Python.h>
#include<numpy/arrayobject.h>
#include<math.h>
#include<omp.h>

#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *mandeltr(PyObject *self, PyObject *args, PyObject *keywd);

static PyObject *mandeltr(PyObject *self, PyObject *args, PyObject *keywd)
{
  PyObject *etc;
  PyArrayObject *t, *y, *params;
  npy_intp dims[1];
  
  static char *kwlist[] = {"params","t","etc",NULL};
  
  if(!PyArg_ParseTupleAndKeywords(args,keywd,"OO|O",kwlist,	\
				  &params,&t,&etc))
    {
      return NULL;
    }
 
  double midpt,rprs,cosi,ars,flux;
  double p,z,k0,k1,mod;

  midpt = IND(params,0);
  rprs  = IND(params,1);
  cosi  = IND(params,2);
  ars   = IND(params,3);
  flux  = IND(params,4);
  p     = IND(params,5);
  
  dims[0] = t->dimensions[0];
 
  y = (PyArrayObject *) PyArray_SimpleNew(1,dims,PyArray_DOUBLE);
  
  int i;
  #pragma omp parallel for private(z,k0,k1)
  for(i=0;i<dims[0];i++)
    {
      IND(y,i) = 1;
      mod = (IND(t,i)-midpt)-floor((IND(t,i)-midpt)/p)*p;
      if((mod > p/4.) && (mod < 3*p/4.))
	{
	  z = ars;
	}
      else
	{
	  z = ars*sqrt(pow(sin(2*M_PI*(IND(t,i)-midpt)/p),2)+pow((cosi*cos(2*M_PI \
	     *(IND(t,i)-midpt)/p)),2));
	}
      if(z<=(1-rprs))
	{
	  IND(y,i) = 1-pow(rprs,2);
	}
      if(z>(1-rprs)&&z<=(1+rprs))
	{
	  k0       = acos((rprs*rprs+z*z-1)/2/rprs/z);
	  k1       = acos((1-rprs*rprs+z*z)/2/z);
	  IND(y,i) = 1-1/M_PI*(k0*rprs*rprs+k1-sqrt((4*z*z-\
				     pow((1+z*z-rprs*rprs),2))/4));
	}

      IND(y,i) *= flux;

    }
  return PyArray_Return(y);
}

static char mandeltr_doc[] ="\
   This function computes the primary transit shape using equations provided by Mandel & Agol (2002)\n\
\n\
  Parameters\n\
  ----------\n\
    midpt:  Center of eclipse\n\
    rprs:   Planet radius / stellar radius\n\
    cosi:   Cosine of the inclination\n\
    ars:    Semi-major axis / stellar radius\n\
    flux:   Flux offset from 0\n\
    t:	    Array of phase/time points\n\
    p:      Period in same units as t\n\
\n\
  Returns\n\
  -------\n\
    This function returns the flux for each point in t.\n\
\n\
  Revisions\n\
  ---------\n\
  2010-11-27	Kevin Stevenson, UCF\n\
                kevin218@knights.ucf.edu\n\
                Original version\n\
  2010-12-19    Nate Lust, UCF\n\
                natelust at linux dot com\n\
                converted function to c\n\
";

static PyMethodDef mandeltr_methods[]={
  {"mandeltr",mandeltr,METH_VARARGS|METH_KEYWORDS,mandeltr_doc}, \
  {NULL}};

void initmandeltr(void)
{
  Py_InitModule("mandeltr",mandeltr_methods);
  import_array();
}
