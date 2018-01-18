/*
	This code computes the separation of centers between the occulting and occulted objects.  
	It takes as arguments an array of times t,
	eccentricity e, semi-major axis a (in units of stellar radii), inclination angle i,
	stellar radius r_s, longitude of periastron w, orbital period P, time of periastron t0, and
	an error tolerance eps for computing the eccentric anomaly.

	LK 9/14/12 
*/
#include <Python.h>
#include<numpy/arrayobject.h>

#define TWOPI 6.28318531
#define PI 3.14159265
#define IND(a,i) *((double *)(a->data+i*a->strides[0]))

static PyObject *rsky(PyObject *self, PyObject *args);

double getE(double M, double e)	//determines the eccentric anomaly (Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
{
	double E = M, eps = 1.0e-7;
	
	while((E - e*sin(E) - M) > eps)
	{
		E = E - (E - e*sin(E) - M)/(1.0 - e*cos(E));
	}
	return E;
}

static PyObject *rsky(PyObject *self, PyObject *args)
{
	double e, E, i, a, r, r_s, d, f, w, P, M, n, t0, eps;
	int ii;
	npy_intp dims[1];
	PyArrayObject *z, *t;
  	if(!PyArg_ParseTuple(args,"ddddddddO", &e, &a, &i, &r_s, &w, &P, &t0, &eps, &t))
   	 {
      		return NULL;
   	 }
	dims[0] = t->dimensions[0];

	z = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_DOUBLE);

	n = TWOPI/P;	// mean motion
	#pragma omp parallel for private(M,E,r,f,d)
	for(ii = 0; ii < dims[0]; ii++)
	{
//		t = t0 + 0.005*P*(double)ii/(double)npoints;
		if(e > eps)
		{
			M = n*(IND(t,ii) - t0);
			E = getE(M, e);
			r = a*(1.0 - e*cos(E));
			f = acos(a*(1.0 - e*e)/(r*e) - 1.0/e);
		}
		else f = ((IND(t, ii)-t0)/P - (int)((IND(t,ii)-t0)/P))*TWOPI;
		d = a*(1.0-e*e)/(1.0+e*cos(f))*sqrt(1.0 - sin(w+f)*sin(w+f)*sin(i)*sin(i));
		IND(z, ii) = d/r_s;
	}
	return PyArray_Return(z);
} 

static char rsky_doc[] = "\
This code computes the distance between the centers of the\n\
star and the planet in the plane of the sky.  This parameter is\n\
denoted r_sky = sqrt(x^2 + y^2) in the Seager Exoplanets book\n\
(see the section by Murray, and Winn eq. 5).  In the Mandel & Agol (2002) paper,\n\
this quantity is denoted d.\n\
K 4/27/12 ";

static PyMethodDef rsky_methods[] = {
  {"rsky", rsky,METH_VARARGS,rsky_doc},{NULL}};

void initrsky(void)
{
  Py_InitModule("rsky", rsky_methods);
  import_array();
}

