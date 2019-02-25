import numpy as np
import spiderman


def spiderman_sph(rampparams, t, etc = []):
   """
  This function creates a model that fits spherical harmonics 

  Parameters
  ----------
    	t0:		time of conjunction 
	per:		orbital period
	a_abs:		semi-major axis (AU)
	cos(i):	        cosine of the orbital inclination	
	ecc:		eccentricity
	w:		arg of periastron (deg)
	rp:		planet radius (stellar radii)
	a:		semi-major axis (stellar radii)
	p_u1:		planet linear limb darkening parameter
	p_u2:		planet quadratic limb darkening
	T_s:		stellar Teff
	l1:		short wavelength (m)
	l2:		long wavelength (m)
	degree:		maximum degree of harmonic (typically no more than 2)
	la0:		latitudinal offset of coordinate center from substellar point (degrees)
	lo0:		latitudinal offset of coordinate center from substellar point (degrees)
	sph0:		coefficient1
	sph1:		coefficient2	
	sph2:		coefficient3
	sph3:		coefficient4
        npoints:        number of phase bins for light curve interpolation


	**Note: 4 harmonic coefficients are needed for a model of degree 2

  Returns
  -------
	This function returns an array of planet/star flux values 

  Revisions
  ---------
  2016-11-19 	Laura Kreidberg	
                laura.kreidberg@gmail.com 
                Original version
  2019-02-24	update interpolation, add to github version 
  TODO          add response function to etc
   """
   p = spiderman.ModelParams(brightness_model =  'spherical', thermal = True, stellar_model = 'blackbody')

   p.t0    	    = rampparams[0]
   p.per       	    = rampparams[1]
   p.a_abs 	    = rampparams[2]
   p.inc	    = np.arccos(rampparams[3])*180./np.pi
   p.ecc	    = rampparams[4]
   p.w	   	    = rampparams[5]
   p.rp	    	    = rampparams[6]
   p.a	   	    = rampparams[7]
   p.p_u1	    = rampparams[8]
   p.p_u2	    = rampparams[9]
   p.T_s	    = rampparams[10]
   p.l1	   	    = rampparams[11]
   p.l2	    	    = rampparams[12]
   p.degree 	    = int(rampparams[13])
   p.la0	    = rampparams[14]
   p.lo0	    = rampparams[15]
   sph0	    = rampparams[16]
   sph1	    = rampparams[17]
   sph2	    = rampparams[18]
   sph3	    = rampparams[19]
   p.n_layers = 5

   sph = []
   for i in range(0, p.degree**2): sph.append(rampparams[16+i])
   p.sph = sph

   npoints = int(rampparams[20])

 
   #p.filter = "/Users/lkreidberg/Desktop/Util/Throughput/spitzer_irac_ch2.txt"


   #calculate light curve over npoints phase bins
   phase = (t - p.t0)/p.per
   phase -= np.round(phase)

   phase_bin = np.linspace(phase.min(), phase.max(), npoints)
   t_bin = phase_bin*p.per + p.t0

   lc_bin = spiderman.web.lightcurve(t_bin, p)

   #interpolate the binned light curve to the original t array
   lc = np.interp(phase, phase_bin, lc_bin)

   return lc 


