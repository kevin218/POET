import numpy as np
import spiderman

def spiderman_zhang(rampparams, t, etc = []):
   """
  This function creates a model that fits a physical motivated model based on Zhang et al. 2017, ApJ, 836, 73

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
	xi:		radiative/advective timescale
	T_n:		nightside temperature
	delta_T:	day-night temperature contrast
        npoints:        number of phase bins for light curve interpolation
	

  Returns
  -------
	This function returns an array of y values...

  Revisions
  ---------
  2016-11-19 	Laura Kreidberg	
                laura.kreidberg@gmail.com 
                Original version
  2019-02-24	update interpolation, add to github version 
  TODO          add response function to etc
   """
   p = spiderman.ModelParams(brightness_model =  'zhang', stellar_model = 'blackbody')

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
   p.xi	   	    = rampparams[13]
   p.T_n	    = rampparams[14]
   p.delta_T	    = rampparams[15]
   npoints          = int(rampparams[16])

   #TODO: add filter path to etc
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


