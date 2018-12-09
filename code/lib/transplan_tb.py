"""
 NAME:
    TRANSPLAN_TB

 PURPOSE:
    This function calculates the brightness temperature of a
    transiting planet given the observed ratios of areas and
    fluxes, the filter transmission function, stellar gravity and
    temperature, and a grid of Kurucz models in which to
    interpolate the intrinsic stellar brightness.

 CATEGORY:
    Astronomy

 CALLING SEQUENCE:

    Tb = TRANSPLAN_TB(Arat, Frat, Logg, Tstar)

 INPUTS:
    Arat:    Area ratio between planet and star (Rp^2/R*^2).  If
        this parameter is a 4-element array, they are assigned
        as [Arat, Frat, Logg, Tstar], and the other three
        positional parameters are ignored if present.  This is
        useful when calling this function from a Monte-Carlo
        error calculator.
    Frat:    Flux ratio between planet and star (Fp/F*)

    Logg:    Stellar log(gravity) (cm s-2).  NOTE: NOT MKS!

    Tstar:    Stellar temperature (K).

 OPTIONAL INPUTS:
    If GUESS is set and fmfreq and fmstar are given, Freq,
    Trans, and Fstar are all optional.  If Arat is a 4-element
    array, Frat, Ffreq, and Ftrans are ignored, and if present are
    assigned to the corresponding value of Arat IN THE CALLER.

 KEYWORD PARAMETERS:
    FFREQ:    [nfreq] array of frequencies (Hz) at which Trans and
        Fstar are tabulated.  MUST BE IN ASCENDING ORDER!
    FTRANS:    [nfreq] array giving the fractional filter
        transmission at each frequency in Freq.
    GMULT:    [3] array of multipliers of TBG that produces input to
        FX_ROOT_JH, default [1d, 0.6, 1.6].
    GUESS:    Don't solve the integral, just do the guess.
    DOUBLE, ITMAX, STOP, TOL, VERBOSE: Passed to FX_ROOT_JH().
    FMFREQ:    (input or returned) Representative frequency (Hz) of
        the filter.  If set, is used in calculating brightness
        temperature guess.  If not set, is calculated from the
        transmission function and returned.
    FMSTAR:    (input or returned) stellar brightness at
            FMFREQ (W m-2 sr-1 Hz-1).  If set, is used in
            calculating brightness temperature guess.  If
            not set, is calculated from the transmission
            function and returned.
    FSTAR:    (returned) [nffreq] interpolated stellar brightness at
            FFREQ (W m-2 sr-1 Hz-1).  If set, is used in
            calculating brightness temperature.  If
            not set, is calculated from the Kurucz
            parameters and returned.
    TBG:    (returned) Brightness temperature guess (K).
    STATUS:     (returned) Set to 0 if the algorithm did not
         converge, 1 if it did IN THE CALLER.

 OUTPUTS:
    This function returns the brightness temperature of a
    transiting planet in Kelvin.

 SIDE EFFECTS:
    Sets the values of the returned variables IN THE CALLER.

 PROCEDURE:
    See Eq. 6 and 10 of 2006-02-01-planetbrightnesstemp.pdf.  To
    get a double-precision calculation, use double-precision
    inputs.

 EXAMPLE:

 MODIFICATION HISTORY:
     Written by:    Joseph Harrington, Cornell.  2006-02-11
            jh@oobleck.astro.cornell.edu
    2006-02-12 jh    Made 4-element vector as first argument (opt,
            for Monte Carlo).  Added STATUS.
"""

def transplan_tb(arat, frat, logg, tstar, kfreq=None, kgrav=None, ktemp=None, kinten=None, ffreq=None, ftrans=None, fmfreq=None, fstar=None, fmstar=None, guess=None, gmult=None, itmax=None, stop=None, tol=None, verbose=None, status=None):

   import numpy as np
   import transplan_tb
   import kurucz_inten
   import integrate

   # unpack arat if array
   if np.array(arat).size == 4:
        iarat = arat[0]
        frat  = arat[1]
        logg  = arat[2]
        tstar = arat[3]
   else:
        iarat = arat
   
   status = 1 # guess always succeeds
   
   # default guess array
   if gmult == None:   
          gmult = np.array([1.0, 0.6, 1.6])
   
   # constants, MKS
   c = 2.99792458e8    # m/s, speed of light
   h = 6.6260693e-34   # J s, Planck's constant, Wikipedia
   k = 1.3806503e-23   # J/K, Boltzmann's constant, Google
   
   # find fstar
   if fstar is None:   
          kstar = kurucz_inten.interp(kinten, kgrav, ktemp, logg, tstar)
          fstar = np.interp(ffreq, kfreq, kstar)
   
   # calculate filter mean and brightness there
   if fmfreq is None:
    fmfreq = (integrate.integrate(ffreq, ffreq*ftrans, min(ffreq), max(ffreq)) /
          integrate.integrate(ffreq, ftrans,       min(ffreq), max(ffreq)))
    
   if fmstar is None:
          fmstar = np.interp([fmfreq], ffreq, fstar)
   
   # guess Tb at the weighted band center
   tbg = (h*fmfreq)/k/np.log((2.*h*iarat*fmfreq**3)/(c**2*frat*fmstar) + 1.)
   if (guess is not None):   
      return -1, tbg[0], fmfreq, fmstar
   
   # fx_root parameters
   guess3 = tbg * gmult
   func = tbfunc
   extra = [iarat, frat, ffreq, ftrans, fstar]
   #Convert line below
   tb = transplan_tb.root(guess3, func, extra=extra)
   
   return tb, tbg, fmfreq, fmstar   #, fmfreq, fstar, fmstar


#
# Copyright (c) 1994-2005, Research Systems, Inc.  All rights reserved.
#       Unauthorized reproduction prohibited.
# Modifications by Joseph Harrington are in the public domain.
#+
# NAME:
#       ROOT
#
# PURPOSE:
#       This function computes real and complex roots (zeros) of
#       a univariate nonlinear function.  This version improves on
#       that in the IDL release by offering _EXTRA, STATUS, and a
#       sanity check on TOL.
#
# CATEGORY:
#       Nonlinear Equations/Root Finding
#
# CALLING SEQUENCE:
#       Result = ROOT(X, Func)
#
# INPUTS:
#       X :      A 3-element initial guess vector of type real or complex.
#                Real initial guesses may result in real or complex roots.
#                Complex initial guesses will result in complex roots.
#
#       Func:    A scalar string specifying the name of a user-supplied IDL
#                function that defines the univariate nonlinear function.
#                This function must accept the vector argument X.
#
# KEYWORD PARAMETERS:
#       DOUBLE:  If set to a non-zero value, computations are done in
#                double precision arithmetic.
#
#       ITMAX:   Set this keyword to specify the maximum number of iterations
#                The default is 100.
#
#       STOP:    Set this keyword to specify the stopping criterion used to
#                judge the accuracy of a computed root, r(k).
#                STOP = 0 implements an absolute error criterion between two
#                successively-computed roots, |r(k) - r(k+1)|.
#                STOP = 1 implements a functional error criterion at the
#                current root, |Func(r(k))|. The default is 0.
#
#       TOL:     Set this keyword to specify the stopping error tolerance.
#                If the STOP keyword is set to 0, the algorithm stops when
#                |x(k) - x(k+1)| < TOL.
#                If the STOP keyword is set to 1, the algorithm stops when
#                |Func(x(k))| < TOL. The default is 1.0e-4.
#         Tol is limited to machine precision.  If set below
#         precision, it will be reset to precision IN THE
#         CALLER.
#
#    STATUS:     (returned) Set to 0 if the algorithm did not
#         converge, 1 if it did IN THE CALLER.
#
#    EXTRA:     Structure containing parameters to pass to FUNC.
#
# EXAMPLE:
#       Define an IDL function named FUNC.
#         function FUNC, x
#           return, exp(sin(x)^2 + cos(x)^2 - 1) - 1
#         end
#
#       Define a real 3-element initial guess vector.
#         x = [0.0, -!pi/2, !pi]
#
#       Compute a root of the function using double-precision arithmetic.
#         root = FX_ROOT(x, 'FUNC', /double)
#
#       Check the accuracy of the computed root.
#         print, exp(sin(root)^2 + cos(root)^2 - 1) - 1
#
#       Define a complex 3-element initial guess vector.
#         x = [complex(-!pi/3, 0), complex(0, !pi), complex(0, -!pi/6)]
#
#       Compute a root of the function.
#         root = FX_ROOT(x, 'FUNC')
#
#       Check the accuracy of the computed root.
#         print, exp(sin(root)^2 + cos(root)^2 - 1) - 1
#
# PROCEDURE:
#       FX_ROOT implements an optimal Muller's method using complex
#       arithmetic only when necessary.
#
# SIDE EFFECTS:
#    Sets STATUS and may set TOL IN THE CALLER.
#
# REFERENCE:
#       Numerical Recipes, The Art of Scientific Computing (Second Edition)
#       Cambridge University Press
#       ISBN 0-521-43108-5
#
# MODIFICATION HISTORY:
#       Written by:  GGS, RSI, March 1994
#       Modified:    GGS, RSI, September 1994
#                    Added support for double-precision complex inputs.
#    2005-02-07 jh    Added _extra.
#    2005-02-12 jh    Added status, tol protection.  Fixed indentation.
#    2008-08-15 kevin
#            Kevin Stevenson, UCF
#            Converted to Python
#-

def root(xi, func, itmax=None, stop=None, tol=None, status=None, extra=None):

   import numpy as np
   import transplan_tb

   status = 0
   
   x  = xi + 0.0   #Create an internal floating-point variable, x.
   if x.size!=3:
      print('x must be a 3-element initial guess vector.')
      return np.nan
   
   if itmax is None:   
      itmax = 100
   if stop is None:   
      stop = 0
   if tol is None:   
      tol = 1.0e-4

   #Initialize stopping criterion and iteration count.
   cond = 0
   it   = 0
   
   #Begin to iteratively compute a root of the nonlinear function.
   while (it < itmax and cond != 1):
      q    = (x[2] - x[1]) / (x[1] - x[0])
      pls  = (1 + q)
      f    = func(x, extra)
      a    = q * f[2] - q * pls * f[1] + q ** 2 * f[0]
      b    = (2 * q + 1) * f[2] - pls ** 2 * f[1] + q ** 2 * f[0]
      c    = pls * f[2]
      disc = b ** 2 - 4 * a * c
      roc  = disc.size  #Real or complex discriminant?
      if (roc != 6 and roc != 9):    #Proceed toward real root.
         if disc < 0:    #Switch to complex root.
            r0 = b + complex(0, np.sqrt(abs(disc)))
            r1 = b - complex(0, np.sqrt(abs(disc)))
            if abs(r0) > abs(r1):   
               div = r0
            else:   
               div = r1
         else:           # real root
            rr0 = b + np.sqrt(disc)
            rr1 = b - np.sqrt(disc)
            div = ((abs(rr0) >= abs(rr1)) and [rr0] or [rr1])[0]
      else:               #Proceed toward complex root.
         c0 = b + np.sqrt(disc)
         c1 = b - np.sqrt(disc)
         if abs(c0) > abs(c1):   
            div = c0
         else:   
            div = c1
      root = x[2] - (x[2] - x[1]) * (2 * c / div)
      #Absolute error tolerance.
      if (stop == 0 and abs(root - x[2]) <= tol):   
         cond = 1
      else:   
         evalfunc = func(root, extra)
         #Functional error tolerance.
         if (stop != 0 and abs(evalfunc) <= tol):   
            cond = 1
         else:   
            if evalfunc == 0:   
               cond = 1
      x = np.array((x[1], x[2], root))
      it = it + 1

   if (it >= itmax and cond == 0):   
      print('Algorithm failed to converge within given parameters.')
   else:   
      status = 1
   
   return root

"""
 NAME:
    TBFUNC

 PURPOSE:
    This function is used by fx_root_jh to calculate the
    brightness temperature of a transiting planet.  When
    this function returns 0, a root has been found and Tb is the
    brightness temperature of the planet in Kelvin.

 CATEGORY:
    Astronomy.

 CALLING SEQUENCE:

    Result = tbfunc(Tb)

 INPUTS:
    Tb:    Brightness temperature (Kelvin)

 KEYWORDS:
    ARAT:    Area ratio between planet and star (Rp^2/R*^2)
    FRAT:    Flux ratio between planet and star (Fp/F*)
    FREQ:    [nfreq] array of frequencies (Hz) at which Trans and
        Istar are tabulated.  MUST BE IN ASCENDING ORDER!
    TRANS:    [nfreq] array giving the fractional filter
        transmission at each frequency in Freq.
    ISTAR:    [nfreq] array giving the stellar brightness at each
        frequency in Freq in units of W m-2 sr-1 Hz-1.
    VERBOSE:    Print input tb and function values.

 OUTPUTS:
    This function returns the vertical distance from the Y axis
    for the given input parameters.  If the distance is 0, a root
    has been found and Tb is the brightness temperature of the
    planet in Kelvin.

 PROCEDURE:
    See Eq. 6 of 2006-02-01-planetbrightnesstemp.pdf.  To get a
    double-precision calculation, use double-precision inputs.

 EXAMPLE:
    ; Note: fx_root_jh just adds _extra=e to the function
    ; definition and the internal named function call, to allow
    ; passing of extra parameters.
    func = 'transplan_tbfunc'
    tb = fx_root_jh(guess, func, $
         arat=arat, frat=frat, freq=freq, trans=trans, istar=istar, $
         itmax=itmax, stop=stop, tol=tol)

 MODIFICATION HISTORY:
     Written by:    Joseph Harrington, Cornell.  2006 Feb 7
            jh@oobleck.astro.cornell.edu
    2006 02 10 jh    Fix c, add VERBOSE.
    2008-08-15 kevin
            Kevin Stevenson, UCF
            Converted to Python
"""


def tbfunc(tb, extra, verbose=None):

   import numpy as np
   import integrate

   arat  = extra[0]
   frat  = extra[1]
   freq  = extra[2]
   trans = extra[3]
   istar = extra[4]

   c = 2.99792458e8    # m/s, speed of light
   h = 6.6260693e-34   # J s, Planck's constant, Wikipedia
   k = 1.3806503e-23   # J/K, Boltzmann's constant, Google
   
   if (verbose is not None):   
          print(tb)
   
   ntb = np.array(tb).size
   if ntb == 1:
    ret    = c**2 * frat / (2.*h*arat) * integrate.integrate(freq, trans * istar) - \
           integrate.integrate(freq, trans*freq**3 / (np.exp(h*freq/(k*tb)) - 1.))
   else:
     ret = np.zeros(ntb)
     for i in range(ntb):
          ret[i] = c**2 * frat / (2.*h*arat) * integrate.integrate(freq, trans * istar) - \
           integrate.integrate(freq, trans*freq**3 / (np.exp(h*freq/(k*tb[i])) - 1.))

   if (verbose is not None):   
          print(ret)
   
   return ret


