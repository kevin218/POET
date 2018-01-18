# $Author: patricio $
# $Revision: 301 $
# $Date: 2010-07-10 03:33:44 -0400 (Sat, 10 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/psffit.py $
# $Id: psffit.py 301 2010-07-10 07:33:44Z patricio $

import numpy as np
from scipy.optimize              import leastsq
from scipy.ndimage.interpolation import zoom, shift

def scalepsf(psf, pars, shape, psfc, order, norm):

  """
    Shifts, reshapes and scales psf image.

    Parameters
    ----------
    psf     : 2D ndarray
              array containing a psf, image that will be fitted to data. 
              Doesn't need to have the same size of data.
    pars    : 4-elements tuple containing the parameters to fit.
              pars[0], pars[1] are the x,y coordinates of the center.
              pars[2] is the aperture flux.
              pars[3] is the sky level flux.
    shape   : The shape of the output psf, rescaled using 
              scipy.ndimage.interpolation.zoom
    psfc    : 2 elements tuple
              Contains the y,x pair of the center of the psf image.
    order   : Perform a bilinear interpolation of the psf (order=1), or a 
              nearest neighbor interpolation (order=0).
    norm    : If set to 1, normalizes the psf to sum(psf) = 1.0 before 
              scaling it. Default 0.


    Returns
    -------
    scaledpsf : 2D ndarray
                The psf shifted, rescaled and scaled according to pars.

    Examples
    --------
    import psffit as pf
    reload(pf)

    pars = [10.202, 9.78, 2012, 987]
    fit = pf.scalepsf(psf[:,:,0], pars, [21,21],
                      psfc, 1, 1)
    plt.figure(1)
    plt.clf()
    plt.imshow(fit, interpolation='nearest', origin='ll')

    # go psffit and try to fit it

    Revisions
    ---------
    2010-07-06  Patricio  Initial version, based     pcubillos@fulbrightmail.org
                          on CF_MIPS1_PSFCTR code 
                          written by jh 2004-03-16.
  """

  psfshape = np.shape(psf)

#   # Shift of the image:
#   shifty = psfc[0] - pars[0]
#   shiftx = psfc[1] - pars[1]

#   expandy = (psfshape[0]-1.0)/(shape[0]-1.0)
#   expandx = (psfshape[1]-1.0)/(shape[1]-1.0)
#   rolly   = int(-np.round(shifty*expandy))
#   rollx   = int(-np.round(shiftx*expandx))

#   resize = psfshape[0]/shape[0]

#   # Use roll to center the psf at [pars[0], pars[1]]
#   scaledpsf = np.roll(np.roll(psf, rolly, axis=0), rollx, axis=1)

  expandy = (psfshape[0]-1.0)/(shape[0]-1.0)
  expandx = (psfshape[1]-1.0)/(shape[1]-1.0)

  # Shift of the image:
  shifty = (pars[0] - psfc[0])*expandy
  shiftx = (pars[1] - psfc[1])*expandx

  # Use roll to center the psf at [pars[0], pars[1]]
  scaledpsf = shift(psf, [shifty, shiftx], order=order, mode='wrap')


  # Reshape the psf:
  resize = psfshape[0]/shape[0]
  scaledpsf = zoom(scaledpsf, 1.0/resize, order=order)

  if norm:
    scaledpsf /= np.sum(scaledpsf)

  # Scale the psf: using [pars[2], pars[3]
  scaledpsf = scaledpsf*pars[2] + pars[3]

  return scaledpsf


def residuals(pars, data, psf, weights, psfc, order, norm):
  """
    Returns a weighted and flattened to 1D array of the squared 
    differences between data and psf image.

    Parameters
    ----------
    data    : 2D ndarray
              array giving a data image.
    weights : 2D ndarray
              array giving the weights of data. same size of data.
    psf     : 2D ndarray
              array containing a psf, image that will be fitted to data. 
              Doesn't need to have the same size of data.
    pars    : parameters to fit.
    psfc    : Passed to scalepsf.
    order   : Passed to scalepsf. Perform a bilinear (1), or nearest 
              neighbor interpolation (0). 
    norm    : Passed to scalepsf. If set to 1 normalizes the psf.


    Returns
    -------
    params : Best fitting parameters.
    sigma  : Params' variances.
    niter  : Number of iterations.
    cspdof : Chi-square value.
    status : An integer flag. If it is equal to 1, 2, 3 or 4, the
             solution was found. Otherwise, the solution was not found.

    Notes
    -----

    Examples
    --------

    Revisions
    ---------
    2010-07-06  Patricio  Initial version.  pcubillos@fulbrightmail.org

  """
  shape = np.shape(data)
  scaledpsf = scalepsf(psf, pars, shape, psfc, order, norm)

  # Trim 2 pixels of the edges
  resid = ((data - scaledpsf)*np.sqrt(weights))[2:-2,2:-2]

  return resid.flatten()


def psffit(data, weights, psf, pars, psfc, order=1, norm=1):

  """
    Use least squares to fit the psf to the data.

    Parameters
    ----------
    data    : 2D ndarray
              array giving a data image.
    weights : 2D ndarray
              array giving the weights of data. same size of data.
    psf     : 2D ndarray
              array containing a psf, image that will be fitted to data. 
              Doesn't need to have the same size of data.
    pars    : parameters to fit.
    psfc    : Center coordinates of the psf.
    order   : Passed to scalepsf. Perform a bilinear interpolation of the 
              psf (order=1), or a nearest neighbor interpolation (order=0).
    norm    : Passed to scalepsf. If norm=1 normalizes the psf to sum(psf)=1

    Returns
    -------
    A tuple with 
    params : Best fitting parameters.
    sigma  : Params' variances.
    niter  : Number of iterations.
    cspdof : Chi-square value.
    status : An integer flag. If it is equal to 1, 2, 3 or 4, the
             solution was found. Otherwise, the solution was not found.

    Notes
    -----
    
    Examples
    --------

    weights = 1.0 / subun**2.0
    weights[np.where(weights == 0)] = np.mean(weights)/100.0 
    thefit = pf.psffit(fit, weights, psf[:,:,0], fitpars, psfc, 1, 1)
    pars    = [10.202, 9.58, 2000, 1000]
    pars    = [10.202, 9.78, 2012, 987]
    fitpars = [10.065, 8.48, 1500, 700]

tini = time.time()
plsq = leastsq(pf.residuals, fitpars, args=(fit, psf[:,:,0], weights, psfc, 
                                         1, 1), full_output=True)
print(time.time()-tini)

    reslt = pf.scalepsf(psf[:,:,0], plsq[0], [21,21], psfc, 1, 1)
    plt.figure(2)
    plt.clf()
    plt.imshow(reslt, interpolation='nearest', origin='ll') 

    Revisions
    ---------
    2010-07-06  Patricio  Initial version.  pcubillos@fulbrightmail.org

  """

  plsq = leastsq(residuals, pars, args=(data, psf, weights, psfc, order, 
                                        norm), full_output=True)  

  # The fitted parameters
  params = plsq[0]

  # Errors of the fitted parameters
  sigma  = np.diagonal(np.sqrt(plsq[1]))

  # Number of iterations
  niter  = plsq[2]['nfev']

  chisq  = np.sum( (plsq[2]['fvec'])**2.0 )
  dof    = len(plsq[2]["fvec"]) - len(plsq[0])
  # Chi-square per degree of freedom
  cspdof = chisq / dof

  # Status, 1 thru 4 means that the solution was found
  status = plsq[4]

  return params, sigma, niter, cspdof, status
