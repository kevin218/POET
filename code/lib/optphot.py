# $Author: patricio $
# $Revision: 301 $
# $Date: 2010-07-10 03:33:44 -0400 (Sat, 10 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/optphot.py $
# $Id: optphot.py 301 2010-07-10 07:33:44Z patricio $

import numpy as np

def optphot(image, psf, var=None, mask=None, nochecks=False):

  """
    This function performs optimally-weighted photometry on an image,
    analogous to the optimally-weighted spectroscopy of Horne (1986).

    Parameters:
    ----------
    image   : 2D floating image array containing object to measure.
              The image should have all corrections including sky
              subtraction.  Bad pixel correction is unnecessary.
    psf     : 2D floating image, same size, giving point-spread function.

    var     : 2D floating image, same size, giving variance.
    mask    : 2D byte array giving status of corresponding
              pixel in Image: bad pixel=0, good pixel=1.  All pixels
              are good if the mask is not specified.
    nochecks: Set to True to skip checks of input sanity.

    Outputs:
    -------
    immean  : This function returns the optimal estimate of the number of
              counts in the image, according to the PSF and variance.
    uncert  : Uncertainty in the derived measurement.
    profile : Normalized, positive-definite, non-zero version of
              PSF used in function.
    weights : Weights used in calculating measurement.

    Example:
    -------

    Modification History:
    --------------------
        Written by: Joseph Harrington <jh@oobleck.astro.cornell.edu>
        2004-03-16  jh        first version
        2005-01-21  jh        minor header tweak
        2010-06-29  patricio  converted to python   pcubillos@fulbrightmail.org
  """

  # Check inputs
  if not nochecks:
    if np.ndim(image) != 2:
      print('image must be a 2D array')
      return None

    if mask == None:
      mask = np.ones(np.shape(image), byte)
  
    if np.shape(image) != np.shape(mask): 
      print('image and mask sizes differ')
      return None

    if np.shape(image) != np.shape(psf):
      print('image and PSF sizes differ')
      return None

  # Calculate profile
  small = 1e-20 / np.size(image)
  profile = psf / np.sum(psf)

  profile = profile * mask
  profile[np.where(profile < small)] =  small

  smallim = np.ones(np.shape(profile)) * small
  profile = np.amax(np.array([[profile * mask],[smallim]]), axis=0)
  #profile = (profile * mask) > small
  
  # Estimate the variance (should really be done by caller to get sky
  # and number of electrons per DN)
  if var == None:
    var = image
    var[np.where(var < 0)] = 0.0

  # In each pixel, estimate the total flux
  ftotest   = image / profile
  vartotest = var   / profile**2.0

  # Do weighting (Bevington, p. 57)
  weights  = mask / vartotest
  immean   = np.sum(weights * ftotest) / np.sum(weights)
  imvar    = 1.0 / np.sum(weights)
  uncert   = np.sum(imvar)
  
  return immean, uncert, profile, weights
