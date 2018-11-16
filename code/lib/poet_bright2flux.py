# $Author: carthik $
# $Revision: 267 $
# $Date: 2010-06-08 22:33:22 -0400 (Tue, 08 Jun 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_bright2flux.py $
# $Id: poet_bright2flux.py 267 2010-06-09 02:33:22Z carthik $

import numpy as np
from scipy.constants import arcsec

def poet_bright2flux(data, uncd, posscl):
  """
    This function converts the data and uncertainty arrays from
    brightness units (MJy/sr) to flux units (uJy/pix).

    Parameters:
    -----------
    data:    ndarray 
             data array of shape ([nx, ny, nimpos, npos]) in units of MJy/sr.
    uncd:    ndarray 
             uncertainties of data (same shape and units).
    posscl:  ndarray
             Pixel scale (arcsec/pix) in each position
	     of shape [2, npos].

    Return:
    -------
    This procedure returns the input arrays Data and Uncd into
    flux units (uJy/pix), if they are defined in the input.

    Notes:
    ------
    The input arrays Data and Uncd are changed in place.

    Modification History:
    ---------------------
    2005-06-20 statia    Written by  Statia Luszcz, Cornell.
                         shl35@cornell.edu
    2005-10-13 jh        Renamed, modified doc, removed posmed, fixed
		         nimpos default bug (was float rather than int).
    2005-10-28 jh        Updated header to give units being converted
		         from/to, made srperas value a calculation
		         rather than a constant, added Allen reference.
    2005-11-24 jh        Eliminated NIMPOS.
    2008-06-28 jh        Allow npos=1 case.
    2010-01-29 patricio  Converted to python. pcubillos@fulbrightmail.org
    2010-11-01 patricio  Documented, and incorporated scipy.constants.
  """
  # steradians per square arcsecond
  srperas = arcsec**2.0

  npos = data.shape[3]
  for pos in np.arange(npos):
    data[:, :, :, pos] *= srperas * 1e12 * posscl[0, pos] * posscl[1, pos] 
    uncd[:, :, :, pos] *= srperas * 1e12 * posscl[0, pos] * posscl[1, pos] 

  return data, uncd

