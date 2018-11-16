# $Author: patricio $
# $Revision: 291 $
# $Date: 2010-07-01 13:31:46 -0400 (Thu, 01 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_targpos.py $
# $Id: poet_targpos.py 291 2010-07-01 17:31:46Z patricio $

import numpy     as np
import asym      as asy
import centroid  as cen
import fitgaussprofile as fgp
from   scipy.ndimage import center_of_mass

def poet_targpos(posim, srcest, rad=0, method='fgc', larad=4):

  """
    This function estimates the target center at each position and
    returns a 2xN array giving the center in each position.

    Parameters
    ----------
    posim : 3D ndarray
            [nx, ny, npos] array giving the nominal images in each 
            set (position).
    srcest: 2 elements tuple
            An initial guess of the position of the star.  Has the form 
            (y, x) of the guess center. If None, this routine defines it.
    rad   : Scalar (positive)
            If rad!=0, trim the image in a box of 2*rad pixels around 
            the guess center. Must be !=0 for 'col' method.
    method: String
            The centering method. Can be:
            'fgc': Fit gauss function,        'col': Center of light fitting,
            'lag': Least asymmetry gaussian,  'lac': Least asymmetry col.
    larad : Scalar (positive)
            Least asymmetry radius.

    Returns
    -------
    This function returns a 2xN array giving the center in each position.

    Example:
    -------
    >>>

    Revisions
    ---------
    2005-10-12  Written by: Joseph Harrington      jh@oobleck.astro.cornell.edu
                            and Statia Luszcz      shl35@cornell.edu
                        Based on earlier code by both
    2005-10-13 jh	Added doc header, converted to function, cleaned up.
    2005-10-26 jh       Added Lim keyword, renamed pxmed to posim.
    2005-10-27 jh       Tweaked header.
    2007-06-28 jh       Handle npos=1 data.
    2010-07-12 patricio Implemented in python.     pcubillos@fulbrightmail.org
  """

  # get sizes
  ny, nx, npos = np.shape(posim)

  ctr = np.zeros((2, npos))

#   # If a guess is not supplied, find the greates value in posim.
#   if srcest == None:
#     srcest = np.zeros((2*npos))
#     for pos in np.arange(npos):
#       srcest[2*pos:2*pos+2] = (cen.ctrguess(posim[:,:,pos]))[1]

  for pos in np.arange(npos):
    # trim if requested
    if rad != None:
      ceny, cenx  =  srcest[:, pos] # srcest[2*pos:2*pos+2],
      image = np.copy(posim[ceny-rad:ceny+rad, cenx-rad:cenx+rad, pos])
      loc   = (rad, rad)
    else:
      image = posim[:,:,pos]
      loc   = srcest[:, pos]

    # Find the center
    if   method == 'fgc':
      ctr[:,pos] = cen.ctrgauss(image, loc)
    elif method == 'col':
      ctr[:,pos]  = center_of_mass(image) 
      #ctr[:,pos] += 0.5
    elif method in ['lag', 'lac']:
      ctr[:,pos] = asy.actr(image, loc, resize=True, method=method, rad=larad)
    elif method == 'fgp':
      ctr[:,pos] = fgp.center(image, loc)
    elif method == 'fap':
      pass

    # Do trim corrections
    ctr[:,pos] += srcest[:, pos] - loc
  return ctr
