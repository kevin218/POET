# $Author: patricio $
# $Revision: 313 $
# $Date: 2010-07-15 16:34:13 -0400 (Thu, 15 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/fitgaussprofile.py $
# $Id: fitgaussprofile.py 313 2010-07-15 20:34:13Z patricio $

#------------------------------------------------------------------
# PHASE 1

from numpy.random import normal
import numpy as np
import matplotlib.pyplot as plt
import time
import numpy as np

def gaussimg(pars, shape, sample):

  """
    Produce a gaussian image.

    Parameters
    ----------
    pars    : 5-elements tuple containing the parameters to fit.
              pars[0] : y coordinate of the center.
              pars[1] : x coordinate of the center.
              pars[2] : radial stretching factor.
              pars[3] : flux of the gaussian profile.
              pars[4] : sky level flux.
    shape   : The shape of the output image.
    n       : Scalar
              Number of points to calculate the cumulative probability 
              distribution.
    rad     : Scalar
              The radius whithin which the distribution of points are
              going to lay.
    nd      : Scalar
              Number of poits to sample the obscured airy pattern image.

    Returns
    -------
    gaussimage: 2D ndarray
                The psf shifted, rescaled and scaled according to pars.

    Examples
    --------

sys.path.append('/home/patricio/ast/esp01/convert/lib/python')
from imageedit import *
import fitgaussp as fp
reload(fp)

# go to an event:
event = loadevent('wa018bs41_pht')
os.chdir('..')
updateevent(event, 'wa018bs41_ctr', ['data', 'uncd', 'mask'])
os.chdir(event.photdir)
pos = 0
id = 350
image = event.data[id,:,:,pos]
mask  = event.mask[id,:,:,pos]
uncd  = event.uncd[id,:,:,pos]
print(event.fp.y[pos,id], event.fp.x[pos,id], event.fp.aplev[pos,id], 
      event.fp.skylev[pos,id])
center = np.round((event.fp.y[pos,id], event.fp.x[pos,id]))
subim, mask, uncd = trimimage(image, center, (10,10), mask=mask, uncd=uncd)

plt.figure(1, (6.5,5.5))
plt.clf()
plt.imshow(subim, interpolation='nearest', origin='ll', cmap=plt.cm.gray)
plt.colorbar()
plt.title('wa018bs41 data image')

appars = [22.3104223, 24.189664, 0.0, 37246.19938, 79.6434958]
sample = fp.gsample(1e4)
pars = [10.3, 10.2, 0.8, 35000, 60]
model = fp.gaussimg(pars, np.shape(subim), sample)

plt.figure(2)
plt.clf()
plt.imshow(model, interpolation='nearest', origin='ll', cmap=plt.cm.gray)
plt.colorbar()


    Revisions
    ---------
    2010-07-14  Patricio  initial version          pcubillos@fulbrightmail.org
  """

  y = pars[0] + sample[0]*pars[2] + 0.5  # [pars[0], pars[2]]
  x = pars[1] + sample[1]*pars[2] + 0.5  # [pars[1]]

  gaussprofile, xe, ye = np.histogram2d(y, x, bins=shape, 
                                        range=[[0,shape[0]],[0,shape[1]]])

  gaussprofile  /= np.sum(gaussprofile)

  # Scale the image:
  gaussprofile = gaussprofile*pars[3] + pars[4]  # [pars[3], pars[4]]

  return gaussprofile



def gausschisq(pars, data, weights, sample, mask=None):
  """
    Returns a weighted and flattened (1D) array of the squared 
    differences between data and psf image.

    Parameters
    ----------
    pars    : 6 scalar elements tuple 
              Parameters to fit.
    data    : 2D ndarray
              array giving a data image.
    weights : 2D ndarray
              array giving the weights of data. Same size of data.


    Returns
    -------
    This function return the weighted residuals of the difference
    between data and a obscured airy pattern image.


    Examples
    --------

    Revisions
    ---------
    2010-07-14  Patricio  Initial version.  pcubillos@fulbrightmail.org
  """

  shape = np.shape(data)
  if mask == None:
    mask = np.ones(shape, int)
  # get a airy image with parameters pars
  gaussp = gaussimg(pars, shape, sample)

  # Return the weighted chi squared
  residuals = ((data - gaussp)*np.sqrt(weights)*mask).flatten()
  dof = np.sum(mask) - len(pars)

  return np.sum(residuals**2.0) / dof


def gsample(n):
  y = normal(loc=0, scale=1, size=n)
  x = normal(loc=0, scale=1, size=n)
  return y, x


# Simplex method
def simplex(pars, data, weights, sample, ftol, mask=None,
            dels=[1.0, 1.0, 0.1, 100.0, 5.0]):

  """
    Use simplex method to minimize chi square

    Parameters
    ----------
    pars    : parameters to fit.

    Returns
    -------

    Examples
    --------

weights = subim * 0.0 + 1.0
weights = 1.0 / uncd**2.0
weights[np.where(weights == 0)] = np.mean(weights)/100.0
mask[16,10] = 0
shape = np.shape(subim)
sample = fp.gsample(1e4)
pars = [10.3, 10.2, 0.8, 35000.0, 60.0]
dels = [ 1.0,  1.0, 0.1,  1000.0, 10.0]

tini = time.time()
sol, niter, nfeval, chisq = fp.simplex(pars, subim, weights, sample, 1e-3, 
                                       dels=dels, mask=mask)
print(time.time()-tini)
sol, niter, nfeval, chisq

tini = time.time()
sol, niter, nfeval, chisq = fp.simplex(sol, subim, weights, sample, 1e-5, 
                                       dels=dels, mask=mask)
print(time.time()-tini)
sol, niter, nfeval, chisq

fit = fp.gaussimg(sol, np.shape(subim), sample)
plt.figure(6, (6.5,5.5))
plt.clf()
plt.imshow(fit, interpolation='nearest', origin='ll', cmap=plt.cm.gray)
plt.colorbar()
plt.title('Gauss profile fitting No weights')

plt.figure(4)
plt.clf()
plt.imshow(fit-subim, interpolation='nearest', origin='ll', cmap=plt.cm.gray)
plt.colorbar()

plt.figure(4)
plt.clf()
#plt.imshow(mask, interpolation='nearest', origin='ll', cmap=plt.cm.gray)
plt.imshow(uncd, interpolation='nearest', origin='ll', cmap=plt.cm.gray)
plt.colorbar()



    Revisions
    ---------
    2010-07-14  Patricio  Initial version.  pcubillos@fulbrightmail.org

  """

  itmax  = 5e3
  tiny   = 1e-10
  # indices of lowest, highest and next highest values 
  ilo, ihi, inhi = 0, 0, 0
  # pars dims, points
  npoints = len(pars) + 1
  npars   = len(pars)

  # initialize the amoeba
  amoeba = np.zeros((npoints, npars))
  y = np.zeros(npoints)
  for   i in np.arange(npoints):
    amoeba[i,:] = pars
    if i != 0:
      amoeba[i, i-1] += dels[i-1]
    y[i] = gausschisq(amoeba[i,:], data, weights, sample, mask=mask)

  psum  = np.sum(amoeba, axis=0)
  iter  = 0
  nfunc = 0

  while True:
    # get the highest, next hi and lowest index
    ihi  = np.argmax(y)
    ymask = np.ones(npoints)
    ymask[ihi] = 0
    inhi = np.argmax(y*ymask)
    ilo  = np.argmin(y)

    # Evaluate return condition
    rtol = 2.0*np.abs(y[ihi]-y[ilo])/(np.abs(y[ihi]) + np.abs(y[ilo]) + tiny)
    if rtol < ftol:
      return amoeba[ilo,:], iter, nfunc, y[ilo]

    if iter >= itmax:
      print(" ITMAX exceeded")
      return amoeba[ilo,:], iter, nfunc, y[ilo]

    # Begin new iteration
    # Reflect across the highest point
    ytry  = amotry(amoeba, y, psum, ihi, -1.0, data, weights, sample, mask=mask)
    nfunc += 1

    if ytry <= y[ilo]:
      # new point is better than best, try additional extrapolation
      ytry = amotry(amoeba, y, psum, ihi, 2.0, data, weights, sample, mask=mask)
      nfunc += 1

    elif ytry >= y[inhi]:
      # new point is worse than next hi, contract looking intermed. point
      ysave = y[ihi]
      ytry = amotry(amoeba, y, psum, ihi, 0.5, data, weights, sample, mask=mask)
      if ytry >= ysave:
        for i in np.arange(npoints):
          if i != ilo:
            amoeba[i,:] = psum = 0.5 * (amoeba[i,:] + amoeba[ilo,:])
            y[i] = gausschisq(psum, data, weights, sample, mask=mask)
        nfunc += npars + 1
        psum = np.sum(amoeba, axis=0)
    iter += 1


def amotry(amoeba, y, psum, ihi, fac, data, weights, sample, mask=None):
  ndim = len(psum)
  fac1 = (1.0-fac) / ndim
  fac2 = fac1 - fac
  ptry = psum*fac1 - amoeba[ihi,:]*fac2

  # Evaluate the function at the trial point
  ytry = gausschisq(ptry, data, weights, sample, mask=mask)
  # if trial is better than hi, replace it
  if ytry < y[ihi]:
    y[ihi] =  ytry
    psum += ptry - amoeba[ihi,:]
    amoeba[ihi,:] = ptry
  return ytry


def center(data, loc, mask=None, uncd=None):
  """
  returns the best fitting y, x parameters.
  """
  shape = np.shape(data)

  if mask == None:
    mask = np.ones(shape)
  if uncd == None:
    weights = np.ones(shape)
  else:
    weights = 1.0 / uncd**2.0
    weights[np.where(weights == 0)] = np.mean(weights)/100.0

  sample = gsample(1e4)

  # initial guess
  sky  = np.median(data)
  flux = np.sum((data-sky)*mask)
  pars = [loc[0], loc[1], 1.0, flux, sky]

  sol, niter, nfeval, chisq = simplex(pars, data, weights, sample, 1e-5, 
                                      mask=mask)

  return sol[0], sol[1]
