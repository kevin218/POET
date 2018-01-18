# $Author: patricio $
# $Revision: 301 $
# $Date: 2010-07-10 03:33:44 -0400 (Sat, 10 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/psfbuild.py $
# $Id: psfbuild.py 301 2010-07-10 07:33:44Z patricio $

import numpy    as np
import time
import matplotlib.pyplot as plt
import imageedit  as ie
import interp2d as i2d
from scipy.ndimage.interpolation import zoom

def psfbuild(data, mask, fp, center, nimpos, resize=30, trim=10):

  """

  The resampling uses scipy.ndimage.interpolation.zoom, the size of
  the result will have size = resize*size, the value of the corners
  coincide in both arrays.

  Work plan:

  + trim the images to a rectangular area around the star
  + bilinear interpolate images and masks
  + clip mask interpolation so that any mask value < 1 becomes 0
  + get centers from frameparameters (= from aperture photometry)
  + shift each image to align with center of first image
  + shift mask likewise
  + sum images
  + sum masks
  + multiply images by masks
  + psf = sum(image) / sum(mask)   (pixel-wise sum)
  + psf = psf - np.median(psf)    (to subtract sky)

import sys
sys.path.append('/home/patricio/ast/esp01/convert/lib/python')
import psfbuild as pb

import event      as ev
from loadevent import loadEvent
event = loadEvent('hd209bs51_ctr', load=['data','uncd','mask','event'])

psf = pb.psfbuild(event.data, event.mask, event.fp, event.srcest, 
                  event.nimpos, resize=20, trim=10)

fig = plt.figure(11)
plt.clf()
plt.imshow(psf[:,:,0], interpolation='nearest', origin='ll')
#plt.xlim(200,400)
#plt.ylim(200,400)
plt.colorbar()

  """

  tini = time.time()
  print('\nPSF Building:')

  mnpos, ny, nx, npos = np.shape(data)

  # Extract only a rectangular area around the star:
  subdata = np.zeros((mnpos, 2*trim+1, 2*trim+1, npos))
  submask = np.zeros((mnpos, 2*trim+1, 2*trim+1, npos))

  # Get shape of the sub-images
  sny, snx = np.shape(subdata[0,:,:,0])

  # Coordinates of the centers
  yc, xc = center[0::2], center[1::2]

  # Trim the image:
  for pos in np.arange(npos):
    for i in np.arange(mnpos):
      subdata[i,:,:,pos], submask[i,:,:,pos] = \
                   ie.trimimage( data[i,:,:,pos], (yc[pos],xc[pos]), 
                                 (trim, trim),    mask=mask[i,:,:,pos] )


  # Shape of resampled images:
  coory    = zoom(np.arange(sny)*1.0, resize, order=1)
  coorx    = zoom(np.arange(snx)*1.0, resize, order=1)
  rny, rnx = np.size(coory), np.size(coorx)


  # Get shifts according to centering:
  shifty = fp.y - fp.y[:,0:1]
  shiftx = fp.x - fp.x[:,0:1] 

  # Amount by which the images should be rolled to align them
  expandy = (rny-1.0)/(sny-1.0)
  expandx = (rnx-1.0)/(snx-1.0)
  rolly   = -np.round(shifty*expandy)
  rollx   = -np.round(shiftx*expandx)

  # Stupid python gets confused with -0 values!
  rolly[np.where(rolly==0)] = 0
  rollx[np.where(rollx==0)] = 0

  verbose = False
  if verbose:
    # Rounded shifts
    roundsy = -rolly/expandy
    roundsx = -rollx/expandx
    
    plt.figure(1,(10,5))
    plt.clf()
    plt.plot(shifty [3,:], 'r')
    plt.plot(roundsy[3,:], 'b')


    '''
    im  = 0
    pos = 0
    resize = 25
    z = zoom(subdata[i,:,:,pos], resize, order=1)
    z = zoom(subdata[i,:,:,pos], 501.0/sny, order=1)

    rny, rnx = (np.array([sny, snx]) - 1) * resize + 1

    print(np.shape(z))
    print(rny)

    q = zoom(z, 1./resize, order=1)

    resize=3
    z = zoom(x, 51./3,    order=1)
    q = zoom(z, 1./resize, order=1)
    q

    rny, rnx = (np.array([3, 3]) - 1) * (resize) + 1
    np.linspace(0, 3-1, rny)
    np.size(np.linspace(0, 3-1, rny))
    i2d.interp2d(x, expand=resize+1, order=0)

    coor = zoom(np.arange(sny)*1.0, resize, order=1)

    '''


  verbose = False
  if verbose:
    fig = plt.figure(5)
    plt.clf()
    plt.imshow(subdata[ 0,:,:,3], interpolation='nearest', origin='ll')
    plt.colorbar()

    fig = plt.figure(6)
    plt.clf()
    plt.imshow(subdata[36,:,:,3], interpolation='nearest', origin='ll') 
    plt.colorbar()

  nfig = 7

  # Allocate space for the resampled data
  rsdata = np.zeros((mnpos, rny, rnx, npos))
  rsmask = np.zeros((mnpos, rny, rnx, npos))

#  for i in np.arange(mnpos):
#    for pos in np.arange(npos):
  # Resample and align:
  for pos in np.arange(npos):
    for i in np.arange(nimpos[pos]):

      # Resample using linear interpolation:
      rsdata[i,:,:,pos] = zoom(subdata[i,:,:,pos], resize, order=1)
      rsmask[i,:,:,pos] = zoom(submask[i,:,:,pos], resize, order=1)
      rsmask[i,:,:,pos] = (rsmask[i,:,:,pos] == 1)

#       rsdata[i,:,:,pos] = i2d.interp2d(subdata[i,:,:,pos], expand=resize, 
#                                        y=y, x=x, yi=yi, xi=xi)
#       rsmask[i,:,:,pos] = i2d.interp2d(submask[i,:,:,pos], expand=resize, 
#                                        y=y, x=x, yi=yi, xi=xi)

#       fig = plt.figure(nfig)
#       nfig += 1
#       plt.clf()
#       plt.imshow(rsdata[i,:,:,pos], interpolation='nearest', origin='ll')
#       plt.xlim(250,350)
#       plt.ylim(250,350)
#       plt.colorbar()

      # Align the data using roll:
      rsdata[i,:,:,pos] = np.roll(np.roll(rsdata[i,:,:,pos], 
                                          int(rolly[pos,i]),axis=0), 
                                          int(rollx[pos,i]),axis=1)

      rsmask[i,:,:,pos] = np.roll(np.roll(rsmask[i,:,:,pos], 
                                          int(rolly[pos,i]),axis=0), 
                                          int(rollx[pos,i]),axis=1)
#      if i%15 == 0:  # progress
#        print(i)


  print('resizing time: %f seconds'%(time.time()-tini))

  verbose = False
  if verbose:
    fig = plt.figure(2)
    plt.clf()
    plt.imshow(rsmask[0,:,:,3], interpolation='nearest', origin='ll')
    plt.colorbar()

    fig = plt.figure(3)
    plt.clf()
    plt.imshow(rsmask[36,:,:,3], interpolation='nearest', origin='ll')
#    plt.xlim(250,350)
#    plt.ylim(250,350)
    plt.colorbar()

  verbose = False
  if verbose:
    fig = plt.figure(9)
    plt.clf()
    plt.imshow(rsdata[36,:,:,3], interpolation='nearest', origin='ll')
    #plt.xlim(250,350)
    #plt.ylim(250,350)
    plt.colorbar()


  # Multiply images by masks
  rsdata *= rsmask

  # Number of good values in each pixel
  ngood = np.sum(rsmask, axis=0)
  # Avoid dividing by 0
  loc = np.where(ngood == 0)
  ngood[loc] = 1.0

  # The PSF
  psf = np.sum(rsdata, axis=0) / ngood

  # Subtract the sky level:
  psf -= np.median(psf)

  # Fix bad pixels
  psf[loc] = 0.0

  # Normalize to sum(psf) = 1.0
  # psf /= np.sum(psf)

  # Plots:
  # plot(shifty[3,:])
  # resize  = 30


  return psf
