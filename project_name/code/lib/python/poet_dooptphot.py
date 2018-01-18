# $Author: patricio $
# $Revision: 301 $
# $Date: 2010-07-10 03:33:44 -0400 (Sat, 10 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_dooptphot.py $
# $Id: poet_dooptphot.py 301 2010-07-10 07:33:44Z patricio $

import numpy     as np
import psffit    as pf
import psfbuild  as pb
import optphot   as op
import imageedit as ie
import logedit   as le
import time
import matplotlib.pyplot as plt

def dooptphot(data, uncd, mask, fp, center, nimpos,
              rejlim=[1.45, 1, 0.005], itmax=20, trim=10, order=1,
              ftol=1e-8, resize=30, norm=1, log=None, psf=None):

  """
   Perform optimal photometry on a set of images.  

   Parameters
   ----------
   data:   4D array_like
           image array (im, x, y, pos)

   uncd:   4D array_like
           Array containing pixel uncertainties (same shape as data)

   mask: Array of same dimensions as data containing the bad
           pixel masks to be used in optimal photometry.  

   sratio: array_like
           Pixel scale ratios for each position in x and y. The
           array should have shape: (2, npos)

   psf:    Passed directly into CURVEFIT the psfdata image (handed
           through to CF_MIPS1_PSFCTR unmodified)

   rejlim: A 3-element array containing the cutoff values for chi
           square, sky-median(sky)(median taken per position), and
           sky error(above minimum sky error for that
           position). Images that produce values outside the
           cutoff range will be marked as such in IOSTATUS. The
           default values for these parameters are: 
           rejlim = [1.45, 1, 0.005]

   ftol:   Relative error desired in the sum of squares. Passed to
           scipy.optimize.leastsq.

   res:    passed into CURVEFIT

   trim:   Define half-width of sub-image

  
   Returns:
   -------
     Array containing a whole bunch of information. 

   OPTIONAL OUTPUTS:

   OSTATUS: An array that labels bad images (according to the
            criterion above- see REJLIM).An image that fails based on
            a high chi-square is labeled with a -1. A high or low sky
            level is labeled with a -10, and a high sky error is
            labeled with a -100.
  
   SIDE EFFECTS:
        This function defines pos, nx, maxnim, and nimpos if they arent
        already defined. 
  
    Examples
    --------
ofp = do.dooptphot(data, uncd, mask,  fp, center, nimpos,
              rejlim=[10.45, 10, 1.5], itmax=20, order=1,
              resize=30, norm=1)

pos = 0
plt.clf()
plt.plot(ofp.time[pos, 0:nimpos[pos]], fp.oskylev [pos, 0:nimpos[pos]], '+')
plt.plot(ofp.time[pos, 0:nimpos[pos]], fp.ophotlev[pos, 0:nimpos[pos]], '+')

   Revisions:
   ---------

   Written by Statia June 20, 2005

   2005-06-21 Statia   Added image rejection part.
   2005-07-14          Added sclpar keyword.
   2005-07-15          Distprf keyword.
   2008-03-07          Changed parinfo to fix xctr, yctr (to use values
                         from apphot).
   2008-06-03 Emily    Rewriting to use curvefit instead of mpfit,
                         as was done in mips1-3b-optphot.pro
   2008-06-05 Emily    (sloppy) MIPS correction for sky median rejection test
   2009-05-18 Kevin    Allowed curvefit to use subimages of data.
   2010-06-24 Patricio Converted to python.
  """

  tini = time.time()

  # pieces of apphot needed (from fp): x, y, aplev, skylev

  # Sizes
  maxnimpos, ny, nx, npos = np.shape(data)

  # initialize skyint and skyslope (used in calculating status and
  # returned if keywords set)
  skyint   = np.zeros(npos)
  skyslope = np.zeros(npos)

  # Allocate space for the results
  fp.oxx        = np.zeros((npos, maxnimpos)) # X position
  fp.oyy        = np.zeros((npos, maxnimpos)) # Y position
  fp.oxerr      = np.zeros((npos, maxnimpos)) # error in X position
  fp.oyerr      = np.zeros((npos, maxnimpos)) # error in Y position
  fp.optlev     = np.zeros((npos, maxnimpos)) # flux from curvefit
  fp.opterr     = np.zeros((npos, maxnimpos)) # flux error curvefit
  fp.oskylev    = np.zeros((npos, maxnimpos)) # sky level
  fp.oskyerr    = np.zeros((npos, maxnimpos)) # sky error
  fp.ophotlev   = np.zeros((npos, maxnimpos)) # flux from optphot
  fp.ophoterr   = np.zeros((npos, maxnimpos)) # flux error from optphot
  fp.oniter     = np.zeros((npos, maxnimpos)) # iteration number for curvefit
  fp.ofitstatus = np.zeros((npos, maxnimpos)) # status output from mpfiot
  fp.ocspdof    = np.zeros((npos, maxnimpos)) # chi squared per dof
                                              # according to curvefit
  fp.ostatus    = np.zeros((npos, maxnimpos)) # optphot return status
                                              # (formally reject)
  fp.ogood      = np.zeros((npos, maxnimpos)) # good flag
  fp.orad       = np.zeros((npos, maxnimpos)) # distance from nearest pix center



  # Get the PSF:
  if psf == None:
    psf = pb.psfbuild(data, mask, fp, center, nimpos, resize=resize, trim=trim)


  # Do photometry

  # The parameters to fit:
  # pars = [Y position, X position, aperture flux, sky flux level]
  pars = np.zeros(4)


  # Pixels of the centers
  yr, xr = center[0::2], center[1::2]  # integer values
  yoff = yr - trim
  xoff = xr - trim

  # Coordinates of the center of the psf in the subimage:
  psfyc, psfxc = fp.y[:,0]-yoff, fp.x[:,0]-xoff

  le.writelog('\nOptimal Photometry\n', log)
  # pos = 0
  # if True:
  for  pos in np.arange(0, npos):
    for im in np.arange(0, nimpos[pos]):
      # calculate shifted, scaled PSF
      # should sky really be a free parameter at this point?

      # sub-image, sub-mask, sub-uncd
      subim, subma, subun = ie.trimimage(data[im,:,:,pos], (yr[pos],xr[pos]), 
                                (trim,trim), mask[im,:,:,pos], uncd[im,:,:,pos])

      weights = 1.0 / subun**2.0
      # weights = 0.0 if uncd is huge, give points tiny value
      weights[np.where(weights == 0)] = np.mean(weights)/100.0

      # The parameters to fit:
      pars[0] = fp.y     [pos, im] - yoff[pos] # Y position in the sub-image
      pars[1] = fp.x     [pos, im] - xoff[pos] # X position
      pars[2] = fp.aplev [pos, im]             # aperture flux
      pars[3] = fp.skylev[pos, im]             # sky flux level

      # center of the psf
      psfc = psfyc[pos], psfxc[pos]

      # Fit the PSF position, sky, and aperture flux level to the frame.
      fit, err, niter, csqpdof, stat = pf.psffit(subim, weights, psf[:,:,pos],
                                             pars, psfc, order=order, norm=norm)

      # Save the results of the fit
      fp.oyy       [pos, im] = fit[0] + yoff[pos]
      fp.oyerr     [pos, im] = err[0]
      fp.oxx       [pos, im] = fit[1] + xoff[pos]
      fp.oxerr     [pos, im] = err[1]
      fp.optlev    [pos, im] = fit[2]
      fp.opterr    [pos, im] = err[2]
      fp.oskylev   [pos, im] = fit[3]
      fp.oskyerr   [pos, im] = err[3]
      fp.oniter    [pos, im] = niter
      fp.ocspdof   [pos, im] = csqpdof
      fp.ofitstatus[pos, im] = stat

      subimpsf = pf.scalepsf(psf[:,:,pos], fit, np.shape(subim), psfc, 
                             order, norm)

      impsf = np.zeros((ny, nx))
      impsf = ie.pasteimage(impsf, 
                            (subimpsf-fp.oskylev[pos,im])/fp.optlev[pos,im],
                            (yr[pos],xr[pos]))


      # save yerror appropriately,
      # documentation for curvefit says yerror is "standard error between
      # yfit and y", what is this, really?

      fp.ophotlev[pos,im], fp.ophoterr[pos, im], profile, wght = \
              op.optphot(data[im,:,:,pos] - fp.oskylev[pos,im],
                         impsf, var=uncd[im,:,:,pos]**2 + fp.oskyerr[pos,im]**2,
                         mask=mask[im,:,:,pos])



      le.writelog('\n'                                    +
                  'frame =%4d      ' %im                  + 
                  'pos =%3d         '%pos                 +
                  'ocspdof =%11.3f  '%fp.ocspdof [pos,im] + 
                  'niter   =%5d   '  %fp.oniter  [pos,im] + '\n' +
                  'y     =%8.3f  '   %fp.oyy     [pos,im] + 
                  'yerr =%7.4f    '  %fp.oyerr   [pos,im] + 
                  'flux    =%11.3f  '%fp.ophotlev[pos,im] + 
                  'fluxerr =%9.3f  ' %fp.ophoterr[pos,im] + '\n' +
                  'x     =%8.3f  '   %fp.oxx     [pos,im] + 
                  'xerr =%7.4f    '  %fp.oxerr   [pos,im] + 
                  'skylev  =%11.3f  '%fp.oskylev [pos,im] +
                  'skyerr  =%9.3f   '%fp.oskyerr [pos,im], log)


    # Finding Bad Frames 
    # Make rejection mask:
    #   0 = good;
    #   1 = chisq too high; 
    #  10 = sky-median(sky) too high; 
    # 100 = skyerror too high 
    # It makes more sense if a bitmask, I don't know how to do that just now.

    # HIGH CHISQ
    loc = np.where(fp.ocspdof[pos, 0:nimpos[pos]] >= rejlim[0])
    if np.size(loc) != 0:
      fp.ostatus[pos, :][loc] -= 1
      le.writelog(str(np.size(loc)) + ' frames rejected for high chisq', log)

    # HIGH/LOW SKY RELATIVE TO SLOPE SUBTRACTED MEDIAN
    #FINDME: for MIPS, ignore first frames (code this more cleanly later)
    if (pos == 0) or (pos == 7):
      frames = np.where((np.arange(nimpos[pos]) % 6) != 0)
    else:
      frames = np.where((np.arange(nimpos[pos]) % 5) != 0)

    # Do a linear fit to the sky level as function of time.
    test = np.polyfit((fp.time   [pos, 0:nimpos[pos]])[frames],
                      (fp.oskylev[pos, 0:nimpos[pos]])[frames], 1)

    skyint[pos]   = test[1]
    skyslope[pos] = test[0]
    line          = test[1] + test[0]*fp.time[pos, 0:nimpos[pos]]
    skydiff       = np.abs(fp.oskylev[pos, 0:nimpos[pos]] - line)
    # print( skyint[pos], skyslope[pos] )
    # print( 'skylev', fp[ioskylev, 0:nimpos[pos]-1,pos]

    print('skydiff med: ' + str(np.median(skydiff)) + 
                ', std: ' + str(np.std(skydiff))    )
    loc = np.where(skydiff >= rejlim[1])
    if np.size(loc) != 0: 
      fp.ostatus[pos, :][loc] -= 10
      le.writelog(str(np.size(loc)) + 
                  ' frames rejected for high/low sky level', log)

    # HIGH SKY ERROR
    y = fp.oskyerr[pos, 0:nimpos[pos]]
    loc = np.where(y >= rejlim[2] + np.amin(y))
    if np.size(loc) != 0:
      fp.ostatus[pos, :][loc] -= 100
      le.writelog(str(np.size(loc)) + ' frames rejected for sky error', log)

    print( 'rejection complete')


  fp.ogood[np.where(fp.ostatus == 0)] = 1

  # pixel R position
  fp.orad = np.sqrt( (fp.oxx % 1.0 - 0.5)**2.0 + (fp.oyy % 1.0 - 0.5)**2.0 )

  #le.writelog('optimal photometry time: %f minutes'%((time.time()-tini)/60), log)
  return fp, psf
