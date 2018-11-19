# $Author: patricio $
# $Revision: 285 $
# $Date: 2010-06-18 17:59:25 -0400 (Fri, 18 Jun 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_5checks.py $
# $Id: poet_5checks.py 285 2010-06-18 21:59:25Z patricio $

#! /usr/bin/env python

"""
  Modificaton History:
  --------------------
  2008-07-02  kevin     Written by Kevin Stevenson, UCF.
                        kevin218@knights.ucf.edu
  2008-09-08  kevin     Finished initial version.
  2009-11-01  kevin     Updated for multi events.
  2010-06-28  kevin     Added ip interpolation.
  2010-07-02  kevin     Added x,y precision calc.
  2010-11-04  patricio  Moved suntime correction here.
                        pcubillos@fulbrightmail.org
"""

import numpy  as np
import matplotlib.pyplot as plt
from astropy.io import fits as pf
import time, os
import manageevent as me
import suntimecorr as stc
import time2phase2 as tp
import logedit     as le
import utc_tt


def checks(eventname, period=None, ephtime=None, cwd=None):
  if cwd == None:
    cwd = os.getcwd()
  os.chdir(cwd)


  # Load the Event
  event = me.loadevent(eventname)

  # Create a log
  oldlogname = event.logname
  logname = event.eventname + "_p5.log"
  log = le.Logedit(logname, oldlogname)
  log.writelog('\nStart Checks: ' + time.ctime())

  # If p5 run after p3: we are using results from PSFfit:
  if not hasattr(event, "phottype"):
    event.phottype = "psffit"
    try:
      os.mkdir("psffit/")
    except:
      pass
    os.chdir("psffit/")

  # Move frame parameters to fit Kevin's syntax:
  # event.fp.param --> event.param
  event.filenames = event.fp.filename
  event.x         = event.fp.x
  event.y         = event.fp.y
  event.sx        = event.fp.sx
  event.sy        = event.fp.sy
  event.time      = event.fp.time
  event.pos       = event.fp.pos
  event.frmvis    = event.fp.frmvis
  event.filename  = event.eventname

  if   event.phottype == "aper":
    event.good      = event.fp.good
    event.aplev     = event.fp.aplev
    event.aperr     = event.fp.aperr
    event.background = event.fp.skylev
    log.writelog('Photometry method is APERTURE')
  elif event.phottype == "psffit":
    event.aplev      = event.fp.psfflux
    event.background = event.fp.psfsky
    # FINDME: do something with aperr and good
    event.aperr      = 0.0025*np.mean(event.fp.psfflux)*(event.aplev*0+1)
    event.good       = np.ones(np.shape(event.aplev))
    log.writelog('Photometry method is PSF FITTING')
  elif event.phottype == "optimal":
    event.good  = event.fp.ogood
    event.aplev = event.fp.ophotlev
    event.aperr = event.fp.ophoterr
    # FINDME: Background from optimal?
    event.background = event.fp.psfsky
    log.writelog('Photometry method is OPTIMAL')


  # UPDATE period AND ephtime
  if period != None:
    event.period     = period[0]
    event.perioderr  = period[1]
  if ephtime != None:
    event.ephtime    = ephtime[0]
    event.ephtimeerr = ephtime[1]

  log.writelog("\nCurrent event = " + event.eventname)
  log.writelog("Kurucz file     = " + event.kuruczfile)
  log.writelog("Filter file     = " + event.filtfile)


  # Light-time correction to BJD:

  # Julian observation date
  #event.juldat = event.jdjf80 + event.fp.time / 86400.0
  event.juldat = event.fp.juldat = event.j2kjd + event.fp.time / 86400.0

  if not event.ishorvec:
    log.writeclose('\nHorizon file not found!')
    return
  print("Calculating BJD correction...")
  event.fp.bjdcor = stc.suntimecorr(event.ra, event.dec, event.fp.juldat,
                                    event.horvecfile)

  # Get bjd times:
  event.bjdcor = event.fp.bjdcor
  #event.bjddat = event.fp.juldat + event.fp.bjdcor / 86400.0
  event.bjdutc = event.fp.juldat + event.fp.bjdcor / 86400.0   # utc bjd date
  event.bjdtdb = np.empty(event.bjdutc.shape)
  for i in range(event.bjdtdb.shape[0]):
    event.bjdtdb[i] = utc_tt.utc_tdb(event.bjdutc[i])   # terrestial bjd date

  # ccampo 3/18/2011: check which units phase should be in
  try:
    if event.tep.ttrans.unit == "BJDTDB":
        event.timestd  = "tdb"
        event.fp.phase = tp.time2phase(event.bjdtdb, event.ephtime,
                                       event.period, event.ecltype)
    else:
        event.timestd  = "utc"
        event.fp.phase = tp.time2phase(event.bjdutc, event.ephtime,
                                       event.period, event.ecltype)
  except:
    event.timestd  = "utc"
    event.fp.phase = tp.time2phase(event.bjdutc, event.ephtime,
                                   event.period, event.ecltype)

  # assign phase variable
  event.phase = event.fp.phase

  # ccampo 3/18/2011: moved this above
  # Eclipse phase, BJD
  #event.fp.phase = tp.time2phase(event.fp.juldat + event.fp.bjdcor / 86400.0,
  #                               event.ephtime, event.period, event.ecltype)

  # verify leapsecond correction
  hfile = event.filenames[0,0]
  try:
    image, event.header = pf.getdata(hfile.decode('utf-8'), header=True)
    dt  = ((event.bjdtdb - event.bjdutc)*86400.0)[0, 0]
    dt2 = event.header['ET_OBS'] - event.header['UTCS_OBS']
    log.writelog('Leap second correction : ' + str(dt) + ' = ' + str(dt2))
  except:
    log.writelog('Could not verify leap-second correction.')

  log.writelog('Min and Max light-time correction: ' +
               np.str(np.amin(event.fp.bjdcor)) + ', ' +
               np.str(np.amax(event.fp.bjdcor)) + ' seconds')

  # Verify light-time correction
  try:
    image, event.header = pf.getdata(hfile.decode('utf-8'), header=True)
    try:
      log.writelog('BJD Light-time correction: ' + str(event.bjdcor[0,0]) +
          ' = ' + str((event.header['BMJD_OBS']-event.header['MJD_OBS'])*86400))
    except:
      log.writelog('HJD Light-time correction: ' + str(event.bjdcor[0,0]) +
          ' = ' + str((event.header['HMJD_OBS']-event.header['MJD_OBS'])*86400))
  except:
    log.writelog('Could not verify light-time correction.')

  # Number of good frames should be > 95%
  log.writelog("Good Frames = %7.3f"%(np.mean(event.good)*100)+ " %")

  log.writelog('\nCentering:     X mean     X stddev  Y mean     Y stddev')
  for pos in np.arange(event.npos):
    log.writelog('position %2d:'%pos+
                 ' %10.5f'%np.mean(event.x[pos, np.where(event.good[pos])]) +
                 ' %9.5f'%np.std(  event.x[pos, np.where(event.good[pos])]) +
                 ' %10.5f'%np.mean(event.y[pos, np.where(event.good[pos])]) +
                 ' %9.5f'%np.std(  event.y[pos, np.where(event.good[pos])]) )

  # COMPUTE RMS POSITION CONSISTENCY
  event.xprecision = np.sqrt(np.median(np.ediff1d(event.x)**2))
  event.yprecision = np.sqrt(np.median(np.ediff1d(event.y)**2))

  log.writelog('RMS of x precision = '    +
               str(np.round(event.xprecision,4)) + ' pixels.')
  log.writelog('RMS of y precision = '    +
               str(np.round(event.yprecision,4)) + ' pixels.')
  if event.phottype == "aper":
    log.writelog('\nCenter & photometry half-width/aperture sizes = ' +
                 str(event.ctrim) + ', ' + str(event.photap) + ' pixels.')
  log.writelog('Period = ' + str(event.period)    + ' +/- ' +
               str(event.perioderr) + ' days')
  log.writelog('Ephemeris = ' + str(event.ephtime)    + ' +/- ' +
               str(event.ephtimeerr) + ' JD')


  fmt1 = ['C0o','C1o','C2o','ro','ko','co','mo','bs','gs','ys','rs','ks','cs','ms']
  fmt2 = ['b,','g,','y,','r,']

  plt.figure(501)
  plt.clf()
  plt.figure(502, figsize=(8,12))
  plt.clf()
  plt.figure(503)
  plt.clf()
  plt.figure(504)
  plt.clf()
  plt.figure(505)
  plt.clf()
  plt.figure(506)
  plt.clf()


  for pos in np.arange(event.npos):
    wheregood = np.where(event.good[pos, :])
    # CHOOSE ONLY GOOD FRAMES FOR PLOTTING
    phase      = event.phase     [pos, :][wheregood]
    aplev      = event.aplev     [pos, :][wheregood]
    jdtime     = event.bjdutc    [pos, :][wheregood]
    background = event.background[pos, :][wheregood]
    # COMPUTE X AND Y PIXEL LOCATION RELATIVE TO ...
    if event.npos > 1:
      # CENTER OF EACH PIXEL
      y = (event.y[pos, :] - np.round(event.y[pos, :]))[wheregood]
      x = (event.x[pos, :] - np.round(event.x[pos, :]))[wheregood]
    else:
      # CENTER OF MEDIAN PIXEL
      y = (event.y[pos, :] - np.round(np.median(event.y)))[wheregood]
      x = (event.x[pos, :] - np.round(np.median(event.x)))[wheregood]

    # SORT aplev BY x, y AND radial POSITIONS
    rad    = np.sqrt(x**2 + y**2)
    xx     = np.sort(x)
    yy     = np.sort(y)
    sxx    = np.sort(event.sx[0])
    syy    = np.sort(event.sy[0])
    rr     = np.sort(rad)
    xaplev = aplev[np.argsort(x)]
    yaplev = aplev[np.argsort(y)]
    raplev = aplev[np.argsort(rad)]

    # BIN RESULTS FOR PLOTTING POSITION SENSITIVITY EFFECT
    nobj      = aplev.size
    nbins     = int(120/event.npos)
    binxx     = np.zeros(nbins)
    binyy     = np.zeros(nbins)
    binsxx    = np.zeros(nbins)
    binsyy    = np.zeros(nbins)
    binrr     = np.zeros(nbins)
    binxaplev = np.zeros(nbins)
    binyaplev = np.zeros(nbins)
    binraplev = np.zeros(nbins)
    binxapstd = np.zeros(nbins)
    binyapstd = np.zeros(nbins)
    binrapstd = np.zeros(nbins)
    binphase  = np.zeros(nbins)
    binaplev  = np.zeros(nbins)
    binapstd  = np.zeros(nbins)
    for i in range(nbins):
        start        = int(1.* i   *nobj/nbins)
        end          = int(1.*(i+1)*nobj/nbins)
        binxx[i]     = np.mean(xx[start:end])
        binyy[i]     = np.mean(yy[start:end])
        binsxx[i]    = np.mean(sxx[start:end])
        binsyy[i]    = np.mean(syy[start:end])
        binrr[i]     = np.mean(rr[start:end])
        binxaplev[i] = np.median(xaplev[start:end])
        binyaplev[i] = np.median(yaplev[start:end])
        binraplev[i] = np.median(raplev[start:end])
        binxapstd[i] = np.std(xaplev[start:end]) / np.sqrt(end-start)
        binyapstd[i] = np.std(yaplev[start:end]) / np.sqrt(end-start)
        binrapstd[i] = np.std(raplev[start:end]) / np.sqrt(end-start)
        binphase[i]  = np.mean(phase[start:end])
        binaplev[i]  = np.median(aplev[start:end])
        binapstd[i]  = np.std(aplev[start:end]) / np.sqrt(end-start)
    
    # PLOT 1: flux
    plt.figure(501)
    plt.errorbar(binphase, binaplev, binapstd, fmt=fmt1[pos],
                 linewidth=1, label=('pos %i'%(pos)))
    plt.title(event.planetname + ' Phase vs. Binned Flux')
    plt.xlabel('Orbital Phase')
    plt.ylabel('Flux')
    plt.legend(loc='best')

    # PLOT 2: position-flux
    plt.figure(502)
    plt.subplot(2,1,1)
    plt.title(event.planetname + ' Position vs. Binned Flux')
    plt.errorbar(binyy, binyaplev, binyapstd, fmt=fmt1[pos],
                 label=('pos %i y'%(pos)))
    plt.ylabel('Flux')
    plt.legend(loc='best')
    plt.subplot(2,1,2)
    plt.errorbar(binxx, binxaplev, binxapstd, fmt=fmt1[pos],
                 label=('pos %i x'%(pos)))
    plt.xlabel('Pixel Postion')
    plt.ylabel('Flux')
    plt.legend(loc='best')

    #PLOT 3: position-phase
    plt.figure(503)

    plt.plot(phase, x, 'b,')
    plt.plot(phase, y, 'r,')
    plt.title(event.planetname + ' Phase vs. Position')
    plt.xlabel('Orbital Phase')
    plt.ylabel('Pixel Position')
    plt.legend('xy')

    #PLOT 4: flux-radial distance
    plt.figure(504)
    plt.errorbar(binrr, binraplev, binrapstd, fmt=fmt1[pos],
                 label=('pos %i'%(pos)))
    plt.title(event.planetname + ' Radial Distance vs. Flux')
    plt.xlabel('Distance From Center of Pixel')
    plt.ylabel('Flux')
    plt.legend(loc='best')

    # ::::::::::: Background setting :::::::::::::::::
    if np.size(background) !=0:
      # number of points per bin:
      npoints = 42
      nbins = int(np.size(background)/npoints)
      medianbg = np.zeros(nbins)
      bphase   = np.zeros(nbins)  # background bin phase
      bintime  = np.zeros(nbins)  # background bin JD time
      for i in range(nbins):
        start        = int(1.0* i   *npoints)
        end          = int(1.0*(i+1)*npoints)
        medianbg[i]  = np.median(background[start:end])
        bphase[i]    = np.mean(  phase     [start:end])
        bintime[i]   = np.mean(  jdtime    [start:end])

      # PLOT 5: background-phase
      day = int(np.floor(np.amin(jdtime)))
      timeunits1 = jdtime  - day
      timeunits2 = bintime - day
      xlabel = 'JD - ' + str(day)
      if event.ecltype == 's':
        timeunits1 = phase
        timeunits2 = bphase
        xlabel = 'Phase'

      plt.figure(505)
      plt.plot(timeunits1, background, color='0.45', linestyle='None',
               marker=',')
      if np.size(background) > 10000:
        plt.plot(timeunits2, medianbg, fmt2[pos], label='median bins')
      plt.title(event.planetname + ' Background level')
      plt.xlabel(xlabel)
      plt.ylabel('Flux')

    # PLOT 6: width-flux
    plt.figure(506)
    plt.subplot(2,1,1)
    plt.title(event.planetname + ' Gaussian Width vs. Binned Flux')
    plt.errorbar(binsyy, binyaplev, binyapstd, fmt=fmt1[pos],
                 label=('width %i y'%(pos)))
    plt.ylabel('Flux')
    plt.legend(loc='best')
    plt.subplot(2,1,2)
    plt.errorbar(binsxx, binxaplev, binxapstd, fmt=fmt1[pos],
                 label=('width %i x'%(pos)))
    plt.xlabel('Gaussian Width')
    plt.ylabel('Flux')
    plt.legend(loc='best')

  figname1 = str(event.eventname) + "-fig501.png"
  figname2 = str(event.eventname) + "-fig502.png"
  figname3 = str(event.eventname) + "-fig503.png"
  figname4 = str(event.eventname) + "-fig504.png"
  figname5 = str(event.eventname) + "-fig505.png"
  figname6 = str(event.eventname) + "-fig506.png"

  plt.figure(501)
  plt.savefig(figname1)
  plt.figure(502)
  plt.savefig(figname2)
  plt.figure(503)
  plt.savefig(figname3)
  plt.figure(504)
  plt.savefig(figname4)
  plt.figure(505)
  plt.plot(timeunits1[0], background[0], color='0.45', linestyle='None',
           marker=',', label='all points')
  plt.legend(loc='best')
  plt.savefig(figname5)
  plt.figure(506)
  plt.savefig(figname6)

  # Saving
  me.saveevent(event, event.eventname + "_p5c")

  cwd = os.getcwd() + "/"
  # Print outputs, end-time, and close log.
  log.writelog("Output files:")
  log.writelog("Data:")
  log.writelog(" " + cwd + event.eventname + "_p5c.dat")
  log.writelog("Log:")
  log.writelog(" " + cwd + logname)
  log.writelog("Figures:")
  log.writelog(" " + cwd + figname1)
  log.writelog(" " + cwd + figname2)
  log.writelog(" " + cwd + figname3)
  log.writelog(" " + cwd + figname4)
  log.writelog(" " + cwd + figname5)
  log.writelog(" " + cwd + figname6)
  log.writeclose('\nEnd Checks: ' + time.ctime())

  return event
