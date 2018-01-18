#! /usr/bin/env python
"""
 MODIFICATION HISTORY:
    Written by:	Kevin Stevenson, UCF  	2008-07-02
                kevin218@knights.ucf.edu
    Finished initial version:   kevin   2008-09-08
    Updated for multi events:   kevin   2009-11-01
    Added ip interpolation:     kevin   2010-06-28
    Added x,y precision calc:   kevin   2010-07-02
"""

import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import os
import pickle
import time
import readeventhdf

def checks(filename, num):
    numfigs = 10

    print('MARK: ' + time.ctime() + ' : Starting Checks')

    #RESTORE HDF5 SAVE FILES FROM IDL
    event = readeventhdf.ReadEventHDF(filename)
    obj   = event.eventname

    print("Current event     = " + obj)
    print("Kurucz file       = " + event.kuruczfile)
    print("Filter file       = " + event.filtfile)

    #IF USING OPTIMAL PHOTOMETRY, MOVE OPHOTLEV TO APLEV
    try:
        event.good  = event.ogood
        event.aplev = event.ophotlev
        event.aperr = event.ophoterr
        print('Photometry method = OPTIMAL')
    except:
        print('Photometry method = APERTURE')

    #Number of good frames should be > 95%
    print("Good Frames [%]   = " + str(np.mean(event.good)*100))
	
    #VERIFY LIGHT-TIME CORRECTION
    #File could be located in /home/jh/ast/esp01/data/ or /home/esp01/data/
    hfile = event.dpref + event.aorname[0] + event.bcddir + event.fpref +           \
            event.aorname[0] + '_' + str(int(event.expid.flat[0])).zfill(4) + '_' + \
            str(int(event.dce.flat[0])).zfill(4) + event.pipev + event.bcdsuf
    try:
        image, event.header = pf.getdata(hfile, header=True)
        print('Light-time correction: ' + str(event.bjdcor.flat[0]) + ' = ' + \
              str((event.header['HMJD_OBS'] - event.header['MJD_OBS'])*86400))
    except:
        print('Could not verify light-time correction.')

    print('Mean, std dev of x center = ' + str(np.mean(event.x)) + ' +/- ' + str(np.std(event.x)))
    print('Mean, std dev of y center = ' + str(np.mean(event.y)) + ' +/- ' + str(np.std(event.y)))
    event.xprecision = [np.mean(np.abs(event.x[:-1] - event.x[1:])), np.std (np.abs(event.x[:-1] - event.x[1:]))]
    event.yprecision = [np.mean(np.abs(event.y[:-1] - event.y[1:])), np.std (np.abs(event.y[:-1] - event.y[1:]))]
    print('Mean, std dev of x precision = ' + 
                    str(np.round(event.xprecision[0],4)) +
          ' +/- ' + str(np.round(event.xprecision[1],4)) + ' pixels.')
    print('Mean, std dev of y precision = ' + 
                    str(np.round(event.yprecision[0],4)) +
          ' +/- ' + str(np.round(event.yprecision[1],4)) + ' pixels.')
    print('Center & photometry aperture sizes = ' + str(event.centap) + ', ' + str(event.photap) + ' pixels.')
    print('Period = ' + str(event.period) + ' +/- ' + str(event.perioderr) + ' days')
    print('Ephemeris = ' + str(event.ephtime) + ' +/- ' + str(event.ephtimeerr) + ' JD')
	
    #CHOOSE ONLY GOOD FRAMES FOR PLOTTING
    phase = event.phase[np.where(event.good == 1)]
    aplev = event.aplev[np.where(event.good == 1)]

    #COMPUTE X AND Y PIXEL LOCATION RELATIVE TO...
    if event.npos > 1:
        #CENTER OF EACH PIXEL
        y = (event.y - np.round(event.y))[np.where(event.good == 1)]
        x = (event.x - np.round(event.x))[np.where(event.good == 1)]
    else:
        #CENTER OF MEDIAN PIXEL
        y = (event.y - np.round(np.median(event.y)))[np.where(event.good == 1)]
        x = (event.x - np.round(np.median(event.x)))[np.where(event.good == 1)]

    #SORT aplev BY x, y AND radial POSITIONS
    rad    = np.sqrt(x**2 + y**2)
    xx     = np.sort(x)
    yy     = np.sort(y)
    rr     = np.sort(rad)
    xaplev = aplev[np.argsort(x)]
    yaplev = aplev[np.argsort(y)]
    raplev = aplev[np.argsort(rad)]

    #BIN RESULTS FOR PLOTTING POSITION SENSITIVITY EFFECT
    nobj      = aplev.size
    nbins     = 120
    binxx     = np.zeros(nbins)
    binyy     = np.zeros(nbins)
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
        start        = int(1.*i*nobj/nbins)
        end          = int(1.*(i+1)*nobj/nbins)
        binxx[i]     = np.mean(xx[start:end])
        binyy[i]     = np.mean(yy[start:end])
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

    #PLOT 
    plt.figure(501+numfigs*num)
    plt.clf()
    #plt.plot(phase, aplev, '.', ms=1)
    plt.errorbar(binphase,binaplev,binapstd,fmt='bo',linewidth=1)
    plt.title(obj + ' Phase vs. Binned Flux')
    plt.xlabel('Orbital Phase')
    plt.ylabel('Flux')
    plt.savefig(str(obj) + "-fig" + str(501+numfigs*num) + ".png")

    #PLOT
    plt.figure(502+numfigs*num, figsize=(8,12))
    plt.clf()
    plt.subplot(2,1,1)
    plt.title(obj + ' Position vs. Binned Flux')
    plt.errorbar(binyy, binyaplev, binyapstd, fmt='ro', label='y')
    #plt.plot(y, aplev, 'r.', ms=1)
    plt.ylabel('Flux')
    plt.legend()
    plt.subplot(2,1,2)
    plt.errorbar(binxx, binxaplev, binxapstd, fmt='bo', label='x')
    #plt.plot(x, aplev, 'b.', ms=1)
    plt.xlabel('Pixel Postion')
    plt.ylabel('Flux')
    plt.legend()
    plt.savefig(str(obj) + "-fig" + str(502+numfigs*num) + ".png")

    #PLOT 
    plt.figure(503+numfigs*num)
    plt.clf()
    plt.plot(phase, x, 'b.', ms=1)
    plt.plot(phase, y, 'r.', ms=1)
    plt.title(obj + ' Phase vs. Position')
    plt.xlabel('Orbital Phase')
    plt.ylabel('Pixel Position')
    plt.legend('xy')
    plt.savefig(str(obj) + "-fig" + str(503+numfigs*num) + ".png")

    #PLOT 
    plt.figure(504+numfigs*num)
    plt.clf()
    #plt.plot(np.sqrt(x**2+y**2), aplev, 'b.', ms=1)
    plt.errorbar(binrr, binraplev, binrapstd, fmt='bo', label='r')
    plt.title(obj + ' Radial Distance vs. Flux')
    plt.xlabel('Distance From Center of Pixel')
    plt.ylabel('Flux')
    plt.legend()
    plt.savefig(str(obj) + "-fig" + str(504+numfigs*num) + ".png")
    
    #PLOT MEANIM
    #plt.figure(104)
    #plt.imshow(event.meanim)
    #plt.title(obj + 'Mean Image')
    #Plot 'x' at mean star center
    #plt.savefig(obj + "-fig104.png")

    plt.show()

    print('MARK: ' + time.ctime() + ' : End p5checks\n')

    return event



 
 
