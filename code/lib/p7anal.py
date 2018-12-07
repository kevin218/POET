# $Author: kevin $
# $Revision: 554 $
# $Date: 2011-08-25 12:06:35 -0400 (Thu, 25 Aug 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/pipeline/trunk/p7anal.py $
# $Id: p7anal.py 554 2011-08-25 16:06:35Z kevin $

"""
 MODIFICATION HISTORY:
    Written by:	Kevin Stevenson, UCF  	2008-07-02
                kevin218@knights.ucf.edu
    Finished initial version:   kevin   2008-09-08
    Updated for multi events:   kevin   2009-11-01
    Added ip interpolation:     kevin   2010-06-28
    Added ip correlation plots: kevin   2010-08-03
    Added bjdutc and bjdtdb     kevin   2011-03-21
    Added transit calcs         kevin   2011-03-31
    Added kstar and arat calc   kevin   2015-05-19
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import time
import hiloerr
import kurucz_inten
import integrate
import transplan_tb
import orbit
from importlib import reload
reload(transplan_tb)
reload(kurucz_inten)

G      = 6.67428e-11    #m^3/kg/s^2
rhosun = 1408.          #kg/m^3, Wikipedia

def stdanal(event, fit, plotnum, printout, isloadallknots=False):
    
    print('MARK: ' + time.ctime() + ' : Starting Standard Analysis', file=printout)
    obj        = event.eventname
    if   hasattr(fit, 'stepsize'):
        stepsize    = fit.stepsize
    elif hasattr(event.params, 'stepsize'):
        stepsize    = event.params.stepsize
    else:
        stepsize    = 100
    if event.kuruczfile.startswith("/Users/"):  #FINDME: Megan changed this to point to her home directory
        kuruczfile = event.kuruczfile
    else:
        kuruczfile = event.ancildir + "kurucz/" + event.kuruczfile
    if event.filtfile.startswith("/Users/"):    #FINDME: Megan changed this to point to her home directory
        filtfile = event.filtfile
    else:
        filtfile   = event.ancildir + "filter/" + event.filtfile
    numfigs    = 10
    
    #LOAD allparams FOR CURRENT FIT
    try:
        print('Loading ' + fit.allparamsfile)
        allparams = np.load(fit.allparamsfile)
    except:
        allparams = fit.allparams
    
    #PRINT RESIDUALS
    fit.bestlinear  = np.polyfit(fit.phase, fit.residuals, 1)
    print("Linear fit to residuals:", file=printout)
    print("Slope = "  + str(fit.bestlinear[0]), file=printout)
    print("Offset = " + str(fit.bestlinear[1]), file=printout)
    
    #LOAD allknots FROM FILE
    if fit.isipmapping and isloadallknots:
        print('Loading ' + fit.allknotfile + ' into memory.')
        xyshape = fit.xygrid[0].shape
        numit   = allparams.shape[1]
        #fit = event[0].fit[0]
        allknotsflat    = np.memmap(fit.allknotfile, dtype='int16', mode='r', offset=0)
        #RESHAPE INTO 2D ARRAY, IGNORING HEADER INFO FROM EACH WRITE STATEMENT
        #allknots        = (np.reshape(allknotsflat, [numit,-1])[:,40:]) \
        #                    .reshape(numit, xyshape[0], xyshape[1])
        #usedknotsflat   = np.reshape(allknots[np.where(allknots > -30000)], [allparams.shape[1],-1])
        usedknotsflat    = (np.reshape(np.reshape(allknotsflat, [numit,-1])[:,40:], [numit,-1])).T
        #usedknotsflat = usedknotsflat/300000.+1
        #del(allknots)
        del(allknotsflat)
        plt.figure(plotnum*numfigs+700)
        plt.clf()
        plt.hist(usedknotsflat[0,::stepsize]/300000.+1, 20)
        plt.xlabel('Knot 0')
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(plotnum*numfigs+700) + "-" + fit.saveext + ".ps")
        a = plt.suptitle(obj + ' Histogram for Knot 0', size=16)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(plotnum*numfigs+700) + "-" + fit.saveext + ".png")
    
    #COMPUTE CORRELATION COEFFICIENTS FOR ALL FREE PARAMS
    fit.paramcorr = np.corrcoef(allparams[fit.nonfixedpars,:])
    if fit.isipmapping and isloadallknots:
        #fit.knotcorr  = np.corrcoef(usedknotsflat.T)
        print('Calculating correlation coefficients.')
        fit.allcorr   = np.corrcoef(np.concatenate((allparams[fit.nonfixedpars,:], usedknotsflat), axis=0))
        del(usedknotsflat)
    else:
        #fit.knotcorr  = 0
        fit.allcorr   = 0

    #PLOT CORRELATIONS B/W KNOTS AND OTHER FREE PARAMETERS
    if fit.isipmapping == True and event.params.isfixipmap == False and isloadallknots:
        palette = plt.matplotlib.colors.LinearSegmentedColormap('YlOrRd2',plt.cm.datad['YlOrRd'],16384)
        palette.set_under(alpha=0.0, color='w')
        yround = fit.yuc[0] - fit.y[0]
        xround = fit.xuc[0] - fit.x[0]
        xmin = fit.xygrid[0].min()+xround
        xmax = fit.xygrid[0].max()+xround
        ymin = fit.xygrid[1].min()+yround
        ymax = fit.xygrid[1].max()+yround
        #PLOT ECLIPSE
        if hasattr(fit.i, 'depth') and fit.i.depth in fit.nonfixedpars:
            eclcorr = np.zeros(xyshape)
            eclcorr[np.where(fit.binfluxmask.reshape(xyshape) == 1)] = \
                    fit.allcorr[fit.nonfixedpars.size:,np.where(fit.nonfixedpars == fit.i.depth)[0]].flatten()
            plt.figure(plotnum*numfigs+701)
            plt.clf()
            a = plt.imshow(np.abs(eclcorr), cmap=palette, vmin=1e-3, vmax=1, aspect='auto', origin='lower', 
                           extent=(xmin,xmax,ymin,ymax), interpolation='nearest')
            a = plt.colorbar(a, pad=0.05, fraction=0.1)
            a = plt.ylabel('Pixel Position in y', size=14)
            a = plt.xlabel('Pixel Position in x', size=14)
            if ymin < -0.5+yround:
                a = plt.hlines(-0.5+yround, xmin, xmax, 'k')
            if ymax >  0.5+yround:
                a = plt.hlines( 0.5+yround, xmin, xmax, 'k')
            if xmin < -0.5+xround:
                a = plt.vlines(-0.5+xround, ymin, ymax, 'k')
            if xmax >  0.5+xround:
                a = plt.vlines( 0.5+xround, ymin, ymax, 'k')
            plt.savefig(event.modeldir + "/" + obj + "-fig" + str(plotnum*numfigs+701) + "-" + fit.saveext + ".ps")
            a = plt.suptitle(obj + ' Correlation Coefficients for Eclipse Depth vs. Position', size=16)
            plt.savefig(event.modeldir + "/" + obj + "-fig" + str(plotnum*numfigs+701) + "-" + fit.saveext + ".png")
        #PLOT MIDPOINT
        if hasattr(fit.i, 'midpt') and fit.i.midpt in fit.nonfixedpars:
            midptcorr = np.zeros(xyshape)
            midptcorr[np.where(fit.binfluxmask.reshape(xyshape) == 1)] = \
                    fit.allcorr[fit.nonfixedpars.size:,np.where(fit.nonfixedpars == fit.i.midpt)[0]].flatten()
            plt.figure(plotnum*numfigs+702)
            plt.clf()
            a = plt.imshow(np.abs(midptcorr), cmap=palette, vmin=1e-3, vmax=1, aspect='auto', origin='lower', 
                           extent=(xmin,xmax,ymin,ymax), interpolation='nearest')
            a = plt.colorbar(a, pad=0.05, fraction=0.1)
            a = plt.ylabel('Pixel Position in y', size=14)
            a = plt.xlabel('Pixel Position in x', size=14)
            if ymin < -0.5+yround:
                a = plt.hlines(-0.5+yround, xmin, xmax, 'k')
            if ymax >  0.5+yround:
                a = plt.hlines( 0.5+yround, xmin, xmax, 'k')
            if xmin < -0.5+xround:
                a = plt.vlines(-0.5+xround, ymin, ymax, 'k')
            if xmax >  0.5+xround:
                a = plt.vlines( 0.5+xround, ymin, ymax, 'k')
            plt.savefig(event.modeldir + "/" + obj + "-fig" + str(plotnum*numfigs+702) + "-" + fit.saveext + ".ps")
            a = plt.suptitle(obj + ' Correlation Coefficients for Eclipse Midpoint vs. Position', size=16)
            plt.savefig(event.modeldir + "/" + obj + "-fig" + str(plotnum*numfigs+702) + "-" + fit.saveext + ".png")
        
    # calculate 2-sided errors for all parameters FROM MEDIAN
    fit.twosidederr = np.zeros([fit.nump, 2])
    for i in range(len(fit.nonfixedpars)):
        fit.twosidederr[fit.nonfixedpars[i],:] = hiloerr.hiloerr( \
            allparams[fit.nonfixedpars[i],0::10], fit.medianp[fit.nonfixedpars[i],0]) 
    
    #Typical error
    fit.typerr = (fit.twosidederr[:,1] - fit.twosidederr[:,0]) / 2.

    #PHOTON-LIMITED S/N CALCULATION
    fit.photsn = np.sqrt(fit.systemflux * 1e-6 / 1e6 / event.fluxconv * event.exptime * event.gain / \
                         event.srperas / np.mean(event.posscl.T[0][0]) / np.mean(event.posscl.T[0][1]))
    fit.sn     = fit.systemflux / np.sqrt(np.mean(fit.newsigma**2))
    print('Photon-limited S/N = ' + str(fit.photsn), file=printout)
    print('Observed S/N = ' + str(fit.sn), file=printout)
    print('Fraction of optimal S/N = ' + str(fit.sn/fit.photsn), file=printout)
        
    #CALCULATE INTERESTING TRANSIT PARAMETERS
    pi   = np.pi
    p    = event.period
    print('Orbital period = ', p, '+/-', event.perioderr, file=printout)
    if hasattr(fit.i, 'rprs'):
        rprs    = fit.bestp[fit.i.rprs]
        rprsmed = fit.medianp[fit.i.rprs,0]
        rprsall = allparams[fit.i.rprs]
        cosi    = fit.bestp[fit.i.cosi]
        cosimed = fit.medianp[fit.i.cosi,0]
        cosiall = allparams[fit.i.cosi]
        ars     = fit.bestp[fit.i.ars]
        arsmed  = fit.medianp[fit.i.ars,0]
        arsall  = allparams[fit.i.ars]
        #Calculate inclination and errors in degrees
        fit.inc     = np.arccos(cosi)*180./np.pi
        incall      = np.arccos(cosiall)*180./np.pi
        fit.incerr  = hiloerr.hiloerr(incall[::stepsize], np.arccos(cosimed)*180./np.pi)
        #Calculate transit depth and error
        fit.trdepth    = rprs*rprs
        trdepthall     = rprsall*rprsall
        fit.trdeptherr = hiloerr.hiloerr(trdepthall[::stepsize], rprsmed*rprsmed)
        #Calculate impact parameter and errors
        fit.trb = ars*cosi
        medtrb  = arsmed*cosimed
        trball  = arsall*cosiall
        fit.trberr = hiloerr.hiloerr(trball[::stepsize], medtrb)
        #Calculate best-fit transit duration and ingress (in days)
        #Using Eqns. 2+3 from Seager & Mallen-Ornelas (2003)
        fit.trt14 = p/pi*np.arcsin(1/ars*np.sqrt(((1+rprs)**2-(ars*cosi)**2)/(1-cosi**2)))  #total transit (t1-t4)
        fit.trt23 = p/pi*np.arcsin(1/ars*np.sqrt(((1-rprs)**2-(ars*cosi)**2)/(1-cosi**2)))  #(t2-t3)
        fit.trt12 = 0.5*(fit.trt14 - fit.trt23)                                             #Ingress (t1-t2)
        #Calculate median transit duration and ingress
        medtrt14 = p/pi*np.arcsin(1/arsmed*np.sqrt(((1+rprsmed)**2-(arsmed*cosimed)**2)/(1-cosimed**2)))
        medtrt23 = p/pi*np.arcsin(1/arsmed*np.sqrt(((1-rprsmed)**2-(arsmed*cosimed)**2)/(1-cosimed**2)))
        medtrt12 = 0.5*(medtrt14 - medtrt23)
        #Calculate errors relative to median values
        trt14all = p/pi*np.arcsin(1/arsall*np.sqrt(((1+rprsall)**2-(arsall*cosiall)**2)/(1-cosiall**2)))
        trt23all = p/pi*np.arcsin(1/arsall*np.sqrt(((1-rprsall)**2-(arsall*cosiall)**2)/(1-cosiall**2)))
        trt12all = 0.5*(trt14all - trt23all)
        fit.trt14err = hiloerr.hiloerr(trt14all[::stepsize], medtrt14)
        fit.trt12err = hiloerr.hiloerr(trt12all[::stepsize], medtrt12)
        #Calculate stellar density in solar units (rho_s/rho_sun)
        #Using Eqn 9 from Seager & Mallen-Ornelas (2003)
        fit.rhostar = 3*pi/(p*86400)**2/G*ars**3        #MKS Units
        medrhostar  = 3*pi/(p*86400)**2/G*arsmed**3
        rhostarall  = 3*pi/(p*86400)**2/G*arsall**3
        fit.rhostarerr = hiloerr.hiloerr(rhostarall[::stepsize], medrhostar)
        #Calculate stellar mass in solar units (m_s/m_sun)
        #Using Eqn 10 from Seager & Mallen-Ornelas (2003)
        x = 0.8     #valid for F-K star on main sequence
        k = 1.0     #valid for F-K star on main sequence
        fit.mstar = (k**3*fit.rhostar/rhosun)**(1/(1-3*x))
        medmstar  = (k**3*medrhostar/rhosun)**(1/(1-3*x))
        mstarall  = (k**3*rhostarall/rhosun)**(1/(1-3*x))
        fit.mstarerr = hiloerr.hiloerr(mstarall[::10], medmstar)
        #Calculate stellar radius in solar units (r_s/r_sun)
        #Using Eqn 11 from Seager & Mallen-Ornelas (2003)
        fit.rstar = k*fit.mstar**x
        medrstar  = k*medmstar**x
        rstarall  = k*mstarall**x
        fit.rstarerr = hiloerr.hiloerr(rstarall[::10], medrstar)
        #Print results
        print('Transit parameters with two-sided errors for trnlldsp model:', file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Inclination (deg):', fit.inc, fit.incerr[0], \
                                          fit.incerr[1]), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Impact Parameter:', fit.trb, fit.trberr[0], \
                                          fit.trberr[1]), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Transit Depth (%):', fit.trdepth*100, fit.trdeptherr[0]*100, \
                                          fit.trdeptherr[1]*100), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Transit Duration (t1-t4, Hrs):', fit.trt14*24, \
                                          fit.trt14err[0]*24, fit.trt14err[1]*24), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Transit Ingress (t1-t2, Hrs):', fit.trt12*24, \
                                          fit.trt12err[0]*24, fit.trt12err[1]*24), file=printout)
        #Divide density by 1000 to convert units from MKS to CGS.
        print('%32s %8.3f %8.3f %8.3f' % ('Stellar Density (g/cm^3):', (fit.rhostar/1000), \
                                          (fit.rhostarerr[0]/1000), (fit.rhostarerr[1])/1000), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Stellar Mass (M_sun):', fit.mstar, \
                                          fit.mstarerr[0], fit.mstarerr[1]), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Stellar Radius (R_sun):', fit.rstar, \
                                          fit.rstarerr[0], fit.rstarerr[1]), file=printout)
    if hasattr(fit.i, 'rprs2'):
        rprs    = fit.bestp[fit.i.rprs2]
        rprsmed = fit.medianp[fit.i.rprs2,0]
        rprsall = allparams[fit.i.rprs2]
        cosi    = fit.bestp[fit.i.cosi2]
        cosimed = fit.medianp[fit.i.cosi2,0]
        cosiall = allparams[fit.i.cosi2]
        ars     = fit.bestp[fit.i.ars2]
        arsmed  = fit.medianp[fit.i.ars2,0]
        arsall  = allparams[fit.i.ars2]
        #Calculate inclination and errors in degrees
        fit.inc2    = np.arccos(cosi)*180./np.pi
        incall      = np.arccos(cosiall)*180./np.pi
        fit.incerr2 = hiloerr.hiloerr(incall[::stepsize], np.arccos(cosimed)*180./np.pi)
        #Calculate transit depth and errors
        fit.trdepth2    = rprs*rprs
        trdepthall      = rprsall*rprsall
        fit.trdeptherr2 = hiloerr.hiloerr(trdepthall[::stepsize], rprsmed*rprsmed)
        #Calculate impact parameter and errors
        fit.trb2= ars*cosi
        medtrb  = arsmed*cosimed
        trball  = arsall*cosiall
        fit.trberr2 = hiloerr.hiloerr(trball[::10], medtrb)
        #Calculate best-fit transit duration and ingress (in days)
        #Using Eqns. 2+3 from Seager & Mallen-Ornelas (2003)
        fit.trt142 = p/pi*np.arcsin(1/ars*np.sqrt(((1+rprs)**2-(ars*cosi)**2)/(1-cosi**2)))  #total transit (t1-t4)
        fit.trt232 = p/pi*np.arcsin(1/ars*np.sqrt(((1-rprs)**2-(ars*cosi)**2)/(1-cosi**2)))  #(t2-t3)
        fit.trt122 = 0.5*(fit.trt142 - fit.trt232)                                           #Ingress (t1-t2)
        #Calculate median transit duration and ingress
        medtrt14 = p/pi*np.arcsin(1/arsmed*np.sqrt(((1+rprsmed)**2-(arsmed*cosimed)**2)/(1-cosimed**2)))
        medtrt23 = p/pi*np.arcsin(1/arsmed*np.sqrt(((1-rprsmed)**2-(arsmed*cosimed)**2)/(1-cosimed**2)))
        medtrt12 = 0.5*(medtrt14 - medtrt23)
        #Calculate errors relative to median values
        trt14all = p/pi*np.arcsin(1/arsall*np.sqrt(((1+rprsall)**2-(arsall*cosiall)**2)/(1-cosiall**2)))
        trt23all = p/pi*np.arcsin(1/arsall*np.sqrt(((1-rprsall)**2-(arsall*cosiall)**2)/(1-cosiall**2)))
        trt12all = 0.5*(trt14all - trt23all)
        fit.trt14err2 = hiloerr.hiloerr(trt14all[::stepsize], medtrt14)
        fit.trt12err2 = hiloerr.hiloerr(trt12all[::stepsize], medtrt12)
        #Calculate stellar density in solar units (rho_s/rho_sun)
        #Using Eqn 9 from Seager & Mallen-Ornelas (2003)
        fit.rhostar2= 3*pi/(p*86400)**2/G*ars**3        #MKS Units
        medrhostar  = 3*pi/(p*86400)**2/G*arsmed**3
        rhostarall  = 3*pi/(p*86400)**2/G*arsall**3
        fit.rhostarerr2 = hiloerr.hiloerr(rhostarall[::stepsize], medrhostar)
        #Calculate stellar mass in solar units (m_s/m_sun)
        #Using Eqn 10 from Seager & Mallen-Ornelas (2003)
        x = 0.8     #valid for F-K star on main sequence
        k = 1.0     #valid for F-K star on main sequence
        fit.mstar2 = (k**3*fit.rhostar2/rhosun)**(1/(1-3*x))
        medmstar   = (k**3*medrhostar/rhosun)**(1/(1-3*x))
        mstarall   = (k**3*rhostarall/rhosun)**(1/(1-3*x))
        fit.mstarerr2 = hiloerr.hiloerr(mstarall[::stepsize], medmstar)
        #Calculate stellar radius in solar units (r_s/r_sun)
        #Using Eqn 11 from Seager & Mallen-Ornelas (2003)
        fit.rstar2 = k*fit.mstar2**x
        medrstar   = k*medmstar**x
        rstarall   = k*mstarall**x
        fit.rstarerr2 = hiloerr.hiloerr(rstarall[::stepsize], medrstar)
        #Print results
        print('Transit parameters with two-sided errors for trnlldsp2 model:', file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Inclination (deg):', fit.inc2, fit.incerr2[0], \
                                          fit.incerr2[1]), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Impact Parameter:', fit.trb2, \
                                          fit.trberr2[0], fit.trberr2[1]), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Transit Depth (%):', fit.trdepth2*100, fit.trdeptherr2[0]*100, \
                                          fit.trdeptherr2[1]*100), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Transit Duration (t1-t4, Hrs):', fit.trt142*24, \
                                          fit.trt14err2[0]*24, fit.trt14err2[1]*24), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Transit Ingress (t1-t2, Hrs):', fit.trt122*24, \
                                          fit.trt12err2[0]*24, fit.trt12err2[1]*24), file=printout)
        #Divide density by 1000 to convert units from MKS to CGS.
        print('%32s %8.3f %8.3f %8.3f' % ('Stellar Density (g/cm^3):', (fit.rhostar2/1000), \
                                          (fit.rhostarerr2[0]/1000), (fit.rhostarerr2[1])/1000), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Stellar Mass (M_sun):', fit.mstar2, \
                                          fit.mstarerr2[0], fit.mstarerr2[1]), file=printout)
        print('%32s %8.3f %8.3f %8.3f' % ('Stellar Radius (R_sun):', fit.rstar2, \
                                          fit.rstarerr2[0], fit.rstarerr2[1]), file=printout)

    #BRIGHTNESS TEMPERATURE CALCULATION
    # filter
    filterf  = np.loadtxt(filtfile, unpack=True)
    filterf  = np.concatenate((filterf[0:2,::-1].T,[filterf[0:2,0]]))

    # stellar spectrum
    kout = kurucz_inten.read(kuruczfile, freq=True)

    # do the calc
    if hasattr(fit.i, 'depth') or hasattr(fit.i, 'depth2') or hasattr(fit.i, 'depth3'):
        print('Starting Monte-Carlo Temperature Calculation, at ' + time.ctime(), file=printout)
        if hasattr(event.tep, 'g'):
            event.logg    = np.log10(event.tep.g.val*100.)
            event.loggerr = np.log10(event.tep.g.uncert*100.)
        bsdata    = np.zeros((3,event.params.numcalc))
        bsdata[1] = np.random.normal(event.logg,event.loggerr,event.params.numcalc)
        bsdata[2] = np.random.normal(event.tstar,event.tstarerr,event.params.numcalc)
        event.tbntrial = bsdata[1].size
        #mandelecl model
        if hasattr(fit.i, 'depth'):
            bsdata[0] = allparams[int(fit.i.depth),0::int(allparams.shape[1]/event.params.numcalc)]
            tb, tbg, numnegf, fmfreq = calcTb(bsdata, kout, filterf, event)
            print('There were ' + str(numnegf) + ' negative flux values for which a temperature was not computed.', file=printout)
            fit.tbm   = np.median(tb[np.where(tb > 0)])
            fit.tbsd  = np.std(tb[np.where(tb > 0)])
            fit.tbgm  = np.median(tbg[np.where(tbg > 0)])
            fit.tbgsd = np.std(tbg[np.where(tbg > 0)])
            print('Band-center brightness temp = ' + str(round(fit.tbgm,2)) + ' +/- ' + str(round(fit.tbgsd,2)) + ' K for mandelecl model.', file=printout)
            print('Integral    brightness temp = ' + str(round(fit.tbm ,2)) + ' +/- ' + str(round(fit.tbsd, 2)) + ' K for mandelecl model.', file=printout)
        #mandelecl2 model
        if hasattr(fit.i, 'depth2'):
            bsdata[0] = allparams[int(fit.i.depth2),0::int(allparams.shape[1]/event.params.numcalc)]
            tb, tbg, numnegf, fmfreq = calcTb(bsdata, kout, filterf, event)
            print('There were ' + str(numnegf) + ' negative flux values for which a temperature was not computed.', file=printout)
            fit.tbm2   = np.median(tb[np.where(tb > 0)])
            fit.tbsd2  = np.std(tb[np.where(tb > 0)])
            fit.tbgm2  = np.median(tbg[np.where(tbg > 0)])
            fit.tbgsd2 = np.std(tbg[np.where(tbg > 0)])
            print('Band-center brightness temp = ' + str(round(fit.tbgm2,2)) + ' +/- ' + str(round(fit.tbgsd2,2)) + ' K for mandelecl2 model.', file=printout)
            print('Integral    brightness temp = ' + str(round(fit.tbm2 ,2)) + ' +/- ' + str(round(fit.tbsd2, 2)) + ' K for mandelecl2 model.', file=printout)
        #mandelecl3 model
        if hasattr(fit.i, 'depth3'):
            bsdata[0] = allparams[fit.i.depth3,0::allparams.shape[1]/event.params.numcalc]
            tb, tbg, numnegf, fmfreq = calcTb(bsdata, kout, filterf, event)
            print('There were ' + str(numnegf) + ' negative flux values for which a temperature was not computed.', file=printout)
            fit.tbm3   = np.median(tb[np.where(tb > 0)])
            fit.tbsd3  = np.std(tb[np.where(tb > 0)])
            fit.tbgm3  = np.median(tbg[np.where(tbg > 0)])
            fit.tbgsd3 = np.std(tbg[np.where(tbg > 0)])
            print('Band-center brightness temp = ' + str(round(fit.tbgm3,2)) + ' +/- ' + str(round(fit.tbgsd3,2)) + ' K for mandelecl3 model.', file=printout)
            print('Integral    brightness temp = ' + str(round(fit.tbm3 ,2)) + ' +/- ' + str(round(fit.tbsd3, 2)) + ' K for mandelecl3 model.', file=printout)
        
        fit.freq  = fmfreq
        print('Finished Monte-Carlo Temperature Calculation, at ' + time.ctime(), file=printout)
    
    #CALCULATE UTC AND TDB EPHEMERIDES
    if hasattr(event, 'timestd'):
        print('Verify that the ephemeris is reported in BJD_' + event.timestd + '!', file=printout)
        offset = event.bjdtdb.flat[0]-event.bjdutc.flat[0]
        if   event.timestd == 'utc':
            ephtimeutc = event.ephtime
            ephtimetdb = event.ephtime + offset
        elif event.timestd == 'tdb':
            ephtimetdb = event.ephtime
            ephtimeutc = event.ephtime - offset
        else:
            print('Assuming that the ephemeris is reported in BJD_UTC!', file=printout)
            ephtimeutc = event.ephtime
            ephtimetdb = event.ephtime + offset
    else:
        print('Assuming that the ephemeris is reported in BJD_UTC!', file=printout)
        offset = event.bjdtdb.flat[0]-event.bjdutc.flat[0]
        ephtimeutc = event.ephtime
        ephtimetdb = event.ephtime + offset
    print('BJD_TDB - BJD_UTC =', offset*86400, 'seconds', file=printout)
    
    #COMPUTE ECLIPSE TIMES
    if hasattr(fit.i,'midpt'):
        fit.ecltimeutc = (np.floor((event.bjdutc.flat[0] - ephtimeutc)/event.period) + \
                          fit.bestp[fit.i.midpt]) * event.period + ephtimeutc
        fit.ecltimetdb = (np.floor((event.bjdtdb.flat[0] - ephtimetdb)/event.period) + \
                          fit.bestp[fit.i.midpt]) * event.period + ephtimetdb
        fit.ecltimeerr = fit.typerr[fit.i.midpt]*event.period
        print('Eclipse time = ', fit.ecltimeutc, '+/-', fit.ecltimeerr, 'BJD_UTC', file=printout)
        print('Eclipse time = ', fit.ecltimetdb, '+/-', fit.ecltimeerr, 'BJD_TDB', file=printout)
        #COMPUTE MINIMUM ECCENTRICITY
        fit.eminerr    = orbit.error_ecosomega(ephtimetdb, fit.ecltimetdb, event.period,
                                     event.ephtimeerr, fit.ecltimeerr, event.perioderr)
        print('Minimum eccentricity = ', fit.eminerr[0], '+/-', fit.eminerr[1], file=printout)
    if hasattr(fit.i,'midpt2'):
        fit.ecltimeutc2 = (np.floor((event.bjdutc.flat[0] - ephtimeutc)/event.period) + \
                          fit.bestp[fit.i.midpt2]) * event.period + ephtimeutc
        fit.ecltimetdb2 = (np.floor((event.bjdtdb.flat[0] - ephtimetdb)/event.period) + \
                          fit.bestp[fit.i.midpt2]) * event.period + ephtimetdb
        fit.ecltimeerr2 = fit.typerr[fit.i.midpt2]*event.period
        print('Eclipse time 2 = ', fit.ecltimeutc2, '+/-', fit.ecltimeerr2, 'BJD_UTC', file=printout)
        print('Eclipse time 2 = ', fit.ecltimetdb2, '+/-', fit.ecltimeerr2, 'BJD_TDB', file=printout)
        #COMPUTE MINIMUM ECCENTRICITY
        fit.eminerr2    = orbit.error_ecosomega(ephtimetdb, fit.ecltimetdb2, event.period,
                                     event.ephtimeerr, fit.ecltimeerr2, event.perioderr)
        print('Minimum eccentricity = ', fit.eminerr2[0], '+/-', fit.eminerr2[1], file=printout)
    if hasattr(fit.i,'midpt3'):
        fit.ecltimeutc3 = (np.floor((event.bjdutc.flat[0] - ephtimeutc)/event.period) + \
                          fit.bestp[fit.i.midpt3]) * event.period + ephtimeutc
        fit.ecltimetdb3 = (np.floor((event.bjdtdb.flat[0] - ephtimetdb)/event.period) + \
                          fit.bestp[fit.i.midpt3]) * event.period + ephtimetdb
        fit.ecltimeerr3 = fit.typerr[fit.i.midpt3]*event.period
        print('Eclipse time 3 = ', fit.ecltimeutc3, '+/-', fit.ecltimeerr3, 'BJD_UTC', file=printout)
        print('Eclipse time 3 = ', fit.ecltimetdb3, '+/-', fit.ecltimeerr3, 'BJD_TDB', file=printout)
        #COMPUTE MINIMUM ECCENTRICITY
        fit.eminerr3    = orbit.error_ecosomega(ephtimetdb, fit.ecltimetdb3, event.period,
                                     event.ephtimeerr, fit.ecltimeerr3, event.perioderr)
        print('Minimum eccentricity = ', fit.eminerr3[0], '+/-', fit.eminerr3[1], file=printout)
    
    #COMPUTE OBSERVED-CALCULATED (O-C) TRANSIT TIMES
    if hasattr(event, 'timestd') and event.timestd == 'tdb':
        ephtime = ephtimetdb
    else:
        ephtime = ephtimeutc
    if hasattr(fit.i,'trspmid'):
        n = np.round(((fit.bestp[fit.i.trspmid] + event.params.tuoffset) - ephtime)/event.period)
        predictedephtime = ephtime + n*event.period
        #Calc O-C and error in minutes
        fit.oc       = (fit.bestp[fit.i.trspmid] + event.params.tuoffset - predictedephtime)*1440
        fit.ocerr    = fit.typerr[fit.i.trspmid]*1440
        if event.timestd == 'tdb':
            fit.ephdat   = fit.bjdtdb   - predictedephtime
            fit.ephdatuc = fit.bjdtdbuc - predictedephtime
        else:
            fit.ephdat   = fit.bjdutc   - predictedephtime
            fit.ephdatuc = fit.bjdutcuc - predictedephtime
        print('Transit O-C (minutes):', fit.oc, '+/-', fit.ocerr, file=printout)
    if hasattr(fit.i,'trspmid2'):
        n = np.round(((fit.bestp[fit.i.trspmid2] + event.params.tuoffset) - ephtime)/event.period)
        predictedephtime = ephtime + n*event.period
        #Calc O-C and error in minutes
        fit.oc2      = (fit.bestp[fit.i.trspmid2] + event.params.tuoffset - predictedephtime)*1440
        fit.ocerr2   = fit.typerr[fit.i.trspmid2]*1440
        if event.timestd == 'tdb':
            fit.ephdat2   = fit.bjdtdb   - predictedephtime
            fit.ephdatuc2 = fit.bjdtdbuc - predictedephtime
        else:
            fit.ephdat2   = fit.bjdutc   - predictedephtime
            fit.ephdatuc2 = fit.bjdutcuc - predictedephtime
        print('Transit 2 O-C (minutes):', fit.oc2, '+/-', fit.ocerr2, file=printout)
    
    #FINDME: Need to compute light-time correction with true omega and eccentricity
    #orbit.light_time(event.semimaj, 0.0, fit.eminerr[0], event.incl, secondary_primary=False)
    
    print('MARK: ' + time.ctime() + ' : Finished Standard Analysis', file=printout)
    
    return

def calcTb(bsdata, kout, filterf, event):

    kinten, kfreq, kgrav, ktemp, knainten, khead = kout
    ffreq     = event.c / (filterf[:,0] * 1e-6)
    ftrans    = filterf[:,1]
    sz        = event.tbntrial
    tb        = np.zeros(sz)
    tbg       = np.zeros(sz)
    numnegf   = 0		#Number of -ve flux values in allparams
    #guess    = 1        #1: Do not compute integral
    complete  = 0
    fmfreq    = None
    kstar     = kurucz_inten.interp2d(kinten, kgrav, ktemp, bsdata[1], bsdata[2])
    fmstar    = None
    if event.tep.rprssq.val!=-1.0:
        arat      = np.random.normal(event.tep.rprssq.val, event.tep.rprssq.uncert, sz)
    else:
        rprssq = event.tep.rprs.val**2.
        rprssquncert = 2.*event.tep.rprs.val*event.tep.rprs.uncert
        arat      = np.random.normal(rprssq, rprssquncert, sz)
    print('\nComputing brightness temperatures...\n')
    for i in range(sz):
        if bsdata[0,i] > 0:
            fstar   = np.interp(ffreq, kfreq, kstar[i])
            tb[i], tbg[i], fmfreq, fmstar = transplan_tb.transplan_tb(arat[i], bsdata[0,i], bsdata[1,i], bsdata[2,i], kfreq=kfreq, kgrav=kgrav, ktemp=ktemp, kinten=kinten, ffreq=ffreq, ftrans=ftrans, fmfreq=fmfreq, fstar=fstar, fmstar=fmstar)
            #tb[i], tbg[i], fmfreq, fstar, fmstar = transplan_tb.transplan_tb(event.arat, bsdata[0,i], bsdata[1,i], bsdata[2,i], kfreq=kfreq, kgrav=kgrav, ktemp=ktemp, kinten=kinten, ffreq=ffreq, ftrans=ftrans, fmfreq=fmfreq, fstar=fstar, fmstar=fmstar)
        else:
            numnegf += 1
        if (i % (sz / 5) == 0): 
            print(str(complete * 20) + "% complete at " + time.ctime())
            complete += 1
    return tb, tbg, numnegf, fmfreq

