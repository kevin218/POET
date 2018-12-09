# $Author: kevin $
# $Revision: 536 $
# $Date: 2011-08-05 13:31:43 -0400 (Fri, 05 Aug 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/pipeline/trunk/p9figs.py $
# $Id: p9figs.py 536 2011-08-05 17:31:43Z kevin $

"""
 MODIFICATION HISTORY:
    Written by:	Kevin Stevenson, UCF  	2008-07-02
                kevin218@knights.ucf.edu
    Finished initial version:   kevin   2008-09-08
    Updated for multi events:   kevin   2009-11-01
    Added ip interpolation:     kevin   2010-06-28
    Fixed correlation plots:    kevin   2010-08-01
    Added bjdutc and bjdtdb     kevin   2011-03-21
"""

import numpy as np
import pylab as plt
import time
import os

def figs(event, fitnum, fignum):
    #Define class
    class FixedOrderFormatter(plt.matplotlib.ticker.ScalarFormatter):
        """Formats axis ticks using scientific notation with a constant order of magnitude"""
        def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
            self._order_of_mag = order_of_mag
            ScalarFormatter.__init__(self, useOffset=useOffset, 
                                     useMathText=useMathText)
        def _set_orderOfMagnitude(self, range):
            """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
            self.orderOfMagnitude = self._order_of_mag
    
    print('MARK: ' + time.ctime() + ' : Plot Figures')
    
    obj             = event.eventname
    fit             = event.fit[fitnum]
    if   hasattr(fit, 'stepsize'):
        stepsize    = fit.stepsize
    elif hasattr(event.params, 'stepsize'):
        stepsize    = event.params.stepsize
    else:
        stepsize    = 100
    period          = round(event.period,2)
    nbins           = event.params.nbins
    width           = [8, 8, 8,   8, 8, 8, 8, 8, 8]
    height          = [6, 6, 6, 7.5, 6, 4, 8, 6, 6]
    numfigs         = 10
    #binst, binend   = min(phase), max(phase)
    #binsz           = (binend - binst) / nbins
    #binphase        = np.arange(nbins+1) * binsz + binst
    
    #LOAD allparams FOR CURRENT FIT
    try:
        print('Loading ' + fit.allparamsfile)
        fit.allparams = np.load(fit.allparamsfile)
    except:
        print('Could not load allparams file.')
    
    
    #BIN DATA
    fit.binres      = np.zeros(nbins)
    fit.binresstd   = np.zeros(nbins)
    flatline        = np.zeros(nbins)
    for j in range(nbins):
        start              = int(1.*j*fit.nobj/nbins)
        end                = int(1.*(j+1)*fit.nobj/nbins)
        fit.binres[j]      = np.mean(fit.residuals[start:end])
        fit.binresstd[j]   = np.std(fit.residuals[start:end])
    
    #Assign abscissa time unit (orbits or days)
    if hasattr(event.params, 'timeunit') ==  False:
        event.params.timeunit = 'orbits'
    if event.params.timeunit == 'orbits':
        tuall      = event.phase.flatten()     #Include all frames
        timeunit   = fit.phase       #Use only clipped frames
        timeunituc = fit.phaseuc
        abscissa   = fit.binphase    #Binned clipped frames
        abscissauc = fit.binphaseuc  #Binned unclipped
        xlabel     = 'Orbital Phase (' + str(round(period, 2)) + '-day period)'
    elif event.params.timeunit == 'days-utc':
        tuall      = event.bjdutc.flatten() - event.params.tuoffset
        timeunit   = fit.bjdutc      - event.params.tuoffset
        timeunituc = fit.bjdutcuc    - event.params.tuoffset
        abscissa   = fit.binbjdutc   - event.params.tuoffset
        abscissauc = fit.binbjdutcuc - event.params.tuoffset
        xlabel     = 'BJD_UTC - ' + str(event.params.tuoffset)
    elif event.params.timeunit == 'days-tdb':
        tuall      = event.bjdtdb.flatten() - event.params.tuoffset
        timeunit   = fit.bjdtdb      - event.params.tuoffset
        timeunituc = fit.bjdtdbuc    - event.params.tuoffset
        abscissa   = fit.binbjdtdb   - event.params.tuoffset
        abscissauc = fit.binbjdtdbuc - event.params.tuoffset
        xlabel     = 'BJD_TDB - ' + str(event.params.tuoffset)
    elif event.params.timeunit == 'days':
        tuall      = event.bjdutc.flatten() - event.params.tuoffset
        timeunit   = fit.bjdutc      - event.params.tuoffset
        timeunituc = fit.bjdutcuc    - event.params.tuoffset
        abscissa   = fit.binbjdutc   - event.params.tuoffset
        abscissauc = fit.binbjdutcuc - event.params.tuoffset
        xlabel     = 'BJD - ' + str(event.params.tuoffset)
    
    #PLOT FULL DATA WITH MEDIAN MCMC ECLIPSE
    plt.figure(fignum*numfigs+901, figsize=(width[0],height[0]))
    plt.clf()
    a = plt.axes([0.15,0.10,0.8,0.80])
    plt.plot(timeunituc, fit.fluxuc, 'k.', ms=1, label='Raw Data')
    #plt.plot(abscissa,fit.binmedianfit,'r--', label='Median Fit', lw=2)
    plt.plot(abscissa,fit.binbestfit,'b-', label='Best Fit', lw=2)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlabel(xlabel,size=14)
    plt.ylabel(r'Flux ($\mu Jy$)',size=14)
    plt.legend(loc='best')
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+901) + "-" + fit.saveext + "-full.ps", dpi=300)
    plt.suptitle(event.params.planetname + ' Data With Eclipse Models',size=18)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+901) + "-" + fit.saveext + "-full.png", dpi=300)
    
    #PLOT BINNED AND MEDIAN MCMC ECLIPSE DATA
    plt.figure(fignum*numfigs+902, figsize=(width[1],height[1]))
    plt.clf()
    a = plt.axes([0.15,0.10,0.8,0.80])
    plt.errorbar(abscissauc,fit.binfluxuc,fit.binstduc,fmt='ko',ms=4,linewidth=1, label='Binned Data')
    plt.plot(abscissa,fit.binnoecl,'k-', label='No Eclipse')
    #plt.plot(abscissa,fit.binmedianfit,'r--', label='Median Fit', lw=2)
    plt.plot(abscissa,fit.binbestfit,'b-', label='Best Fit', lw=2)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlabel(xlabel,size=14)
    plt.ylabel(r'Flux ($\mu Jy$)',size=14)
    plt.legend(loc='best')
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+902) + "-" + fit.saveext + "-bin.ps", dpi=300)
    plt.suptitle('Binned ' + event.params.planetname + ' Data With Eclipse Models',size=18)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+902) + "-" + fit.saveext + "-bin.png", dpi=300)
    #plt.text(min(phaseuc), max(binflux), model[0] + ', ' + model[1] + ', ' + model[2])
    
    #PLOT NORMALIZED BINNED AND MEDIAN MCMC ECLIPSE DATA
    plt.figure(fignum*numfigs+903, figsize=(width[2],height[2]))
    plt.clf()
    a = plt.axes([0.15,0.35,0.8,0.55])
    plt.errorbar(abscissauc,fit.normbinflux,fit.normbinsd,fmt='ko',ms=4,linewidth=1, label='Binned Data')
    #plt.plot(abscissa, fit.normbinmedian,'r--', label='Median Fit', lw=2)
    plt.plot(timeunit, fit.normbestfit,'b-', label='Best Fit', lw=2)
    #plt.xticks(size=16)
    plt.setp(a.get_xticklabels(), visible = False)
    plt.yticks(size=14)
    plt.ylabel('Normalized Flux',size=14)
    plt.legend(loc='best')
    xmin, xmax      = plt.xlim()
    plt.axes([0.15,0.1,0.8,0.2])
    #PLOT RESIDUALS WITH BESTFIT LINE
    fit.binresfit   = fit.bestlinear[0]*fit.binphase + fit.bestlinear[1]
    #plt.errorbar(fit.binphase,fit.binres,fit.binresstd,fmt='ko',ms=4,linewidth=1)
    plt.plot(abscissa,fit.binres/fit.mflux,'ko',ms=4)
    plt.plot(abscissa, flatline,'k-',lw=2)
    #plt.plot(fit.binphase, fit.binresfit/fit.mflux,'k-', label='Linear Fit')
    plt.xlim(xmin,xmax)
    #plt.title(event.params.planetname + ' Binned Residuals With Linear Fit')
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlabel(xlabel,size=14)
    plt.ylabel('Residuals',size=14)
    #plt.legend(loc='best')
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+903) + "-" + fit.saveext + "-norm.ps", dpi=300)
    plt.suptitle('Normalized Binned ' + event.params.planetname + ' Data With Eclipse Models',size=18)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+903) + "-" + fit.saveext + "-norm.png", dpi=300)
    
    #MCMC CORRELATION B/W FREE PARAMETERS
    plt.ioff()
    numfp = fit.nonfixedpars.size
    if numfp <= 7:
        k     = 1
        m     = 1
        plt.figure(fignum*numfigs+904, figsize=(width[3],height[3]))
        plt.clf()
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,hspace=0.15,wspace=0.15)
        for i in fit.nonfixedpars[1:numfp]:
            n     = 0
            for j in fit.nonfixedpars[0:numfp-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                if i > j:
                    a = plt.subplot(numfp-1,numfp-1,k)
                    a.set_facecolor(plt.cm.YlOrRd(np.abs(fit.paramcorr[m,n])))
                    if fit.parname[i].startswith('System Flux'):
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if fit.parname[j].startswith('System Flux'):
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if j == fit.nonfixedpars[0]:
                        plt.yticks(size=11)
                        s = fit.parname[i].replace(',','\n')
                        plt.ylabel(s, size = 12)
                    else:
                        a = plt.yticks(visible=False)
                    if i == fit.nonfixedpars[numfp-1]:
                        plt.xticks(size=11,rotation=90)
                        s = fit.parname[j].replace(',','\n')
                        plt.xlabel(s, size = 12)
                    else:
                        a = plt.xticks(visible=False)
                    plt.plot(fit.allparams[j,0::stepsize],fit.allparams[i,0::stepsize],'b,') #,ms=2)
                k += 1
                n +=1
            m +=1
        a = plt.subplot(numfp-1, numfp-1, numfp-1, frameon=False)
        a.yaxis.set_visible(False)
        a.xaxis.set_visible(False)
        a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
        a = plt.text(1.4, 0.5, '|Correlation Coefficients|', rotation='vertical', ha='center', va='center')
        a = plt.colorbar()
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr.png", dpi=300)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr.ps", dpi=300)
    else:
        #MCMC CORRELATION B/W FREE PARAMETERS
        #THIS VERSION SPLITS THE PARAMETERS INTO 3 FIGURES
        #SHOULD BE USED WHEN THE NUMBER OF FREE PARAMETERS IS LARGE (> 7)
        num1   = int(np.ceil((numfp-1)/2))+1
        num2   = int(np.ceil((numfp-1)/2.))+1
        #Part 1
        k     = 1
        m     = 1
        plt.figure(fignum*numfigs+904, figsize=(width[3],height[3]))
        plt.clf()
        for i in fit.nonfixedpars[1:num1]:
            n = 0
            for j in fit.nonfixedpars[0:num1-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                if i > j:
                    a = plt.subplot(num1-1,num1-1,k)
                    a.set_facecolor(plt.cm.YlOrRd(np.abs(fit.paramcorr[m,n])))
                    if fit.parname[i].startswith('System Flux'):
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if fit.parname[j].startswith('System Flux'):
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if j == fit.nonfixedpars[0]:
                        plt.yticks(size=11)
                        s = fit.parname[i].replace(',','\n')
                        plt.ylabel(s, size = 10)
                    else:
                        a = plt.yticks(visible=False)
                    if i == fit.nonfixedpars[num1-1]:
                        plt.xticks(size=11,rotation=90)
                        s = fit.parname[j].replace(',','\n')
                        plt.xlabel(s, size = 10)
                    else:
                        a = plt.xticks(visible=False)
                    plt.plot(fit.allparams[j,0::int(stepsize)],fit.allparams[i,0::int(stepsize)],'b,') #,ms=2)
                k += 1
                n += 1
            m += 1
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,hspace=0.15,wspace=0.15)
        a = plt.subplot(num1-1, num1-1, num1-1, frameon=False)
        a.yaxis.set_visible(False)
        a.xaxis.set_visible(False)
        a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
        a = plt.text(1.4, 0.5, '|Correlation Coefficients|', rotation='vertical', ha='center', va='center')
        a = plt.colorbar()
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr1.png", dpi=300)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr1.ps", dpi=300)
        #Part 2
        k     = 1
        mprime = m
        nprime = n
        plt.figure(fignum*numfigs+9005, figsize=(width[3],height[3]))
        plt.clf()
        for i in fit.nonfixedpars[num1:numfp]:
            n = 0
            for j in fit.nonfixedpars[0:num1-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                a = plt.subplot(num2-1,num1-1,k)
                a.set_facecolor(plt.cm.YlOrRd(np.abs(fit.paramcorr[m,n])))
                if fit.parname[i].startswith('System Flux'):
                    a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                if fit.parname[j].startswith('System Flux'):
                    a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                if j == fit.nonfixedpars[0]:
                    plt.yticks(size=11)
                    s = fit.parname[i].replace(',','\n')
                    plt.ylabel(s, size = 10)
                else:
                    a = plt.yticks(visible=False)
                if i == fit.nonfixedpars[numfp-1]:
                    plt.xticks(size=11,rotation=90)
                    s = fit.parname[j].replace(',','\n')
                    plt.xlabel(s, size = 10)
                else:
                    a = plt.xticks(visible=False)
                plt.plot(fit.allparams[j,0::int(stepsize)],fit.allparams[i,0::int(stepsize)],'b,') #,ms=1)
                k += 1
                n += 1
            m += 1
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,hspace=0.15,wspace=0.15)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr2.png", dpi=300)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr2.ps", dpi=300)
        #Part 3
        k     = 1
        m     = mprime
        plt.figure(fignum*numfigs+9006, figsize=(width[3],height[3]))
        plt.clf()
        for i in fit.nonfixedpars[num1:numfp]:
            n = nprime
            for j in fit.nonfixedpars[num1-1:numfp-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                if i > j:
                    a = plt.subplot(num2-1,num2-1,k)
                    a.set_facecolor(plt.cm.YlOrRd(np.abs(fit.paramcorr[m,n])))
                    if fit.parname[i].startswith('System Flux'):
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if fit.parname[j].startswith('System Flux'):
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if j == fit.nonfixedpars[num1-1]:
                        plt.yticks(size=11)
                        s = fit.parname[i].replace(',','\n')
                        plt.ylabel(s, size = 10)
                    else:
                        a = plt.yticks(visible=False)
                    if i == fit.nonfixedpars[numfp-1]:
                        plt.xticks(size=11,rotation=90)
                        s = fit.parname[j].replace(',','\n')
                        plt.xlabel(s, size = 10)
                    else:
                        a = plt.xticks(visible=False)
                    plt.plot(fit.allparams[int(j),0::int(stepsize)],fit.allparams[int(i),0::int(stepsize)],'b,') #,ms=1)
                k += 1
                n += 1
            m += 1
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,hspace=0.15,wspace=0.15)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr3.png", dpi=300)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+904) + "-" + fit.saveext + "-corr3.ps", dpi=300)
    plt.ion()

    # PLOT RMS vs. BIN SIZE
    plt.figure(fignum*numfigs+905, figsize=(width[4],height[4]))
    plt.clf()
    plt.axes([0.12, 0.12, 0.82, 0.8])
    a = plt.loglog(fit.binsz, fit.rms, color='black', lw=1.5, label='RMS')
    a = plt.loglog(fit.binsz, fit.stderr, color='red', ls='-', lw=2, label='Std. Err.')
    a = plt.xlim(0, fit.binsz[-1]*2)
    a = plt.ylim(fit.rms[-1]/2., fit.rms[0]*2.)
    a = plt.xlabel("Bin Size", fontsize=14)
    a = plt.ylabel("RMS", fontsize=14)
    a = plt.xticks(size=14)
    a = plt.yticks(size=14)
    a = plt.legend(loc='upper right')
    a.fontsize = 8
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+905) + "-" + fit.saveext + "-rms.ps", dpi=300)
    a = plt.suptitle(event.params.planetname + ' Correlated Noise', size=18)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+905) + "-" + fit.saveext + "-rms.png", dpi=300)
    
    #PLOT POSITION SENSITIVITY AND MODELS
    plt.rcParams.update({'legend.fontsize':11})
    plt.figure(906+fignum*numfigs, figsize=(width[5], height[5]))
    plt.clf()
    yround = fit.yuc[0] - fit.y[0]
    xround = fit.xuc[0] - fit.x[0]
    a = plt.subplot(1,2,1)
    a = plt.errorbar(yround+fit.binyy, fit.binyflux, fit.binyflstd, fmt='ro', label='Binned Flux')
    a = plt.plot(yround+fit.binyy, fit.binybestip, 'k-', lw=2, label='BLISS Map')
    a = plt.xlabel('Pixel Postion in y', size=14)
    a = plt.ylabel('Normalized Flux', size=14)
    a = plt.xticks(rotation=90)
    a = plt.legend(loc='best')
    a = plt.subplot(1,2,2)
    a = plt.errorbar(xround+fit.binxx, fit.binxflux, fit.binxflstd, fmt='bo', label='Binned Flux')
    a = plt.plot(xround+fit.binxx, fit.binxbestip, 'k-', lw=2, label='BLISS Map')
    a = plt.xlabel('Pixel Postion in x', size=14)
    a = plt.xticks(rotation=90)
    #a = plt.ylabel('Normalized Flux', size=14)
    a = plt.legend(loc='best')
    plt.subplots_adjust(left=0.11,right=0.97,bottom=0.20,top=0.96,wspace=0.20)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+906) + "-" + fit.saveext + ".ps", dpi=300)
    #a = plt.suptitle('Normalized Binned Flux vs. Position', size=18)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+906) + "-" + fit.saveext + ".png", dpi=300)

    #HISTOGRAM
    plt.ioff()
    j          = 1
    numfp      = fit.ifreepars.size
    histheight = np.min((int(4*np.ceil(numfp/3.)),height[6]))
    if histheight == 4:
        bottom = 0.23
    elif histheight == 8:
        bottom = 0.13
    else:
        bottom = 0.12
    plt.figure(fignum*numfigs+907, figsize=(width[6],histheight))
    plt.clf()
    for i in fit.ifreepars:
        a = plt.subplot(np.ceil(numfp/3.),3,j)
        if fit.parname[i].startswith('System Flux'):
            a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
        plt.xticks(size=12,rotation=90)
        plt.yticks(size=12)
        #plt.axvline(x=fit.meanp[i,0])
        plt.xlabel(fit.parname[i], size=14)
        a  = plt.hist(fit.allparams[int(i),0::int(stepsize)], 20, label=str(fit.meanp[i,0]))
        j += 1
    plt.subplots_adjust(left=0.07,right=0.95,bottom=bottom,top=0.95,hspace=0.40,wspace=0.25)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+907) + "-" + fit.saveext + "-hist.png", dpi=300)
    plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+907) + "-" + fit.saveext + "-hist.ps", dpi=300)
    plt.ion()
    
    #FLUX VS POSITION CONTOUR PLOT OF INTERPOLATED INTRA-PIXEL
    if hasattr(fit, 'binipflux'):
        if fit.model.__contains__('nnint'):
            interp = 'nearest'
        else:
            interp = 'bilinear'
        palette   = plt.matplotlib.colors.LinearSegmentedColormap('jet3',plt.cm.datad['jet'],16384)
        palette.set_under(alpha=0.0, color='w')
        binipflux = fit.binipflux
        vmin = binipflux[np.where(binipflux > 0)].min()
        vmax = binipflux.max()
        yround = fit.yuc[0] - fit.y[0]
        xround = fit.xuc[0] - fit.x[0]
        xmin = fit.xygrid[0].min()+xround
        xmax = fit.xygrid[0].max()+xround
        ymin = fit.xygrid[1].min()+yround
        ymax = fit.xygrid[1].max()+yround
        plt.figure(908+fignum*numfigs, figsize=(width[7], height[7]))
        plt.clf()
        plt.subplots_adjust(left=0.11,right=0.95,bottom=0.10,top=0.90)
        a = plt.imshow(binipflux, cmap=palette, vmin=vmin, vmax=vmax, origin='lower', 
                       extent=(xmin,xmax,ymin,ymax), aspect='auto', interpolation=interp)
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
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+908) + "-" + fit.saveext + \
                    "-fluxContour.ps", dpi=300)
        a = plt.suptitle('BLISS Map', size=18)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+908) + "-" + fit.saveext + \
                    "-fluxContour.png", dpi=300)
        
        #PLOT # OF POINTS/BIN VS POSITION
        lenbinflux = np.zeros(len(fit.wherebinflux))
        for m in range(len(fit.wherebinflux)):
            lenbinflux[m] = len(fit.wherebinflux[m])
        lenbinflux = lenbinflux.reshape(fit.xygrid[0].shape)
        plt.figure(909+fignum*numfigs, figsize=(width[8], height[8]))
        plt.clf()
        plt.subplots_adjust(left=0.11,right=0.95,bottom=0.10,top=0.90)
        a = plt.imshow(lenbinflux, cmap=palette, vmin=1, vmax=lenbinflux.max(), origin='lower', 
                       extent=(xmin,xmax,ymin,ymax), aspect='auto', interpolation=interp)
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
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+909) + "-" + fit.saveext + \
                    "-densityContour.ps", dpi=300)
        a = plt.suptitle('Pointing Histogram', size=18)
        plt.savefig(event.modeldir + "/" + obj + "-fig" + str(fignum*numfigs+909) + "-" + fit.saveext + \
                    "-densityContour.png", dpi=300)
    
    #CHMOD ALL POSTSCRIPT FILES
    for files in os.listdir('.'):
        if files.endswith('.ps'):
            os.chmod(files, 0o664)   #0664 must be in octal
    
    del fit.allparams
    return

