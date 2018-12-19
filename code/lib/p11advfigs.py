
# $Author: kevin $
# $Revision: 549 $
# $Date: 2011-08-24 12:40:35 -0400 (Wed, 24 Aug 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/pipeline/trunk/p11advfigs.py $
# $Id: p11advfigs.py 549 2011-08-24 16:40:35Z kevin $

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import cPickle as pickle
import time
import orbit
import readeventhdf
modeldir = "/home/esp01/ancil/modelparams/"
sys.path.append(modeldir)
sys.path.append('../')
import run
import p11advfigs as p11

def main(args):
    print('MARK: ' + time.ctime() + ' : Bin and Plot')

    # import param files, if necessary
    #obj      = ['hd149bp41', 'hd149bp42', 'hd149bp43']
    obj      = ['hd149bs11','hd149bs21','hd149bs31','hd149bs32','hd149bs33', \
                'hd149bs41','hd149bs42','hd149bs43','hd149bs44','hd149bs51','hd149bs52']
    for i in range(len(obj)):
        exec 'import ' + obj[i] + 'params'

    #RESTORE ALL SAVE FILES
    event = run.p7Restore(filedir='.',topdir='.')
    #event = run.p6Restore(filedir='.',topdir='.')

    '''
    #MANUALLY RESTORE SAVE FILES
    sys.path.append('lib/')
    loadfile = []
    for fname in os.listdir('.'):
        if (fname.endswith("7anal.dat")):
            loadfile.append(fname)

    loadfile.sort()
    numevents = len(loadfile)
    event    = []
    for i in range(5,len(loadfile)):
        print("Loading " + loadfile[i])
        handle      = open(loadfile[i], 'r')
        event.append(pickle.load(handle))
        handle.close()

    '''
    prefix    = 'hd149b'
    #hd149bp41, hd149bp42, hd149bp43 = event     #Edit line manually
    hd149bs11, hd149bs21, hd149bs31, hd149bs32, hd149bs33, \
    hd149bs41, hd149bs42, hd149bs43, hd149bs44, hd149bs51, hd149bs52 = event
    nbins     = [60, 60, 60, 60, 60, 108, 60, 60, 60, 30, 30]
    numevents = len(event)

    #SELECT WHICH MODELS TO PLOT
    for i in range(numevents):
        print(event[i].params.model)

    #fit    = [hd149bp41.fit[0], hd149bp42.fit[0], hd149bp43.fit[0]] #Edit line manually
    #labels = ['p41', 'p42', 'p43']       #Edit line manually
    fit    = [hd149bs11.fit[0], 
              hd149bs21.fit[0], 
              hd149bs31.fit[0],
              hd149bs32.fit[0],
              hd149bs33.fit[0],
              hd149bs41.fit[0],
              hd149bs42.fit[0],
              hd149bs43.fit[0],
              hd149bs44.fit[0],
              hd149bs51.fit[0],
              hd149bs52.fit[0]] #Edit line manually
    labels = ['s11', 's21', 's31', 's32', 's33', 's41', 's42', 's43', 's44', 's51', 's52']       #Edit line manually
    nfit = len(fit)

    # is this a joint fit? If so, specify which models were run jointly. ENTER MANUALLY
    joint     = True
    jointruns = [[0],
                 [1],
                 [2, 3, 4],
                 [5, 6, 7, 8],
                 [9, 10]]

    period          = str(round(event[0].period,3))

    #RUNNING THE CODE BELOW IS OPTIONAL, BASED ON YOUR NEEDS.
    
    #Calculate time from predicted center of **transit** in days
    '''
    for i in range(nfit):
        n = np.round(((fit[i].bestp[fit[i].i.trspmid] + event[i].params.tuoffset) - event[i].ephtime)/event[i].period)
        predictedephtime = event[i].ephtime + n*event[i].period
        fit[i].ephdat   = fit[i].bjdutc   - predictedephtime
        fit[i].ephdatuc = fit[i].bjdutcuc - predictedephtime
    '''
    #CALCULATE WEIGHTED AVERAGE MIDPOINT
    '''
    weight = np.zeros(nfit)
    midpts = np.zeros(nfit)
    for i in np.array((0,1,2,3)):
	    weight[i] = 1/(fit[i].medianp[fit[i].i.midpt,1])**2
	    midpts[i] = weight[i]*fit[i].bestp[fit[i].i.midpt]

    midptsd = np.sqrt(1/sum(weight))
    midptav = sum(midpts)/sum(weight)
    print('Average midpoint = ' + str(round(midptav,5)) + ' +/- ' + str(round(midptsd,5)))

    #CALCULATE WEIGHTED AVERAGE MIDPOINT
    #SPECIFIC FOR SIMULTANEOUS OBSERVATIONS OF CH13 & CH24
    weight = np.zeros(nfit)
    midpts = np.zeros(nfit)
    for i in (0,2):
	    weight[i] = 1/(fit[i].medianp[fit[i].i.midpt,1])**2
	    midpts[i] = weight[i]*fit[i].bestp[fit[i].i.midpt]

    midptsd = np.sqrt(1/sum(weight))
    midptav = sum(midpts)/sum(weight)
    print('Average midpoint for CH13 = ' + str(round(midptav,5)) + ' +/- ' + str(round(midptsd,5)))
    weight = np.zeros(nfit)
    midpts = np.zeros(nfit)
    for i in (1,3):
	    weight[i] = 1/(fit[i].medianp[fit[i].i.midpt,1])**2
	    midpts[i] = weight[i]*fit[i].bestp[fit[i].i.midpt]

    midptsd = np.sqrt(1/sum(weight))
    midptav = sum(midpts)/sum(weight)
    print('Average midpoint for CH24 = ' + str(round(midptav,5)) + ' +/- ' + str(round(midptsd,5)))

    #CALCULATE WEIGHTED AVERAGE TEMPERATURE
    weight = np.zeros(nfit)
    temps  = np.zeros(nfit)
    for i in range(nfit):
	    weight[i] = 1/(fit[i].tbsd)**2
	    temps[i]  = weight[i]*fit[i].tbm

    tempsd = np.sqrt(1/sum(weight))
    tempav = sum(temps)/sum(weight)
    print('Average temperature = ' + str(round(tempav,0)) + ' +/- ' + str(round(tempsd,0)) + ' K')

    #CALCULATE WEIGHTED AVERAGE MINIMUM ECCENTRICITY
    weight = np.zeros(nfit)
    emins  = np.zeros(nfit)
    for i in range(nfit):
	    #print(fit[i].eminerr[0], fit[i].eminerr[1])
	    weight[i] = 1/(fit[i].eminerr[1])**2
	    emins[i]  = weight[i]*np.abs(fit[i].eminerr[0])

    eminsd  = np.sqrt(1/sum(weight))
    eminav  = sum(emins)/sum(weight)
    print('Average minimum eccentricity = ' + str(round(eminav,5)) + ' +/- ' + str(round(eminsd,5)))
    #eminerr = orbit.error_ecosomega(event[0].ephtime, fit.ecltime, event.period,
    #                                event.ephtimeerr, fit.ecltimeerr, event.perioderr)

    #CALCULATE ECCENTRICITY, GIVEN ARGUMENT OF PERIAPSIS
    omega  = np.array((96., 10))*np.pi/180 #WASP-18b
    ecc    = eminav / np.cos(omega[0])
    eccerr = ecc*np.sqrt((eminsd/eminav)**2 + (omega[1]*np.sin(omega[0])/np.cos(omega[0]))**2)
    print('Given that omega = ' + str(round(omega[0],5)) + ' +/- ' + str(round(eccerr,5)) + ' radians...')
    print('Eccentricity = ' + str(round(ecc,5)) + ' +/- ' + str(round(eccerr,5)))
    '''
    #CALCULATE PERIOD
    '''
    time = np.zeros(nfit-1)
    terr = np.zeros(nfit-1)
    j    = 0
    for i in np.array((0,2,3,4,5)):
	    time[j] = fit[i].ecltime
	    terr[j] = fit[i].ecltimeerr
	    j      += 1
	
    model  = lambda period, t0: np.sin(2*np.pi*(time-t0)/period)/terr
    #output = so.leastsq(model, (event.period, time[2]), full_output=True)
    output = so.leastsq(model, event.period, args=(time[2]), ftol=1e-16, xtol=1e-16, full_output=True)
    period = output[0]
    perioderr = np.sqrt(output[1])[0][0]
    '''
    '''
    period    = (fit[0].ecltime - fit[3].ecltime)/81
    perioderr = np.sqrt(fit[0].ecltimeerr**2 + fit[3].ecltimeerr**2)/81
    '''
    #CALCULATE AND APPLY LEAP-SECOND CORRECTION
    '''
    leapsec = np.zeros(nfit)
    for i in range(nfit):
        leapsec[i] = event[i].header['ET_OBS'] - event[i].header['UTCS_OBS']

    print("Leap seconds = " + str(leapsec))
    '''
    
    # The top comment for the table. It should describe the data.
    topstring = 'This file contains the lightcurve from the paper:\n \
                 Analysis of Exoplanet HD 149026b Using BLISS Mapping\n \
                 by Kevin B. Stevenson, Joseph Harrington, Jonathan Fortney,\n \
                 Thomas J. Loredo, Ryan A. Hardy, Sarah Nymeyer, William C. Bowman,\n \
                 Patricio Cubillos, M. Oliver Bowman, and Matthew Hardin\n \
                 which was published in 2011 in ApJ.\n \
                 The data are from the Infrared Array Camera and the Infrared\n \
                 Spectrograph''s photometric blue peak-up array on the US National\n \
                 Aeronautics and Space Administrations Spitzer Space Telescope,\n \
                 programs 254, 40135, 50517, and are available from the public\n \
                 Spitzer archive (http://spitzer.caltech.edu).\n \
                 The paper cited above and its electronic supplement, of which this\n \
                 file is a part, describe the observations and data analysis.\n \
                 The data below are in IPAC table format. The keywords describe the\n \
                 columns in the table. All series (pixel position, frame number,\n \
                 etc.) are zero-based. Data in a table row are valid if GOOD equals 1,\n \
                 invalid otherwise. ORIGIN originally read "Spitzer Science Center",\n \
                 but has been updated to reflect the origin of the present file.'
    
    import irsa as irsa
    #CREATE IRSA TABLES
    #NOTE: if you get an error, simply repeat a few times until it goes away.
    for i in range(nfit):
        irsa.old_irsa(event[i], fit[i], topstring, obj[i]+'_irsa.tbl')
        #irsa.do_irsa(event[i], fit[i], topstring, obj[i]+'_irsa.tbl')
    
    '''
    #PRINT bestfituc FOR LIGHT CURVE FITS FILES & IRSA TABLE
    for i in range(nfit):
        bestfituc = np.zeros(len(fit[i].phaseuc))
        bestfituc[np.where(fit[i].clipmask)] = fit[i].bestfit
        bestfit   = np.zeros(len(event[i].phase.flat))
        order     = np.argsort(event[i].bjdutc,axis=None)
        phase     = event[i].phase.flat[order]
        bjdtdb    = event[i].bjdtdb.flat[order]
        bjdutc    = event[i].bjdutc.flat[order]
        bestfit[np.where(event[i].good.flat)] = bestfituc
        #bestfituc[np.where(event[i].good.flatten() == 1)] = bestgood
        handle = open(event[i].eventname + '-phase.txt','w')
        np.savetxt(handle, phase)
        handle.close()
        handle = open(event[i].eventname + '-bestfit.txt','w')
        np.savetxt(handle, bestfit)
        handle.close()
        #WRITE LEAP-SECOND-CORRECTED bjdtdb AND bjdutc
        handle = open(event[i].eventname + '-bjdtdb.txt','w')
        np.savetxt(handle, bjdtdb)
        handle.close()
        handle = open(event[i].eventname + '-bjdutc.txt','w')
        np.savetxt(handle, bjdutc)
        handle.close()
        #VERIFY DATA HAS BEEN CONVERTED CORRECTLY
        plt.figure(i+1)
        plt.clf()
        plt.suptitle(obj[i])
        plt.plot(phase,bestfit,'.')
        plt.xlabel('Phase')
        plt.ylabel('Flux')
    '''

    #Uncomment for transits
    '''
    for i in range(nfit):
        fit[i].binephdatuc = np.zeros(nbins[i])    
        fit[i].binephdat   = np.zeros(nbins[i])
        for j in range(nbins[i]):
            startuc               = int(1.*j*fit[i].nobjuc/nbins[i])
            enduc                 = int(1.*(j+1)*fit[i].nobjuc/nbins[i])
            fit[i].binphaseuc[j]  = np.mean(fit[i].phaseuc[startuc:enduc])
            fit[i].binephdatuc[j] = np.mean(fit[i].ephdatuc[startuc:enduc])
            start                 = int(1.*j*fit[i].nobj/nbins[i])
            end                   = int(1.*(j+1)*fit[i].nobj/nbins[i])
            fit[i].binephdat[j]   = np.mean(fit[i].ephdat[start:end])
        fit[i].nbestfit    = fit[i].bestfit    / fit[i].mflux
    '''

    #Bin and normalize data
    for i in range(nfit):
        fit[i].binphaseuc = np.zeros(nbins[i])
        fit[i].binfluxuc  = np.zeros(nbins[i])
        fit[i].binstduc   = np.zeros(nbins[i])
        fit[i].normbinfluxuc = np.zeros(nbins[i])
        fit[i].normbinsduc= np.zeros(nbins[i])
        fit[i].binphase   = np.zeros(nbins[i])
        fit[i].binbestfit = np.zeros(nbins[i])
        fit[i].binramp    = np.zeros(nbins[i])
        fit[i].normbinbest = np.zeros(nbins[i])
        for j in range(nbins[i]):
            startuc               = int(1.*j*fit[i].nobjuc/nbins[i])
            enduc                 = int(1.*(j+1)*fit[i].nobjuc/nbins[i])
            fit[i].binphaseuc[j]  = np.mean(fit[i].phaseuc[startuc:enduc])
            fit[i].binfluxuc[j]   = sum(fit[i].fluxuc[startuc:enduc]/fit[i].sigmauc[startuc:enduc]**2)/ \
             		              sum(1/fit[i].sigmauc[startuc:enduc]**2)
            fit[i].binstduc[j]    = np.sqrt(1 / sum(1/fit[i].sigmauc[startuc:enduc]**2))
            fit[i].normbinfluxuc[j] = sum(fit[i].normfluxuc[startuc:enduc] / \
                                    fit[i].normsigmauc[startuc:enduc]**2) / \
                                    sum(1/fit[i].normsigmauc[startuc:enduc]**2)
            fit[i].normbinsduc[j] = np.sqrt(1 / sum(1/fit[i].normsigmauc[startuc:enduc]**2))
            start                 = int(1.*j*fit[i].nobj/nbins[i])
            end                   = int(1.*(j+1)*fit[i].nobj/nbins[i])
            fit[i].binphase[j]    = np.mean(fit[i].phase[start:end]) 
            fit[i].binbestfit[j]  = np.mean(fit[i].bestfit[start:end])
            fit[i].binramp[j]     = np.mean(fit[i].ramp[start:end])
            fit[i].normbinbest[j] = np.mean(fit[i].normbestfit[start:end])
        fit[i].nbinfluxuc  = fit[i].binfluxuc  / fit[i].mflux
        fit[i].nbinstduc   = fit[i].binstduc   / fit[i].mflux
        fit[i].nbinramp    = fit[i].binramp    / fit[i].mflux
        fit[i].nbinbestfit = fit[i].binbestfit / fit[i].mflux
        fit[i].nfluxuc     = fit[i].fluxuc     / fit[i].mflux
        fit[i].normramp    = np.ones(fit[i].binphase.size)

    fmt = ['bo', 'go', 'ro', 'co', 'mo', 'yo']
    plt.rcParams.update({'legend.fontsize':11})
    #PLOT UNBINNED AND BEST MCMC ECLIPSE DATA
    sep             = 0.075             #Edit line manually
    plt.figure(1, figsize=(6, 12))
    plt.clf()
    a = plt.axes([0.13,0.06,0.82,0.88])
    for i in range(nfit):
        #Secondary eclipses
	    a=plt.plot(fit[i].phaseuc, fit[i].nfluxuc-sep*i, 'k,', ms=1)
	    a=plt.plot(fit[i].binphase,fit[i].nbinbestfit-sep*i, '-', lw=2, label=labels[i])
        #Primary transits
	    #plt.plot(fit[i].ephdatuc, fit[i].nfluxuc-sep*i, 'k.', ms=1)
	    #plt.plot(fit[i].binephdat,fit[i].nbinbestfit-sep*i, '-', lw=2, label=labels[i])
    plt.xlim(0.43, 0.55)                #Edit line manually
    plt.ylim(0.2, 1.0501)               #Edit line manually
    plt.legend(loc='lower left')
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=14)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=14)    #Primary transits
    plt.ylabel('Normalized Flux', size=14)
    plt.savefig(prefix+"-raw.ps")
    plt.suptitle('Raw HD 149026b Data With Transit Models', size=16)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=16)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=16)    #Primary transits
    plt.ylabel('Normalized Flux', size=16)
    plt.savefig(prefix+"-raw.png", dpi=300)
    plt.savefig(prefix+"-raw.pdf", dpi=300)

    #PLOT BINNED AND BEST MCMC ECLIPSE DATA
    sep             = 0.014
    plt.figure(2, figsize=(6, 12))
    plt.clf()
    a = plt.axes([0.15,0.06,0.80,0.88])
    for i in range(nfit):
        #Secondary eclipses
	    plt.errorbar(fit[i].binphaseuc,fit[i].nbinfluxuc-sep*i,fit[i].nbinstduc,fmt='ko',ms=4,linewidth=1)
	    #plt.plot(fit[i].binphase,fit[i].nbinramp-sep*i,'k-',lw=2)
	    plt.plot(fit[i].binphase,fit[i].nbinbestfit-sep*i, lw=2, label=labels[i])
        #Primary transits
	    #plt.errorbar(fit[i].binephdatuc,fit[i].nbinfluxuc-sep*i,fit[i].nbinstduc,fmt=fmt[i],ms=4,linewidth=1)
	    #plt.plot(fit[i].binephdat,fit[i].nbinramp-sep*i,'k-',lw=2)
	    #plt.plot(fit[i].binephdat,fit[i].nbinbestfit-sep*i, lw=2, label=labels[i])
    plt.xlim(0.43, 0.55)
    plt.ylim(0.84, 1.01)
    plt.legend(loc='lower left')
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=14)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=14)    #Primary transits
    plt.ylabel('Normalized Flux', size=14)
    plt.savefig(prefix+"-bin.ps")
    plt.suptitle('Binned HD 149026b Data With Transit Models', size=16)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=16)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=16)    #Primary transits
    plt.ylabel('Normalized Flux', size=16)
    plt.savefig(prefix+"-bin.png", dpi=300)
    plt.savefig(prefix+"-bin.pdf", dpi=300)

    #PLOT NORMALIZED BINNED AND BEST MCMC ECLIPSE DATA
    sep             = 0.0025
    plt.figure(3, figsize=(6, 12))
    plt.clf()
    a = plt.axes([0.15,0.06,0.80,0.88])
    for i in range(5):
        #Secondary eclipses
	    plt.errorbar(fit[i].binphaseuc,fit[i].normbinfluxuc-sep*i,
		         fit[i].normbinsduc,fmt='ko',ms=4,linewidth=1)
	    plt.plot(fit[i].phase, fit[i].normbestfit-sep*i,'-', lw=2, label=labels[i])
        #Primary transits
	    #plt.errorbar(fit[i].binephdatuc,fit[i].normbinfluxuc-sep*i,
	    #	     fit[i].normbinsduc,fmt=fmt[i],ms=4,linewidth=1)
	    #plt.plot(fit[i].ephdat, fit[i].normbestfit-sep*i,'-', lw=2, label=labels[i])
    plt.xlim(0.44, 0.55)
    plt.ylim(0.988,1.0005)
    plt.legend(loc='lower left')
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=14)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=14)    #Primary transits
    plt.ylabel('Normalized Flux', size=14)
    plt.savefig(prefix+"-norm1.ps")
    plt.suptitle('Normalized HD 149026b Data With Transit Models', size=16)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=16)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=16)    #Primary transits
    plt.ylabel('Normalized Flux', size=16)
    plt.savefig(prefix+"-norm1.png", dpi=300)
    plt.savefig(prefix+"-norm1.pdf", dpi=300)

    #PLOT NORMALIZED BINNED AND BEST MCMC ECLIPSE DATA
    sep             = 0.0025
    plt.figure(4, figsize=(6, 12))
    plt.clf()
    a = plt.axes([0.15,0.06,0.80,0.88])
    for i in range(5,nfit):
        #Secondary eclipses
	    plt.errorbar(fit[i].binphaseuc,fit[i].normbinfluxuc-sep*(i-5),
		         fit[i].normbinsduc,fmt='ko',ms=4,linewidth=1)
	    plt.plot(fit[i].phase, fit[i].normbestfit-sep*(i-5),'-', lw=2, label=labels[i])
        #Primary transits
	    #plt.errorbar(fit[i].binephdatuc,fit[i].normbinfluxuc-sep*i,
	    #	     fit[i].normbinsduc,fmt=fmt[i],ms=4,linewidth=1)
	    #plt.plot(fit[i].ephdat, fit[i].normbestfit-sep*i,'-', lw=2, label=labels[i])
    plt.xlim(0.44, 0.55)
    plt.ylim(0.984,1.0015)
    plt.legend(loc='lower left')
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=14)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=14)    #Primary transits
    plt.ylabel('Normalized Flux', size=14)
    plt.savefig(prefix+"-norm2.ps")
    plt.suptitle('Normalized HD 149026b Data With Transit Models', size=16)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=16)        #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=16)    #Primary transits
    plt.ylabel('Normalized Flux', size=16)
    plt.savefig(prefix+"-norm2.png", dpi=300)
    plt.savefig(prefix+"-norm2.pdf", dpi=300)
    
    
    #COMBINE EVENTS TO PLOT SINGLE TRANSIT
    ch3 = p11.combineEvents(fit, jointruns[2], 30)      #Edit line manually
    ch4 = p11.combineEvents(fit, jointruns[3], 30)
    ch5 = p11.combineEvents(fit, jointruns[4], numbins=60)
    
    #Plot combined data for CH3 & CH4
    plt.rcParams.update({'legend.fontsize':12})
    sep = 0.0010
    plt.figure(5, figsize=(6, 6))
    plt.clf()
    a = plt.axes([0.16,0.13,0.80,0.80])
    plt.errorbar(ch3[0], ch3[1], ch3[2],fmt='bo',ms=4,linewidth=1)
    plt.errorbar(ch4[0], ch4[1]-sep, ch4[2],fmt='go',ms=4,linewidth=1)
    #Manually edit which fits to include
    plt.plot(fit[3].phase, fit[3].normbestfit,'b-', lw=2, label='5.8')
    plt.plot(fit[5].phase, fit[5].normbestfit-sep,'g-', lw=2, label='8.0')
    plt.xlim(0.46, 0.545)
    plt.ylim(0.9982,1.0005)
    plt.legend(loc='lower right')
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=14)       #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=14)     #Primary transits
    plt.ylabel('Normalized Flux', size=14)
    plt.savefig(prefix+"3b4-comb.ps")
    plt.suptitle('Combined HD 149026b Data With Eclipse Models', size=16)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.xlabel('Orbital Phase (' + period + '-day period)', size=16)       #Secondary eclipses
    #plt.xlabel('Time From Predicted Center of Transit [Days]', size=16)     #Primary transits
    plt.ylabel('Normalized Flux', size=16)
    plt.savefig(prefix+"3b4-comb.png", dpi=300)
    plt.savefig(prefix+"3b4-comb.pdf", dpi=300)

    #CORRELATION PLOTS OF INTERESTING PARAMETERS
    #Create one plot per joint fit
    #Specify interesting parameters
    #Relabel parameter names with event identifier such as CH13 or CH24
    #CH1
    fit[0].interestingpars = np.array([0,2,5,6],dtype=int)
    ch1parname     = fit[0].parname
    ch1parname[ 0] = 's11 Eclipse Phase'
    ch1parname[ 2] = 's11 Eclipse\nFlux Ratio'
    ch1parname[ 5] = 's11 System\nFlux ($\\mu Jy$)'
    ch1parname[ 6] = 's11 Ramp\nLinear Term'
    #CH2
    fit[1].interestingpars = np.array([1,3,6,8,9],dtype=int)
    ch2parname     = fit[1].parname
    ch2parname[ 1] = 's21 Eclipse Phase'
    ch2parname[ 3] = 's21 Eclipse\nFlux Ratio'
    ch2parname[ 6] = 's21 System\nFlux ($\\mu Jy$)'
    ch2parname[ 8] = 's21 Ramp\nr$_0$'
    ch2parname[ 9] = 's21 Ramp\nr$_1$'
    #CH3
    fit[2].interestingpars = np.array([1,3,6,8,9],dtype=int)
    fit[3].interestingpars = np.array([1,3,6,8,9],dtype=int)
    fit[4].interestingpars = np.array([1,3,6,8,9],dtype=int)
    #ch3parname     = np.concatenate((fit[2].parname, fit[3].parname, fit[4].parname))
    s31parname     = fit[2].parname
    s31parname[ 1] = 's31 Eclipse Phase'
    s31parname[ 3] = 's31 Eclipse\nFlux Ratio'
    s31parname[ 6] = 's31 System\nFlux ($\\mu Jy$)'
    s31parname[ 8] = 's31 Ramp\nr$_0$'
    s31parname[ 9] = 's31 Ramp\nr$_1$'

    s32parname     = fit[3].parname
    s32parname[ 1] = 's32 Eclipse Phase'
    s32parname[ 3] = 's32 Eclipse\nFlux Ratio'
    s32parname[ 6] = 's32 System\nFlux ($\\mu Jy$)'
    s32parname[ 8] = 's32 Ramp\nr$_0$'
    s32parname[ 9] = 's32 Ramp\nr$_1$'

    s33parname     = fit[4].parname
    s33parname[ 1] = 's33 Eclipse Phase'
    s33parname[ 3] = 's33 Eclipse\nFlux Ratio'
    s33parname[ 6] = 's33 System\nFlux ($\\mu Jy$)'
    s33parname[ 8] = 's33 Ramp\nr$_0$'
    s33parname[ 9] = 's33 Ramp\nr$_1$'
    '''
    ch3parname = []
    for i in jointruns[2]:
        for j in range(len(fit[i].parname)):
            ch3parname.append(fit[i].parname[j])
    ch3parname[ 0] = 's31 Eclipse Phase'
    #ch3parname[ 2] = 'Eclipse Flux Ratio'
    ch3parname[ 5] = 's31 System\nFlux ($\\mu Jy$)'
    ch3parname[10] = 's32 Eclipse Phase'
    ch3parname[15] = 's32 System\nFlux ($\\mu Jy$)'
    ch3parname[19] = 's33 Eclipse Phase'
    ch3parname[24] = 's33 System\nFlux ($\\mu Jy$)'
    '''
    #CH4
    fit[5].interestingpars = np.array([1,3,6,8,9,20,21],dtype=int)
    s41parname     = fit[5].parname
    s41parname[ 1] = 's41 Eclipse Phase'
    s41parname[ 3] = 's41 Eclipse\nFlux Ratio'
    s41parname[ 6] = 's41 System\nFlux ($\\mu Jy$)'
    s41parname[ 8] = 's41 Ramp\nr$_0$'
    s41parname[ 9] = 's41 Ramp\nr$_1$'
    s41parname[20] = 's41 Visit\nSensitivity v$_0$'
    s41parname[21] = 's41 Visit\nSensitivity v$_1$'
    
    fit[6].interestingpars = np.array([1,3,6,7,8],dtype=int)
    s42parname     = fit[6].parname
    s42parname[ 1] = 's42 Eclipse Phase'
    s42parname[ 3] = 's42 Eclipse\nFlux Ratio'
    s42parname[ 6] = 's42 System\nFlux ($\\mu Jy$)'
    s42parname[ 8] = 's42 Ramp\nr$_2$'
    s42parname[ 7] = 's42 Ramp\nr$_3$'
    
    fit[7].interestingpars = np.array([0,2,5,6],dtype=int)
    s43parname     = fit[7].parname
    s43parname[ 0] = 's43 Eclipse Phase'
    s43parname[ 2] = 's43 Eclipse\nFlux Ratio'
    s43parname[ 5] = 's43 System\nFlux ($\\mu Jy$)'
    s43parname[ 6] = 's43 Ramp\nr$_2$'

    fit[8].interestingpars = np.array([1,3,6,8,9],dtype=int)
    s44parname     = fit[8].parname
    s44parname[ 1] = 's44 Eclipse Phase'
    s44parname[ 3] = 's44 Eclipse\nFlux Ratio'
    s44parname[ 6] = 's44 System\nFlux ($\\mu Jy$)'
    s44parname[ 8] = 's44 Ramp\nr$_0$'
    s44parname[ 9] = 's44 Ramp\nr$_1$'
    
    #CH5
    fit[9].interestingpars = np.array([1,3,6,7,8,11,12,13],dtype=int)
    s51parname     = fit[9].parname
    s51parname[ 1] = 's51 Eclipse Phase'
    s51parname[ 3] = 's51 Eclipse\nFlux Ratio'
    s51parname[ 6] = 's51 System\nFlux ($\\mu Jy$)'
    s51parname[ 8] = 's51 Ramp\nr$_2$'
    s51parname[ 7] = 's51 Ramp\nr$_3$'
    s51parname[11] = 's51 Position\nSensitivity p$_0$'
    s51parname[12] = 's51 Position\nSensitivity p$_1$'
    s51parname[13] = 's51 Position\nSensitivity p$_2$'

    fit[10].interestingpars = np.array([1,3,6,8,9,11,12,13],dtype=int)
    s52parname     = fit[10].parname
    s52parname[ 1] = 's52 Eclipse Phase'
    s52parname[ 3] = 's52 Eclipse\nFlux Ratio'
    s52parname[ 6] = 's52 System\nFlux ($\\mu Jy$)'
    s52parname[ 8] = 's52 Ramp\nr$_0$'
    s52parname[ 9] = 's52 Ramp\nr$_1$'
    s52parname[11] = 's52 Position\nSensitivity p$_0$'
    s52parname[12] = 's52 Position\nSensitivity p$_1$'
    s52parname[13] = 's52 Position\nSensitivity p$_2$'
    
    #Plotting parameters
    #cbarpos=[-0.7, 3.9, 0.8, 0.565, 0.1, 0.3]
    fit = p11.corrplots(fit, [0], fignum=0, label='hd149bs11', parname=ch1parname, stepsize=500)
    fit = p11.corrplots(fit, [1], fignum=1, label='hd149bs21', parname=ch2parname, stepsize=400)
    fit = p11.corrplots(fit, [2], fignum=2, label='hd149bs31', parname=s31parname, stepsize=400)
    fit = p11.corrplots(fit, [3], fignum=3, label='hd149bs32', parname=s32parname, stepsize=400)
    fit = p11.corrplots(fit, [4], fignum=4, label='hd149bs33', parname=s33parname, stepsize=400)
    fit = p11.corrplots(fit, [5], fignum=5, label='hd149bs41', parname=s41parname, stepsize=500)
    fit = p11.corrplots(fit, [6], fignum=6, label='hd149bs42', parname=s42parname, stepsize=1000)
    fit = p11.corrplots(fit, [7], fignum=7, label='hd149bs43', parname=s43parname, stepsize=500)
    fit = p11.corrplots(fit, [8], fignum=8, label='hd149bs44', parname=s44parname, stepsize=1000)
    fit = p11.corrplots(fit, [9], fignum=9, label='hd149bs51', parname=s51parname, stepsize=2000)
    fit = p11.corrplots(fit,[10],fignum=10, label='hd149bs52', parname=s52parname, stepsize=2000)
    
    fit = p11.histplot(fit, [0], fignum=0, label='hd149bs11', parname=ch1parname, stepsize=500)
    fit = p11.histplot(fit, [1], fignum=1, label='hd149bs21', parname=ch2parname, stepsize=400)
    fit = p11.histplot(fit, [2], fignum=2, label='hd149bs31', parname=s31parname, stepsize=400)
    fit = p11.histplot(fit, [3], fignum=3, label='hd149bs32', parname=s32parname, stepsize=400)
    fit = p11.histplot(fit, [4], fignum=4, label='hd149bs33', parname=s33parname, stepsize=400)
    fit = p11.histplot(fit, [5], fignum=5, label='hd149bs41', parname=s41parname, stepsize=500)
    fit = p11.histplot(fit, [6], fignum=6, label='hd149bs42', parname=s42parname, stepsize=1000)
    fit = p11.histplot(fit, [7], fignum=7, label='hd149bs43', parname=s43parname, stepsize=500)
    fit = p11.histplot(fit, [8], fignum=8, label='hd149bs44', parname=s44parname, stepsize=1000)
    fit = p11.histplot(fit, [9], fignum=9, label='hd149bs51', parname=s51parname, stepsize=2000)
    fit = p11.histplot(fit,[10],fignum=10, label='hd149bs52', parname=s52parname, stepsize=2000)
    
    '''
    fit = p11.corrplots(fit, jointruns[0], fignum=0, label='CH1', parname=ch1parname, stepsize=100)
    fit = p11.corrplots(fit, jointruns[1], fignum=1, label='CH2', parname=ch2parname, stepsize=100)
    fit = p11.corrplots(fit, jointruns[2], fignum=2, label='CH3', parname=ch3parname, stepsize=300)
    fit = p11.corrplots(fit, jointruns[3], fignum=3, label='CH4', parname=ch4parname, stepsize=800)
    fit = p11.corrplots(fit, jointruns[4], fignum=4, label='CH5', parname=ch5parname, stepsize=800)
    
    fit = p11.histplot(fit, jointruns[0], fignum=0, label='CH1', parname=ch1parname, stepsize=100)
    fit = p11.histplot(fit, jointruns[1], fignum=1, label='CH2', parname=ch2parname, stepsize=100)
    fit = p11.histplot(fit, jointruns[2], fignum=2, label='CH3', parname=ch3parname, stepsize=300)
    fit = p11.histplot(fit, jointruns[3], fignum=3, label='CH4', parname=ch4parname, stepsize=800)
    fit = p11.histplot(fit, jointruns[4], fignum=4, label='CH5', parname=ch5parname, stepsize=800)
    '''
    #Plot RMS vs. Bin Size
    plt.rcParams.update({'legend.fontsize':10})
    height = 2*np.ceil(nfit/3.0)
    plt.figure(6, figsize=(6,height))
    plt.clf()
    for i in range(nfit):
        a = plt.subplot(np.ceil(nfit/3.0),3,i+1)
        a = plt.loglog(fit[i].binsz, fit[i].rms, color='black', lw=1, label=labels[i])
        a = plt.loglog(fit[i].binsz, fit[i].stderr, color='red', ls='-', lw=2)
        a = plt.xlim(0, fit[i].binsz[-1]*2)
        a = plt.ylim(fit[i].rms[-1]/2., fit[i].rms[0]*2.)
        a = plt.yticks(size=11)
        a = plt.xticks(size=11)
        if i % 3 == 0:
            a = plt.ylabel("RMS", fontsize=11)
        if i >= (nfit-3):
            a = plt.xlabel("Bin Size", fontsize=12)
        a = plt.legend(loc='upper right')
    plt.subplots_adjust(left=0.11,right=0.95,bottom=0.08,top=0.95,hspace=0.20,wspace=0.28)
    plt.savefig(prefix+"s-rms.ps", dpi=300)
    a = plt.suptitle(event[0].params.planetname + ' Correlated Noise', size=16)
    plt.savefig(prefix+"s-rms.png", dpi=300)
    
    
    #CHMOD ALL POSTSCRIPT FILES
    for files in os.listdir('.'):
        if files.endswith('.ps'):
            os.chmod(files, 0664)   #0664 must be in octal
    
    return 0

#HISTOGRAMS
def histplot(fit, fitlist, fignum, label, parname=[], stepsize=100, dim=[8,None]):
    freepars  = []
    allparams = []
    numfp     = 0
    nump      = 0
    for i in fitlist:
        if (hasattr(fit[i], 'allparamsfile') == True) and (hasattr(fit[i], 'allparams') == False):
            #LOAD allparams FOR CURRENT FIT
            try:
                allparamsfile = fit[i].allparamsfile.rpartition('/')[2]
                fit[i].allparams = np.load(allparamsfile)
                print('Loaded ' + allparamsfile)
            except:
                print('Could not load ' + allparamsfile)
                return 0
        if len(allparams) == 0:
            #First iteration
            allparams = fit[i].allparams
            freepars  = fit[i].interestingpars
        else:
            #Remaining iterations
            allparams = np.concatenate((allparams,fit[i].allparams))
            freepars  = np.concatenate((freepars, fit[i].interestingpars + nump))
        parname   = np.concatenate((parname, fit[i].parname))
        numfp    += fit[i].interestingpars.size
        nump     += fit[i].nump
    
    if dim[1] == None:
        dim[1] = int(4*np.ceil(numfp/3.))
        dim[1] = np.min((dim[1],11))
    if dim[1] == 4:
        bottom = 0.21
        hspace = 0.60
    elif dim[1] == 8:
        bottom = 0.15
        hspace = 0.60
    else:
        bottom = 0.15
        hspace = 0.60
    j          = 1
    plt.figure(14+10*fignum, figsize=(dim[0],dim[1]))
    plt.clf()
    for i in freepars:
        a = plt.subplot(np.ceil(numfp/3.),3,j)
        if parname[i].endswith('Flux ($\\mu Jy$)'):
            a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
        if parname[i].find('Position') != -1:
            a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
        plt.xticks(size=12,rotation=90)
        plt.yticks(size=12)
        #plt.axvline(x=fit.meanp[i,0])
        plt.xlabel(parname[i], size=12)
        a  = plt.hist(allparams[i,0::stepsize], 20)
        j += 1
    plt.subplots_adjust(left=0.07,right=0.97,bottom=bottom,top=0.97,hspace=hspace,wspace=0.26)
    plt.savefig(label + "-fig" + str(14+10*fignum) + "-hist.png", dpi=300)
    plt.savefig(label + "-fig" + str(14+10*fignum) + "-hist.ps", dpi=300)
    return fit

#CORRELATION PLOTS OF INTERESTING PARAMETERS
def corrplots(fit, fitlist, fignum, label, parname=[], stepsize=100, dim=[8,7.5], cbarpos=False):
    freepars  = []
    allparams = []
    numfp     = 0
    nump      = 0
    for i in fitlist:
        if (hasattr(fit[i], 'allparamsfile') == True) and (hasattr(fit[i], 'allparams') == False):
            #LOAD allparams FOR CURRENT FIT
            try:
                allparamsfile = fit[i].allparamsfile.rpartition('/')[2]
                fit[i].allparams = np.load(allparamsfile)
                print('Loaded ' + allparamsfile)
            except:
                print('Could not load ' + allparamsfile)
                return 0
        if len(allparams) == 0:
            #First iteration
            allparams = fit[i].allparams
            freepars  = fit[i].interestingpars
        else:
            #Remaining iterations
            allparams = np.concatenate((allparams,fit[i].allparams))
            freepars  = np.concatenate((freepars, fit[i].interestingpars + nump))
        parname   = np.concatenate((parname, fit[i].parname))
        numfp    += fit[i].interestingpars.size
        nump     += fit[i].nump
    
    #Compute correlation coefficients
    paramcorr = np.corrcoef(allparams[freepars,:])
    
    #Create plots
    if numfp <= 6:
        k     = 1
        m     = 1
        plt.figure(10+10*fignum, figsize=(dim[0],dim[1]))
        plt.clf()
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.16,top=0.95,hspace=0.15,wspace=0.15)
        for i in freepars[1:numfp]:
            n     = 0
            for j in freepars[0:numfp-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                if i > j:
                    a = plt.subplot(numfp-1,numfp-1,k)
                    a.set_axis_bgcolor(plt.cm.YlOrRd(np.abs(paramcorr[m,n])))
                    if parname[i].endswith('Flux ($\\mu Jy$)'):
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if parname[j].endswith('Flux ($\\mu Jy$)'):
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if parname[i].find('Position') != -1:
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                    if parname[j].find('Position') != -1:
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                    if j == freepars[0]:
                        plt.yticks(size=11)
                        s = parname[i].replace(',','\n')
                        plt.ylabel(s, size = 12)
                    else:
                        a = plt.yticks(visible=False)
                    if i == freepars[numfp-1]:
                        plt.xticks(size=11,rotation=90)
                        s = parname[j].replace(',','\n')
                        plt.xlabel(s, size = 12)
                    else:
                        a = plt.xticks(visible=False)
                    plt.plot(allparams[j,0::stepsize],allparams[i,0::stepsize],'b,') #,ms=2)
                k += 1
                n +=1
            m +=1
        a = plt.subplot(numfp-1, numfp-1, numfp-1, frameon=False)
        a.yaxis.set_visible(False)
        a.xaxis.set_visible(False)
        a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
        a = plt.text(1.4, 0.5, '|Correlation Coefficients|', rotation='vertical', ha='center', va='center')
        a = plt.colorbar()
        plt.savefig(label + "-fig" + str(10+10*fignum) + "-corr.png", dpi=300)
        plt.savefig(label + "-fig" + str(10+10*fignum) + "-corr.ps", dpi=300)
    else:
        #MCMC CORRELATION B/W FREE PARAMETERS
        #THIS VERSION SPLITS THE PARAMETERS INTO 3 FIGURES
        #SHOULD BE USED WHEN THE NUMBER OF FREE PARAMETERS IS LARGE (> 7)
        num1   = int(np.ceil((numfp-1)/2))+1
        num2   = int(np.ceil((numfp-1)/2.))+1
        #Part 1
        k     = 1
        m     = 1
        plt.figure(11+10*fignum, figsize=(dim[0],dim[1]))
        plt.clf()
        for i in freepars[1:num1]:
            n = 0
            for j in freepars[0:num1-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                if i > j:
                    a = plt.subplot(num1-1,num1-1,k)
                    a.set_axis_bgcolor(plt.cm.YlOrRd(np.abs(paramcorr[m,n])))
                    if parname[i].endswith('Flux ($\\mu Jy$)'):
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if parname[j].endswith('Flux ($\\mu Jy$)'):
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if parname[i].find('Position') != -1:
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                    if parname[j].find('Position') != -1:
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                    if j == freepars[0]:
                        plt.yticks(size=11)
                        s = parname[i].replace(',','\n')
                        plt.ylabel(s, size = 12)
                    else:
                        a = plt.yticks(visible=False)
                    if i == freepars[num1-1]:
                        plt.xticks(size=11,rotation=90)
                        s = parname[j].replace(',','\n')
                        plt.xlabel(s, size = 12)
                    else:
                        a = plt.xticks(visible=False)
                    plt.plot(allparams[j,0::stepsize],allparams[i,0::stepsize],'b,') #,ms=2)
                k += 1
                n += 1
            m += 1
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,hspace=0.15,wspace=0.15)
        if cbarpos == False:
            a = plt.subplot(num1-1, num1-1, num1-1, frameon=False)
            a.yaxis.set_visible(False)
            a.xaxis.set_visible(False)
            a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
            a = plt.text(1.4, 0.5, '|Correlation Coefficients|', rotation='vertical', ha='center', va='center')
            a = plt.colorbar()
        else:
            a = plt.subplot(numfp-1, numfp-1, 2*(numfp-1), frameon=False)
            a.yaxis.set_visible(False)
            a.xaxis.set_visible(False)
            a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
            a = plt.text(cbarpos[0], cbarpos[1], '|Correlation Coefficients|', rotation='vertical', ha='center', va='center', size=12)
            a = plt.axes(cbarpos[2:6],visible=False)
            a.yaxis.set_visible(False)
            a.xaxis.set_visible(False)
            a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
            a = plt.colorbar()
        
        plt.savefig(label + "-fig" + str(11+10*fignum) + "-corr1.png", dpi=300)
        plt.savefig(label + "-fig" + str(11+10*fignum) + "-corr1.ps", dpi=300)
        #Part 2
        k     = 1
        mprime = m
        nprime = n
        plt.figure(12+10*fignum, figsize=(dim[0],dim[1]))
        plt.clf()
        for i in freepars[num1:numfp]:
            n = 0
            for j in freepars[0:num1-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                a = plt.subplot(num2-1,num1-1,k)
                a.set_axis_bgcolor(plt.cm.YlOrRd(np.abs(paramcorr[m,n])))
                if parname[i].endswith('Flux ($\\mu Jy$)'):
                    a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                if parname[j].endswith('Flux ($\\mu Jy$)'):
                    a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                if parname[i].find('Position') != -1:
                    a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                if parname[j].find('Position') != -1:
                    a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                if j == freepars[0]:
                    plt.yticks(size=11)
                    s = parname[i].replace(',','\n')
                    plt.ylabel(s, size = 12)
                else:
                    a = plt.yticks(visible=False)
                if i == freepars[numfp-1]:
                    plt.xticks(size=11,rotation=90)
                    s = parname[j].replace(',','\n')
                    plt.xlabel(s, size = 12)
                else:
                    a = plt.xticks(visible=False)
                plt.plot(allparams[j,0::stepsize],allparams[i,0::stepsize],'b,') #,ms=1)
                k += 1
                n += 1
            m += 1
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,hspace=0.15,wspace=0.15)
        plt.savefig(label + "-fig" + str(12+10*fignum) + "-corr2.png", dpi=300)
        plt.savefig(label + "-fig" + str(12+10*fignum) + "-corr2.ps", dpi=300)
        #Part 3
        k     = 1
        m     = mprime
        plt.figure(13+10*fignum, figsize=(dim[0],dim[1]))
        plt.clf()
        for i in freepars[num1:numfp]:
            n = nprime
            for j in freepars[num1-1:numfp-1]:
                #plt.suptitle(event.params.planetname + ' Correlation Between Free Parameters',size=24)
                if i > j:
                    a = plt.subplot(num2-1,num2-1,k)
                    a.set_axis_bgcolor(plt.cm.YlOrRd(np.abs(paramcorr[m,n])))
                    if parname[i].endswith('Flux ($\\mu Jy$)'):
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if parname[j].endswith('Flux ($\\mu Jy$)'):
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.0f'))
                    if parname[i].find('Position') != -1:
                        a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                    if parname[j].find('Position') != -1:
                        a.xaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.4f'))
                    if j == freepars[num1-1]:
                        plt.yticks(size=11)
                        s = parname[i].replace(',','\n')
                        plt.ylabel(s, size = 12)
                    else:
                        a = plt.yticks(visible=False)
                    if i == freepars[numfp-1]:
                        plt.xticks(size=11,rotation=90)
                        s = parname[j].replace(',','\n')
                        plt.xlabel(s, size = 12)
                    else:
                        a = plt.xticks(visible=False)
                    plt.plot(allparams[j,0::stepsize],allparams[i,0::stepsize],'b,') #,ms=1)
                k += 1
                n += 1
            m += 1
        plt.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.95,hspace=0.15,wspace=0.15)
        if cbarpos == False:
            a = plt.subplot(num1-1, num1-1, num1-1, frameon=False)
            a.yaxis.set_visible(False)
            a.xaxis.set_visible(False)
            a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
            a = plt.text(1.4, 0.5, '|Correlation Coefficients|', rotation='vertical', ha='center', va='center')
            a = plt.colorbar()
        else:
            a = plt.subplot(numfp-1, numfp-1, 2*(numfp-1), frameon=False)
            a.yaxis.set_visible(False)
            a.xaxis.set_visible(False)
            a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
            a = plt.text(cbarpos[0], cbarpos[1], '|Correlation Coefficients|', rotation='vertical', ha='center', va='center', size=12)
            a = plt.axes(cbarpos[2:6],visible=False)
            a.yaxis.set_visible(False)
            a.xaxis.set_visible(False)
            a = plt.imshow([[0,1],[0,0]], cmap=plt.cm.YlOrRd, visible=False)
            a = plt.colorbar()
        plt.savefig(label + "-fig" + str(13+10*fignum) + "-corr3.png", dpi=300)
        plt.savefig(label + "-fig" + str(13+10*fignum) + "-corr3.ps", dpi=300)
    
    return fit

#COMBINE EVENTS TO PLOT SINGLE TRANSIT/ECLIPSE
def combineEvents(fit, fitlist, numbins, eclipse=True):
    time     = []
    timeuc   = []
    flux     = []
    fluxuc   = []
    sigma    = []
    sigmauc  = []
    for i in fitlist:
        if eclipse == True:
            time     = np.concatenate((time,    fit[i].phase))       #Secondary eclipses
            timeuc   = np.concatenate((timeuc,  fit[i].phaseuc))
        else:
            time     = np.concatenate((time,    fit[i].ephdat))    #Primary transits
            timeuc   = np.concatenate((timeuc,  fit[i].ephdatuc))
        flux     = np.concatenate((flux,    fit[i].normflux))
        fluxuc   = np.concatenate((fluxuc,  fit[i].normfluxuc))
        sigma    = np.concatenate((sigma,   fit[i].normsigma))
        sigmauc  = np.concatenate((sigmauc, fit[i].normsigmauc))
    
    #Sort in time
    args     = np.argsort(time)
    argsuc   = np.argsort(timeuc)
    time     = time[args]
    timeuc   = timeuc[argsuc]
    flux     = flux[args]
    fluxuc   = fluxuc[argsuc]
    sigma    = sigma[args]
    sigmauc  = sigmauc[argsuc]
    
    #Bin
    nobjuc        = len(fluxuc)
    bintimeuc   = np.zeros(numbins)
    normbinfluxuc = np.zeros(numbins)
    normbinsduc   = np.zeros(numbins)
    for j in range(numbins):
        startuc             = int(1.*j*nobjuc/numbins)
        enduc               = int(1.*(j+1)*nobjuc/numbins)
        bintimeuc[j]      = np.mean(timeuc[startuc:enduc])
        normbinfluxuc[j]    = sum(fluxuc[startuc:enduc] / \
                                sigmauc[startuc:enduc]**2) / \
                                sum(1/sigmauc[startuc:enduc]**2)
        normbinsduc[j]      = np.sqrt(1 / sum(1/sigmauc[startuc:enduc]**2))
    return (bintimeuc, normbinfluxuc, normbinsduc)

if __name__ == '__main__':
    sys.exit(main(sys.argv))

