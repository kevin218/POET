'''
# $Author: kevin $
# $Revision: 551 $
# $Date: 2011-08-24 13:51:07 -0400 (Wed, 24 Aug 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_denoise.py $
# $Id: poet_denoise.py 551 2011-08-24 17:51:07Z kevin $
'''

'''
import sys, os
r = os.getcwd().split("/")
maindir = "/".join(r[:r.index("run")])
sys.path.append(maindir + '/lib/')
'''
import numpy  as np
#import pyfits as pf
from astropy.io import fits as pf
import matplotlib.pyplot as plt
import sys, time, os, shutil, copy
import reader3      as rd
import logedit      as le
import manageevent  as me
import centerdriver as cd
import imageedit    as ie
import timer        as t
import multiprocessing as mp
import pywt


#This function performs image denoising using BayesShrink soft thresholding.
#Ref: Chang et al. "Adaptive Wavelet Thresholding for Image Denoising and Compression", 2000
def bayesshrink(frames, wavelet, numlvls, loc):
    sys.stdout.write('\r'+wavelet+' '+str(loc[2])+' '+str(loc[1])+' '+str(loc[0])+'     ')
    sys.stdout.flush()
    nframes = len(frames)
    #Perform wavelet decomposition
    dec = pywt.wavedec(frames,wavelet)
    #Estimate noise variance
    noisevar = np.inf
    for i in range(-1,-numlvls-1,-1):
        noisevar = np.min([(np.median(np.abs(dec[i]))/0.6745)**2,noisevar])
        #noisevar = (np.median(np.abs(dec[-1]))/0.6745)**2
    #At each level of decomposition...
    for i in range(-1,-numlvls-1,-1):
        #Estimate variance at level i then compute the threshold value
        sigmay2 = np.mean(dec[i]*dec[i])
        sigmax  = np.sqrt(np.max([sigmay2-noisevar,0]))
        if sigmax == 0:
            threshold = np.max(np.abs(dec[i]))
        else:
            threshold = noisevar/sigmax
        #Compute less noisy coefficients by applying soft thresholding
        dec[i] = map (lambda x: pywt.thresholding.soft(x,threshold), dec[i])

    frames = pywt.waverec(dec,wavelet)[:nframes]
    return [frames, loc]

# Update event.data with denoised values
def writedata(arg):
    [frames, [j,i,pos]] = arg
    (event.data[:,j,i,pos])[np.where(event.mask[:,j,i,pos])] = frames
    return

# Plot histogram of wavelet coefficients at various levels
def histwc(event, wavelet, numlvls=4, pos=0, log=None, denoised=True, ylim=None):
    if denoised == True:
        fignum = 1
        title = "After Denoising"
    else:
        fignum = 0
        title  = "Before Denoising"
    # Perform wavelet decomposition on pixel srcest
    ipixel  = np.array((int(event.srcest[0,pos]),int(event.srcest[1,pos])))
    pixflux = event.data[:,ipixel[0],ipixel[1],pos][np.where(event.mask[:,ipixel[0],ipixel[1],pos])]
    dec = pywt.wavedec(pixflux,wavelet)
    if denoised == False and log != None:
        #Estimate noise variance
        noisevar = np.inf
        for i in range(-1,-numlvls-1,-1):
            noisevar = np.min([(np.median(np.abs(dec[i]))/0.6745)**2,noisevar])
        log.writelog("Denoising results for pixel " + str(ipixel[0]) +","+ str(ipixel[1]) +":")
        log.writelog("Noise uncertainty: " + str(int(np.sqrt(noisevar))))
    # Plot histogram of wavelet coefficients at pixel srcest
    plt.figure(300 + 10*pos + fignum)
    plt.clf()
    for i in range(-1,-numlvls-1,-1):
        # Estimate variance at level i assuming zero mean
        sigmay2 = np.mean(dec[i]*dec[i])
        if denoised == False and log != None:
            # then compute the threshold value
            sigmax  = np.sqrt(np.max([sigmay2-noisevar,0]))
            if sigmax == 0:
                threshold = np.max(np.abs(dec[i]))
            else:
                threshold = noisevar/sigmax
            log.writelog("Level " +str(-i)+ " threshold: " +str(int(threshold)))
        # Set hist range to ignore outliers
        meanwc  = np.mean(dec[i])
        stdwc   = np.std(dec[i])
        wcrange = (meanwc-3*stdwc, meanwc+3*stdwc)
        #log.writelog("Pos " +str(pos)+ ", pixel " + str(ipixel[0]) +','+ str(ipixel[1]) + ":")
        #log.writelog("  Mean pixel flux: " + str(int(np.round(meanpixflux))))
        #log.writelog("  Standard deviation: " + str(int(np.round(stdpixflux))))
        a = plt.hist(dec[i], bins=20, range=wcrange, histtype='step', lw=2,
                     label="L"+str(-i)+" $\sigma$="+str(int(np.round(np.sqrt(sigmay2)))))
        plt.suptitle(event.denoisedir + ' Histogram at Pixel ' + str(ipixel[0]) +','+ str(ipixel[1]))
        plt.xlabel('Values of Wavelet Coefficients')
        plt.ylabel('Count')
        if ylim == None:
            ylim = plt.ylim()
        else:
            plt.ylim(ylim)
        plt.title(title, size=11)
        plt.legend(loc='best')
        plt.savefig("fig" + str(300+10*pos+fignum) + ".png")

    return ylim

# Plot first 'length' frames of lightcurve at pixel srcest
def plotlc(event, pos=0, length=200, denoised=True):
    # Obtain frames at pixel srcest
    ipix0, ipix1  = int(event.srcest[0,pos]),int(event.srcest[1,pos])
    frames = (event.data[:,ipix0,ipix1,pos])[np.where(event.mask[:,ipix0,ipix1,pos])]
    plt.figure(302 + 10*pos)
    if denoised == False:
        plt.clf()
        #label  = "None"
        plt.plot(frames[:length], 'bo', label="None")
    else:
        #label  = "L"+str(event.numlvls)
        plt.plot(frames[:length], 'sg-', ms=4, lw=2, label="L"+str(event.numlvls))
        plt.suptitle(event.denoisedir + ' Lightcurve at Pixel ' + str(ipix0) +','+ str(ipix1))
        plt.xlabel('Frame Number')
        plt.ylabel('Flux ($\mu$Jy)')
        plt.legend(loc='best')
        plt.savefig("fig" + str(302+10*pos) + ".png")

    return

def denoise(pcf, denoisedir):
    tini = time.time()

    # Create denoising log
    logname = event.logname
    log = le.Logedit(denoisedir + "/" + logname, logname)
    log.writelog("\nStart " + denoisedir + " denoising: " + time.ctime())

    os.chdir(denoisedir)

    # copy denoise.pcf in denoisedir
    pcf.make_file("denoise.pcf")

    # Parse the attributes from the control file to the event:
    attrib = vars(pcf)
    keys = attrib.keys()
    for key in keys:
        if key != 'srcest':
            setattr(event, key, attrib.get(key).get())

    for pos in range(event.npos):
        # Plot histogram of noisy wavelet coefficients
        ylim = histwc(event, event.wavelet, event.numlvls+1, pos, log=log, denoised=False)
        # Plot first 'length' frames of noisy lightcurve at pixel srcest
        plotlc(event, pos, length=200, denoised=False)

        '''
        maxlvls  = pywt.dwt_max_level(event.nimpos[pos], pywt.Wavelet(event.wavelet))
        # Determine the number of levels to denoise
        for i in range(1,maxlvls+1):
            if (2**i)*event.framtime < event.maxtime:
                numlvls = i
            else:
                break
        '''
        log.writelog("Denoising will occur on the lowest " + str(event.numlvls) + " levels at position " +str(pos)+".")
        # Determine the time resolution of the highled denoised level
        timeres = 2**(event.numlvls)*event.framtime
        log.writelog("Time resolution for position " +str(pos)+ ", level " + str(event.numlvls) + " is " +str(timeres)+ " seconds.")

        # Assess presence of NaNs and Infs in masked data
        print("Checking for NaNs and Infs.")
        data = (event.data[:,:,:,pos])[np.where(event.mask[:,:,:,pos])]
        if (np.sum(np.isnan(data)) + np.sum(np.isinf(data))) > 0:
            log.writelog("***WARNING: Found NaNs and/or Infs in masked data at position " +str(pos)+".")

        del(data)
        pool   = mp.Pool(event.ncpu)
        for i in range(event.nx):
            for j in range(event.ny):
                exec('res = pool.apply_async(' + event.threshold + ',((event.data[:,j,i,pos])[np.where(event.mask[:,j,i,pos])], event.wavelet, event.numlvls, [j,i,pos]),callback=writedata)')

        pool.close()
        pool.join()
        res.wait()

        #Plot histogram of denoised wavelet coefficients
        histwc(event, event.wavelet, event.numlvls+1, pos, log=log, denoised=True, ylim=ylim)
        # Plot first 'length' frames of denoised lightcurve at pixel srcest
        plotlc(event, pos, length=200, denoised=True)

    # Save
    print("\nFinished Denoising. Saving.")
    me.saveevent(event, event.eventname + "_den", save=['data', 'uncd', 'mask'])

    # Print time elapsed and close log:
    cwd = os.getcwd() + "/"
    log.writelog("Output files (" + event.denoisedir + "):")
    log.writelog("Data:")
    log.writelog(" " + cwd + event.eventname + "_den.dat")
    log.writelog(" " + cwd + event.eventname + "_den.h5")
    log.writelog("Log:")
    log.writelog(" " + cwd + event.logname)

    dt = t.hms_time(time.time()-tini)
    log.writeclose("\nEnd Denoising. Time (h:m:s):  %s"%dt  +
                 "  (" + event.denoisedir + ")")
    print("-------------  ------------\n")
    return

def run_denoising(eventname, control):
    """
    Load the event.
    Read the control file.
    Launch a thread for each centering run.
    """
    global event

    pcf = rd.read_pcf(control)
    nruns = len(pcf)
    cwd = os.getcwd()

    if nruns == 1: #, I may be in the denoise dir, to re-run:
        # Get name of denoising dir:
        denoisedir = pcf[0].wavelet.get() +'_'+ pcf[0].threshold.get() +'_L'+ str(pcf[0].numlvls.get())
        if pcf[0].pcfname.get() != "":
            denoisedir += "_" + pcf[0].pcfname.get()
        if cwd[-len(denoisedir):] == denoisedir:
            # Go to dir where poet2 files were saved.
            cwd = cwd[:-len(denoisedir)]
            os.chdir(cwd)

    # Loop over each run:
    for run in np.arange(nruns):
        os.chdir(cwd)

        # Load a fresh event:
        print("Loading " + eventname)
        event = me.loadevent(eventname, load=['data','uncd','mask'])

        # Name of the directory to put the results:
        denoisedir = pcf[run].wavelet.get() +'_'+ pcf[run].threshold.get() +'_L'+ str(pcf[run].numlvls.get())
        if pcf[run].pcfname.get() != "":
            denoisedir += "_" + pcf[run].pcfname.get()
        event.denoisedir = denoisedir

        # Create the denoising directory if it doesn't exist:
        if not os.path.exists(denoisedir):
            os.mkdir(denoisedir)

        # copy the centering and photometry control files to denoisedir
        shutil.copy('center.pcf', denoisedir + '/center.pcf')
        shutil.copy('photom.pcf', denoisedir + '/photom.pcf')

        # Modify source estimate
        if hasattr(pcf[run], 'srcest'):
            if nruns == event.npos:
                event.srcest[:,run] = pcf[run].srcest.get()
            else:
                for pos in range(event.npos):
                    event.srcest[:,pos] = pcf[0].srcest.get()

        #Call denoising on each wavelet:
        denoise(pcf[run], denoisedir)

    return
