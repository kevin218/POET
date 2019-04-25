# $Author: patricio $
# $Revision: 285 $
# $Date: 2010-06-18 17:59:25 -0400 (Fri, 18 Jun 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_6model.py $
# $Id: poet_6model.py 285 2010-06-18 21:59:25Z patricio $


"""
 MODIFICATION HISTORY:
    Written by:	Kevin Stevenson, UCF  	2008-07-02
                kevin218@knights.ucf.edu
    Finished initial version:       kevin    2008-09-08
    Updated for multi events:       kevin    2009-11-01
    Added ip interpolation:         kevin    2010-06-28
    Added minnumpts & fixipmap:     kevin    2010-08-03
    Least-squares w/ multi events:  kevin    2010-08-19
    Modified to a chisq minimizer   patricio 2010-10-22 pcubillos@fulbrightmail.org
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import numexpr as ne
import scipy.optimize as op
import time, os, sys
import models2 as models
reload(models)
import mcmc
reload(mcmc)
import paramedit as pe
reload(pe)
import smoothing
reload(smoothing)
import readeventhdf

"""
 NAME:
	RUNDMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	fit = RUNDMC(fit, event, numit, nbins, chi2flag=0, num=0)

 INPUTS:
	event:    List of event objects
	num:	  Iteration number for runs with multiple models
    printout: File object for directing print statements

 OUTPUTS:
	None

 PROCEDURE:

 EXAMPLE:

 MODIFICATION HISTORY:
    Written by:	Kevin Stevenson, UCF  	2008-07-02
                kevin218@knights.ucf.edu
    Finished initial version:   kevin   2008-09-08
    Updated for multi events:   kevin   2009-11-01
    Added normflux option:      kevin   2010-01-20
    Least-squares fitting:      kevin   2010-03-**
    Added ip interpolation:     kevin   2010-06-28
"""

def read_parameters(file_name):
    # The possible models
    models = np.array(["mandelecl",   "mandelecl2", "mandelgeom",    
                       "mandelorbit", "risingexp",  "expramp", 
                       "re2ramp",     "reqramp",    "relramp",                
                       "fallingexp",  "felramp",    "quadramp",     
                       "linramp",     "logramp",    "log4qramp",     
                       "llramp",      "lqramp",     "sindecay",    
                       "sincos",      "quadip",     "quadip4",   
                       "cubicip",     "sexticip",   "sexticipc", 
                       "medianip",    "nnint",      "bilinint",      
                       "ipspline",    "posflux",    "posfluxlinip",      
                       "vsll",        "vsspline",   "flatfield3",  
                       "not0risingexp", "batman_ecl", "batman_trquad"])

    # Read file
    f = open(file_name, "r")
    lines = f.readlines()
    f.close()
    sline = []
    for line in lines:
        sline.append(line.strip())

    # Output
    mymodels = []
    for i in np.arange(len(sline)):
        # Find a model
        if sline[i] in models:
            npars = len(sline[i+1].split('\t'))
            j=3
            # Count lines until next model or until end
            while(sline[i+j] != '****END FILE****' and 
                  (not sline[i+j] in models)):
                j += 1
            # Fill info
            info = np.zeros((j-2, npars))
            info[0,:] = sline[i+2].split('\t')
            for k in np.arange(1,j-2):
                info[k,:] = sline[i+2+k].split('\t')
            # Update output
            mymodels.append([sline[i], info])

    print('Parameters have been successfully read from ' + file_name)
    return mymodels

def formatpar(arr):
    s = ''
    s += '%7.4e'%arr[0]
    for i in np.arange(1, len(arr)):
        s += '\t%7.4e'%arr[i]

    return s



def minimizer(event, num=0, printout=sys.stdout):
    # When allplots=False: the correlation and histogram figures are 
    #                      not produced.
    # When normflux=Flase: flux is not normalized w.r.t. each position's 
    #                      median value.
    #                      This should only be used when using normflux model.

    chi2flag     = event[0].params.chi2flag
    allplots     = event[0].params.allplots
    leastsq      = event[0].params.leastsq
    nbins        = []
    normflux     = []
    obj          = event[0].filename
    numevents    = len(event)
    nummodels    = np.zeros(numevents, dtype=int)
    fit          = []
    params       = []
    pmin         = []
    pmax         = []
    numleasts    = 0
    fparams      = []
    numparams    = [0]
    iparams      = []
    myfuncs      = []
    myrampfunc   = []
    functype     = []
    parlist      = []
    numfigs      = 10
    pcounter     = 0
    rampidx      = []

    
    for j in range(numevents):
        fit.append(readeventhdf.fits())
        fit[j].model   = event[j].params.model[num]
        fit[j].modelslist = []
        event[j].params.chi2flag = chi2flag
        nbins   .append(event[j].params.nbins)
        normflux.append(event[j].params.normflux)


        # Read in data
        fit[j].mflux = np.mean(event[j].aplev[np.where(event[j].good == 1)])
        # normalize at each position if npos > 1
        if ((event[j].npos > 1) and (normflux == True)):
            npos            = event[j].good.shape[0]
            fit[j].posmflux = np.zeros(npos)
            for i in range(npos):
                posgood     = np.where(event[j].good[i] == 1)
                fit[j].posmflux[i] = np.mean(event[j].aplev[i,posgood])
            fit[j].fluxuc   = ((event[j].aplev[:,:].T/fit[j].posmflux).T)[np.where(event[j].good == 1)] * \
                                  fit[j].mflux
            fit[j].sigmauc  = ((event[j].aperr[:,:].T/fit[j].posmflux).T)[np.where(event[j].good == 1)] * \
                                  fit[j].mflux
        else:
            fit[j].fluxuc   = event[j].aplev[np.where(event[j].good == 1)]
            fit[j].sigmauc  = event[j].aperr[np.where(event[j].good == 1)]

        fit[j].phaseuc  = event[j].phase [np.where(event[j].good == 1)]
        fit[j].yuc      = event[j].y     [np.where(event[j].good == 1)]
        fit[j].xuc      = event[j].x     [np.where(event[j].good == 1)]
        fit[j].time     = event[j].time  [np.where(event[j].good == 1)]
        fit[j].bjdcor   = event[j].bjdcor[np.where(event[j].good == 1)]
        fit[j].juldat   = event[j].juldat[np.where(event[j].good == 1)]
        fit[j].posuc    = event[j].pos   [np.where(event[j].good == 1)]
        fit[j].frmvisuc = event[j].frmvis[np.where(event[j].good == 1)]
        if hasattr(event[j], 'apdata'):
            fit[j].apdatauc = event[j].apdata[np.where(event[j].good == 1)]

        # Sort data by phase
        if event[j].npos > 1:
            isort           = np.argsort(fit[j].phaseuc)
            fit[j].phaseuc  = fit[j].phaseuc[isort]
            fit[j].fluxuc   = fit[j].fluxuc[isort]
            fit[j].sigmauc  = fit[j].sigmauc[isort]
            fit[j].yuc      = fit[j].yuc[isort]
            fit[j].xuc      = fit[j].xuc[isort]
            fit[j].time     = fit[j].time[isort]
            fit[j].posuc    = fit[j].posuc[isort]
            fit[j].frmvisuc = fit[j].frmvisuc[isort]
            if hasattr(event[j], 'apdata'):
                fit[j].apdatauc = fit[j].apdatauc[isort]


        # Obtain initial model parameters
        # Read from file
        if event[j].params.modelfile.__class__ == str:  # IF   params.modelfile is a string
            modelfile = event[j].params.modelfile
        else:                                           # ELSE params.modelfile is a list
            modelfile = event[j].params.modelfile[num]

        if os.path.isfile(modelfile) == False:
            print(modelfile + " does not exist!")
            return

        # Initialize models
        parlist.append(read_parameters(modelfile))
        nummodels[j] = len(parlist[j])
        for i in np.arange(nummodels[j]):
            fit[j].modelslist.append(parlist[j][i][0])

        foo, bar, fit[j].parname, fit[j].i, fit[j].saveext = models.setupmodel(fit[j].modelslist, fit[j].i)
        myfuncs  = np.concatenate(( myfuncs, foo), 0)
        functype = np.concatenate((functype, bar), 0)
        myrampfuncs = myfuncs[np.where(functype == 1)]
        if np.where(bar == 6) != (nummodels[j] - 1):
            print("ERROR: The interpolation model must be the last model listed.")
            return

        # Count number of models/parameters/freeparameters
        npars      = 0  # total number of parameters per event
        indrpars   = []
        ramprange  = [0]
        for i in np.arange(nummodels[j]):
            if functype[i] == 1:
                numleasts += len((parlist[j][i][1])[1:,0])
                indrpars.append([npars,npars+len((parlist[j][i][1])[0])])
                ramprange.append(numleasts)
            npars += len((parlist[j][i][1])[0])
        fparams = np.zeros((numleasts, npars))

        iramp = 0 # ramp index
        for i in range(nummodels[j]):
            pars       = parlist[j][i][1]
            if fit[j].modelslist[i] == 'ipspline':  # This IF is going to fail!
                # read number of knots along x and y
                numipyk, numipxk = pars[0]
                # temporarily populate initial intra-pixel parameters with ones
                temppars = [np.ones(numipxk*numipyk)*1.0,
                            np.ones(numipxk*numipyk)*pars[1][0],
                            np.ones(numipxk*numipyk)*pars[2][0],
                            np.ones(numipxk*numipyk)*pars[3][0]]
                pars     = temppars

            params     = np.concatenate((params,   pars[0]),0)
            if functype[i] != 1:  # non ramp models
                fparams[:,pcounter:pcounter+len(pars[0])] = pars[1]
            else:                 # ramp models
                fparams[ramprange[iramp]:ramprange[iramp+1],pcounter:pcounter+len(pars[0])] = pars[1:]
                iramp = iramp + 1
            numparams  = np.concatenate((numparams, [numparams[-1] + len(pars[0])]),0)
            iparams.append(range(pcounter,pcounter+len(pars[0])))

            # CHECK FOR SHARED PARAMETERS, MAKE APPROPRIATE CHANGES TO PARAMETER VALUES AND LIMITS
            # FINDME: I guess i'm not supporting shared parameters yet
#             for k in range(len(pars[0])):
#                 if pars[3][k] < 0:
#                     params [pcounter+k] = params [-pars[3][k]-1]
#                     pmin   [pcounter+k] = pmin   [-pars[3][k]-1]
#                     pmax   [pcounter+k] = pmax   [-pars[3][k]-1]
#                     iparams[-1][k]      = int(-pars[3][k]-1)

            pcounter += len(pars[0])
        
        q = np.where(functype == 1)[0]
        rampidx = np.zeros(numleasts)
        for r in np.arange(len(q)):
            rampidx[ramprange[r]:ramprange[r+1]] = q[r]
        


    # Get indices of the parameters
    cummodels    = [0]      # Cumulative list of number of models
    funcx        = []
    funcxuc      = []
    myfreepars   = []
    ifreepars    = np.array([], dtype=int)
    nonfixedpars = np.array([], dtype=int)
    text         = ""
    modelnames   = []
    maxnumfp     = 0
    numfreepars  = np.sum(fparams, axis=1) # np.array(np.where(stepsize > 0)).flatten().size
    print("\nCurrent event & model:", file=printout)

    for j in range(numevents):
        print(event[j].eventname, file=printout)
        print(fit[j].modelslist, file=printout)
        cummodels.append(int(cummodels[j] + nummodels[j]))


        # Specify indices of free parameters
        for l in np.arange(numleasts):
            myfreepars.append(np.array(np.where(fparams[l,numparams[cummodels[j]]:numparams[cummodels[j+1]]]!=0)).flatten())
        fit[j].numfreepars = numfreepars # fit[j].myfreepars.size

        # Determine the centroid location relative to the center of the pixel
        # Record in which quadrant the image center falls
        fit[j].nobjuc   = fit[j].fluxuc.size
        fit[j].quadrant = np.zeros(fit[j].nobjuc)
        fit[j].numq     = np.zeros(4)
        if fit[j].modelslist.__contains__('quadip4'):
            fit[j].y        = fit[j].yuc - np.round(np.median(fit[j].yuc))
            fit[j].x        = fit[j].xuc - np.round(np.median(fit[j].xuc))
        else:
            fit[j].y        = fit[j].yuc - np.round(np.median(fit[j].yuc))
            fit[j].x        = fit[j].xuc - np.round(np.median(fit[j].xuc))
        for i in range(fit[j].nobjuc):
            if fit[j].y[i] > event[j].params.ydiv:
                if fit[j].x[i] < event[j].params.xdiv:
                    fit[j].quadrant[i] = 0   #Upper left
                    fit[j].numq[0]    += 1
                else:
                    fit[j].quadrant[i] = 1   #Upper right
                    fit[j].numq[1]    += 1
            else:
                if fit[j].x[i] < event[j].params.xdiv:
                    fit[j].quadrant[i] = 3   #Lower left
                    fit[j].numq[3]    += 1
                else:
                    fit[j].quadrant[i] = 2   #Lower right
                    fit[j].numq[2]    += 1
        print("Number of frames per pixel quadrant:", file=printout)
        print(fit[j].numq, file=printout)

        # Clip first and last few points of data
        if len(event[j].params.preclip) == len(fit[j].modelslist):  # if preclip for each model-set
            preclip     = event[j].params.preclip[num]
            postclip    = fit[j].nobjuc - event[j].params.postclip[num]
        else:
            preclip     = event[j].params.preclip[0]
            postclip    = fit[j].nobjuc - event[j].params.postclip[0]
        # Define clipmask: 1=kept, 0=clipped
        fit[j].clipmask = np.zeros(fit[j].nobjuc)
        fit[j].clipmask[preclip:postclip] = 1
        # Use interclip (IF IT EXISTS) to clip points in the middle
        if hasattr(event[j].params, 'interclip'):
            interclip   = event[j].params.interclip
            for i in range(len(interclip)):
                fit[j].clipmask[interclip[i][0]:interclip[i][1]] = 0

        # Define ipmask: 1=intrapixel interpolation model, 0=not
        fit[j].ipmaskuc = np.ones(fit[j].nobjuc)
        if hasattr(event[j].params, 'ipclip'):
            ipclip = event[j].params.ipclip
            for i in range(len(ipclip)):
                fit[j].ipmaskuc[ipclip[i][0]:ipclip[i][1]] = 0
        
        # Define temporary parameters
        fit[j].ipmask   = np.copy(fit[j].ipmaskuc[np.where(fit[j].clipmask)])
        fit[j].position = np.array([fit[j].y[np.where(fit[j].clipmask)], 
                                    fit[j].x[np.where(fit[j].clipmask)], 
                                    fit[j].quadrant[np.where(fit[j].clipmask)]])
        fit[j].nobj     = fit[j].position[0].size


        # Calculate minnumptsmask = 1 FOR ATLEAST minnumpts IN EACH BIN
        # FINDME: I haaven't look this block yet
        fit[j].minnumptsmask = np.ones(fit[j].nobj, dtype=int)
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 6:    #IP MAPPING
                # Read bin sizes
                if len(event[j].params.ystep) == len(event[j].params.modelslist):
                    ystep = event[j].params.ystep[num]
                    xstep = event[j].params.xstep[num] 
                else:
                    ystep = event[j].params.ystep[0]
                    xstep = event[j].params.xstep[0]
                yfactor  = 10.
                xfactor  = 10.
                ymin     = np.floor(yfactor*fit[j].y.min())/yfactor - ystep
                ymax     = np.ceil (yfactor*fit[j].y.max())/yfactor + ystep
                xmin     = np.floor(xfactor*fit[j].x.min())/xfactor - xstep
                xmax     = np.ceil (xfactor*fit[j].x.max())/xfactor + xstep
                # Number of bins in y,x dimensions
                ysize    = int((ymax-ymin)/ystep + 1)
                xsize    = int((xmax-xmin)/xstep + 1)
                ygrid, ystep  = np.linspace(ymin, ymax, ysize, retstep=True)
                xgrid, xstep  = np.linspace(xmin, xmax, xsize, retstep=True)
                # Minimum number of acceptable points in a bin
                if event[j].params.minnumpts.__class__ == int:
                    minnumpts = event[j].params.minnumpts
                elif len(event[j].params.minnumpts) == len(event[j].params.modelslist):
                    minnumpts = event[j].params.minnumpts[num]
                else:
                    minnumpts = event[j].params.minnumpts[0]
                # CALCULATE BINS FOR 2D BINNING
                for m in range(ysize):
                    wbftemp   = np.where(np.abs(fit[j].position[0]-ygrid[m]) < (ystep/2.))[0]
                    for n in range(xsize):
                        wbf       = wbftemp[np.where((np.abs(fit[j].position[1,[wbftemp]]-xgrid[n]) < (xstep/2.))[0])]
                        wbfipmask = wbf    [np.where(fit[j].ipmask[wbf] == 1)]
                        if len(wbfipmask) < minnumpts:
                            fit[j].minnumptsmask[wbf] = 0

        #REDEFINE CLIPPED VARIABLES BASED ON minnumpts FOR IP MAPPING
        fit[j].clipmask[np.where(fit[j].clipmask)] *= fit[j].minnumptsmask
        fit[j].phase    = fit[j].phaseuc[np.where(fit[j].clipmask)]
        fit[j].flux     = np.copy(fit[j].fluxuc  [np.where(fit[j].clipmask)])
        fit[j].sigma    = np.copy(fit[j].sigmauc [np.where(fit[j].clipmask)])
        fit[j].pos      = np.copy(fit[j].posuc   [np.where(fit[j].clipmask)])
        fit[j].frmvis   = np.copy(fit[j].frmvisuc[np.where(fit[j].clipmask)])
        fit[j].ipmask   = np.copy(fit[j].ipmaskuc[np.where(fit[j].clipmask)])
        if hasattr(event[0], 'apdata'):
            fit[j].apdata   = np.copy(fit[j].apdatauc[np.where(fit[j].clipmask)])
        fit[j].position = np.array([fit[j].y[np.where(fit[j].clipmask)], 
                                    fit[j].x[np.where(fit[j].clipmask)], 
                                    fit[j].quadrant[np.where(fit[j].clipmask)]])
        fit[j].nobj     = fit[j].flux.size
        fit[j].positionuc = np.array([fit[j].y, fit[j].x, fit[j].quadrant])


        #DETERMINE BIN LOCATION FOR EACH POSITION
        fit[j].isipmapping = False
        fit[j].numknots = 0
        for i in range(cummodels[j],cummodels[j+1]):
            # FINDME: I'm skipping this also!
            if functype[i] == 6:    #nnint, bilinint
                fit[j].isipmapping    = True    #Using mapping for IP sensitivity
                fit[j].wherebinflux   = []      #List of size = # of bins, def which points fall into each bin
                fit[j].wherebinfluxuc = []      #Un-clipped version of above
                fit[j].wbfipmask      = []      #Intrapixel masked version of wherebinflux
                fit[j].wbfipmaskuc    = []      #Un-clipped version of above
                fit[j].binloc         = np.zeros((2, fit[j].nobj),   dtype=int) - 1
                fit[j].binlocuc       = np.zeros((2, fit[j].nobjuc), dtype=int) - 1
                #Read bin sizes
                if len(event[j].params.ystep) == len(event[j].params.modelslist):
                    ystep = event[j].params.ystep[num]
                    xstep = event[j].params.xstep[num] 
                else:
                    ystep = event[j].params.ystep[0]
                    xstep = event[j].params.xstep[0]  
                #yfactor  = 0.2/ystep
                #xfactor  = 0.2/xstep
                yfactor  = 10.
                xfactor  = 10.
                ymin     = np.floor(yfactor*fit[j].y.min())/yfactor - ystep
                ymax     = np.ceil (yfactor*fit[j].y.max())/yfactor + ystep
                xmin     = np.floor(xfactor*fit[j].x.min())/xfactor - xstep
                xmax     = np.ceil (xfactor*fit[j].x.max())/xfactor + xstep
                #Number of bins in y,x dimensions
                ysize    = int((ymax-ymin)/ystep + 1)
                xsize    = int((xmax-xmin)/xstep + 1)
                ygrid, ystep  = np.linspace(ymin, ymax, ysize, retstep=True)
                xgrid, xstep  = np.linspace(xmin, xmax, xsize, retstep=True)
                fit[j].xygrid = np.meshgrid(xgrid, ygrid)
                fit[j].binfluxmask   = np.zeros((ysize, xsize), dtype=int)
                fit[j].binfluxmaskuc = np.zeros((ysize, xsize), dtype=int)
                fit[j].numpts        = np.zeros((ysize, xsize)) #Number of points per bin
                fit[j].numptsuc      = np.zeros((ysize, xsize))
                #Minimum number of acceptable points in a bin
                if event[j].params.minnumpts.__class__ == int:
                    minnumpts = event[j].params.minnumpts
                elif len(event[j].params.minnumpts) == len(event[j].params.modelslist):
                    minnumpts = event[j].params.minnumpts[num]
                else:
                    minnumpts = event[j].params.minnumpts[0]
                print('Step size in y = ' + str(ystep), file=printout)
                print('Step size in x = ' + str(xstep), file=printout)
                print('Ignoring bins with < ' + str(minnumpts) + ' points.', file=printout)
                print('Computing bin for each position.')
                #ASSIGN BINS FOR 2D BINNING
                for m in range(ysize):
                    wbftemp   = np.where(np.abs(fit[j].position[0]-ygrid[m]) < (ystep/2.))[0]
                    wbftempuc = np.where(np.abs(fit[j].y          -ygrid[m]) < (ystep/2.))[0]
                    for n in range(xsize):
                        wbf   = wbftemp[np.where((np.abs(fit[j].position[1,[wbftemp]]-xgrid[n]) < (xstep/2.))[0])]
                        wbfuc = wbftempuc[np.where(np.abs(fit[j].x[wbftempuc]-xgrid[n]) < (xstep/2.))]
                        wbfipmask   = wbf  [np.where(fit[j].ipmask  [wbf  ] == 1)]
                        wbfipmaskuc = wbfuc[np.where(fit[j].ipmaskuc[wbfuc] == 1)]
                        if len(wbfipmask) >= minnumpts:
                            fit[j].binfluxmask  [m,n] = 1
                            fit[j].numpts       [m,n] = len(wbfipmask)
                            fit[j].binloc  [0, wbf  ] = m*xsize + n
                            fit[j].wherebinflux.append(wbf)
                            fit[j].wbfipmask.  append(wbfipmask)
                        else:
                            fit[j].wherebinflux.append([])
                            fit[j].wbfipmask.  append([])
                        if len(wbfipmaskuc) > 0:
                            fit[j].binfluxmaskuc[m,n] = 1
                            fit[j].numptsuc     [m,n] = len(wbfipmaskuc)
                            fit[j].binlocuc[0, wbfuc] = m*xsize + n
                            fit[j].wherebinfluxuc.append(wbfuc)
                            fit[j].wbfipmaskuc.append(wbfipmaskuc)
                        else:
                            fit[j].wherebinfluxuc.append([])
                            fit[j].wbfipmaskuc.append([])
                
                fit[j].numknots = np.sum(fit[j].binfluxmask)
                #Read smoothing paramters
                if len(event[j].params.nx) == len(event[j].params.modelslist):
                    ny = event[j].params.ny[num]
                    nx = event[j].params.nx[num] 
                else:
                    ny = event[j].params.ny[0]
                    nx = event[j].params.nx[0] 
                if len(event[j].params.sx) == len(event[j].params.modelslist):
                    sy = event[j].params.sy[num]
                    sx = event[j].params.sx[num] 
                else:
                    sy = event[j].params.sy[0]
                    sx = event[j].params.sx[0]
                #Adjust weighting based on number of points in bin
                #weightbinfluxmask   = np.sqrt(fit[j].numpts)   * fit[j].binfluxmask
                #weightbinfluxmaskuc = np.sqrt(fit[j].numptsuc) * fit[j].binfluxmaskuc
                weightbinfluxmask   = fit[j].binfluxmask
                weightbinfluxmaskuc = fit[j].binfluxmaskuc
                #Calculate smoothing kernel
                fit[j].smoothingp = [ny,nx,sy,sx]
                fit[j].kernel     = smoothing.gauss_kernel_mask((ny,nx), (sy,sx), weightbinfluxmask  )
                fit[j].kerneluc   = smoothing.gauss_kernel_mask((ny,nx), (sy,sx), weightbinfluxmaskuc)

                fit[j].binfluxmask   = fit[j].binfluxmask  .flatten()
                fit[j].binfluxmaskuc = fit[j].binfluxmaskuc.flatten()

                #DETERMINE DISTANCES TO FOUR NEAREST GRID POINTS
                #USED FOR BILINEAR INTERPOLATION
                #ORDER [grid #, (y-y1)/ystep, (y2-y)/ystep, (x-x1)/xstep, (x2-x)/xstep)
                print('Computing distances to four nearest bins.')
                fit[j].griddist       = np.ones((4, fit[j].nobj  ))
                fit[j].griddistuc     = np.ones((4, fit[j].nobjuc))
                for m in range(ysize-1):
                    wherey = np.where(np.bitwise_and(fit[j].position[0] > ygrid[m  ],
                                                     fit[j].position[0] < ygrid[m+1]))[0]
                    for n in range(xsize-1):
                        wherexy = wherey[np.where(np.bitwise_and(fit[j].position[1,[wherey]] > xgrid[n  ],
                                                                 fit[j].position[1,[wherey]] < xgrid[n+1])[0])[0]]
                        if len(wherexy) > 0:
                            fit[j].binloc[1, wherexy] = gridpt = m*xsize + n
                            #IF THERE ARE NO POINTS IN ONE OR MORE BINS...
                            if (len(fit[j].wbfipmask[gridpt        ]) == 0) or \
                               (len(fit[j].wbfipmask[gridpt      +1]) == 0) or \
                               (len(fit[j].wbfipmask[gridpt+xsize  ]) == 0) or \
                               (len(fit[j].wbfipmask[gridpt+xsize+1]) == 0):
                                #SET griddist = NEAREST BIN (USE NEAREST NEIGHBOR INTERPOLATION)
                                for loc in wherexy:
                                    if   loc in fit[j].wherebinflux[gridpt        ]:
                                        fit[j].griddist[0, loc] = 0
                                        fit[j].griddist[2, loc] = 0
                                    elif loc in fit[j].wherebinflux[gridpt      +1]:
                                        fit[j].griddist[0, loc] = 0
                                        fit[j].griddist[3, loc] = 0
                                    elif loc in fit[j].wherebinflux[gridpt+xsize  ]:
                                        fit[j].griddist[1, loc] = 0
                                        fit[j].griddist[2, loc] = 0
                                    elif loc in fit[j].wherebinflux[gridpt+xsize+1]:
                                        fit[j].griddist[1, loc] = 0
                                        fit[j].griddist[3, loc] = 0
                            else:
                                #CALCULATE griddist NORMALLY FOR BILINEAR INTERPOLATION
                                fit[j].griddist[0, wherexy] = np.array((fit[j].position[0][wherexy]-ygrid[m])  /ystep)
                                fit[j].griddist[1, wherexy] = np.array((ygrid[m+1]-fit[j].position[0][wherexy])/ystep)
                                fit[j].griddist[2, wherexy] = np.array((fit[j].position[1][wherexy]-xgrid[n])  /xstep)
                                fit[j].griddist[3, wherexy] = np.array((xgrid[n+1]-fit[j].position[1][wherexy])/xstep)
                #REPEAT FOR uc
                #print("50% complete...")
                for m in range(ysize-1):
                    wherey = np.where(np.bitwise_and(fit[j].y > ygrid[m  ],
                                                     fit[j].y < ygrid[m+1]))[0]
                    for n in range(xsize-1):
                        wherexy = wherey[np.where(np.bitwise_and(fit[j].x[wherey] > xgrid[n  ],
                                                                 fit[j].x[wherey] < xgrid[n+1]))[0]]
                        if len(wherexy) > 0:
                            fit[j].binlocuc[1, wherexy] = gridpt = m*xsize + n
                            #IF THERE ARE NO POINTS IN ONE OR MORE BINS...
                            if (len(fit[j].wbfipmaskuc[gridpt        ]) == 0) or \
                               (len(fit[j].wbfipmaskuc[gridpt      +1]) == 0) or \
                               (len(fit[j].wbfipmaskuc[gridpt+xsize  ]) == 0) or \
                               (len(fit[j].wbfipmaskuc[gridpt+xsize+1]) == 0):
                                #SET griddist = NEAREST BIN (USE NEAREST NEIGHBOR INTERPOLATION)
                                for loc in wherexy:
                                    if loc in fit[j].wherebinfluxuc[gridpt        ]:
                                        fit[j].griddistuc[0, loc] = 0
                                        fit[j].griddistuc[2, loc] = 0
                                    if loc in fit[j].wherebinfluxuc[gridpt      +1]:
                                        fit[j].griddistuc[0, loc] = 0
                                        fit[j].griddistuc[3, loc] = 0
                                    if loc in fit[j].wherebinfluxuc[gridpt+xsize  ]:
                                        fit[j].griddistuc[1, loc] = 0
                                        fit[j].griddistuc[2, loc] = 0
                                    if loc in fit[j].wherebinfluxuc[gridpt+xsize+1]:
                                        fit[j].griddistuc[1, loc] = 0
                                        fit[j].griddistuc[3, loc] = 0
                            else:
                                #CALCULATE griddist NORMALLY FOR BILINEAR INTERPOLATION
                                fit[j].griddistuc[0, wherexy] = np.array((fit[j].y[wherexy]-ygrid[m])  /ystep)
                                fit[j].griddistuc[1, wherexy] = np.array((ygrid[m+1]-fit[j].y[wherexy])/ystep)
                                fit[j].griddistuc[2, wherexy] = np.array((fit[j].x[wherexy]-xgrid[n])  /xstep)
                                fit[j].griddistuc[3, wherexy] = np.array((xgrid[n+1]-fit[j].x[wherexy])/xstep)
                                
                fit[j].posflux  = [fit[j].y[np.where(fit[j].clipmask)], 
                                   fit[j].x[np.where(fit[j].clipmask)], 
                                   fit[j].flux,
                                   fit[j].wbfipmask,
                                   fit[j].binfluxmask,
                                   fit[j].kernel,
                                   fit[j].smoothingp,
                                   fit[j].binloc,
                                   fit[j].griddist,
                                   fit[j].xygrid[0].shape,
                                   event[j].params.issmoothing]
                fit[j].posfluxuc= [fit[j].y,
                                   fit[j].x,
                                   fit[j].fluxuc,
                                   fit[j].wbfipmaskuc,
                                   fit[j].binfluxmaskuc,
                                   fit[j].kerneluc,
                                   fit[j].smoothingp,
                                   fit[j].binlocuc,
                                   fit[j].griddistuc,
                                   fit[j].xygrid[0].shape,
                                   event[j].params.issmoothing]


    # RESTART LOOP: Define funcx, fit[j].etc
    # Assigns phase or x,y-position as independent variable for model fitting
    for j in range(numevents):
        k = 0
        fit[j].etc = []
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 0:
                funcx.  append(fit[j].phase)
                funcxuc.append(fit[j].phaseuc)
                fit[j].etc.append([])
            elif functype[i] == 1:
                funcx.  append(fit[j].phase)
                funcxuc.append(fit[j].phaseuc)
                fit[j].etc.append([])
            elif functype[i] == 2:
                if fit[j].modelslist[k] == 'ipspline':
                    funcx.  append(fit[j].position)
                    funcxuc.append(fit[j].positionuc)
                    #CREATE KNOTS FOR IPSPLINE
                    fit[j].etc.append(np.meshgrid(
                                  np.linspace(fit[j].position[1].min(),fit[j].position[1].max(),numipxk), 
                                  np.linspace(fit[j].position[0].min(),fit[j].position[0].max(),numipyk)))
                if fit[j].modelslist[k] == 'posfluxlinip':
                    fit[j].wherepos   = []
                    fit[j].whereposuc = []
                    fit[j].meanypos   = []
                    fit[j].meanxpos   = []
                    for i in range(event[j].npos):
                        wherepos = np.where(fit[j].pos   == i)[0]
                        fit[j].wherepos  .append(wherepos)
                        fit[j].meanypos  .append(np.mean(fit[j].position[0][wherepos]))
                        fit[j].meanxpos  .append(np.mean(fit[j].position[1][wherepos]))
                        fit[j].whereposuc.append(np.where(fit[j].posuc == i)[0])
                    funcx.  append([fit[j].position,   fit[j].nobj,   fit[j].wherepos])
                    funcxuc.append([fit[j].positionuc, fit[j].nobjuc, fit[j].whereposuc])
                    fit[j].etc.append([fit[j].meanypos, fit[j].meanxpos])
                else:
                    funcx.  append(fit[j].position)
                    funcxuc.append(fit[j].positionuc)
                    fit[j].etc.append([])
            elif functype[i] == 3:
                fit[j].wherepos   = []
                fit[j].whereposuc = []
                for i in range(event[j].npos):
                    fit[j].wherepos  .append(np.where(fit[j].pos   == i)[0])
                    fit[j].whereposuc.append(np.where(fit[j].posuc == i)[0])
                funcx.  append([fit[j].nobj,   fit[j].wherepos])
                funcxuc.append([fit[j].nobjuc, fit[j].whereposuc])
                fit[j].etc.append([])
            elif functype[i] == 4:
                funcx.  append((fit[j].frmvis,   event[j].params.vsknots))
                funcxuc.append((fit[j].frmvisuc, event[j].params.vsknots))
                fit[j].etc.append([])
            elif functype[i] == 5:
                funcx.  append(fit[j].apdata)
                funcxuc.append(fit[j].apdatauc)
                fit[j].etc.append([])
            elif functype[i] == 6:
                funcx.  append(fit[j].posflux)
                #funcxuc.append([fit[j].y, fit[j].x, fit[j].fluxuc, fit[j].whereltraduc])   #medianip
                funcxuc.append(fit[j].posfluxuc)
                fit[j].etc.append([])
            text   += fit[j].modelslist[k] + " "
            modelnames = np.append(modelnames, fit[j].modelslist[k])
            k      += 1

    domodel = np.zeros((numleasts, len(functype)))
    for i in np.arange(numleasts):
        domodel[i,:] = functype != 1
        domodel[i,rampidx[i]] = 1

    # **** LEAST-SQUARES FIT *****
    flux      = []
    phase     = []
    sigma     = []
    for j in range(numevents):
        flux. append(fit[j].flux)
        phase.append(fit[j].phase)
        sigma.append(fit[j].sigma)
    
    print('\nParams before minimizer:')
    print(params)
    print('\n')
    initialpars = np.copy(params)

    bic       = np.zeros(numleasts)
    bic2      = np.zeros(numleasts)
    sdnr      = np.zeros(numleasts)
    redchisq  = np.zeros(numleasts)
    residuals = []
    bestfit   = []

    fout = [] # Output file(s)
    for j in np.arange(numevents):
        if event[j].params.modelfile.__class__ == str:
            fout.append(open("output_" + event[j].params.modelfile, "w"))
            print("output_" + event[j].params.modelfile)
        else:
            fout.append(open("output_" + event[j].params.modelfile[num], "w"))
            print( "output_" + event[j].params.modelfile[num])

    # Loop for each ramp model:
    for m in np.arange(numleasts):
        params = np.copy(initialpars)
        # ================== A DEF inside the code?? ============================
        def callmodelfunc(freepars, params, sigma):
            params[myfreepars[m]] = freepars
            residuals = []
            for j in range(numevents):
                fit0 = np.ones(fit[j].nobj)
                k    = 0
                for i in range(cummodels[j],cummodels[j+1]):
                    if domodel[m,i]:
                        if functype[i] == 6:
                            fit[j].etc[k] = fit0
                        fit0 *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
                        k    += 1
                residuals = np.concatenate((residuals,(fit0 - flux[j])/sigma[j]),axis=0)
            return residuals
        # =================== A DEF inside the code?? ===========================

        print("Calculating least-squares fit.")
        modelfunc = lambda freepars, params, sigma: callmodelfunc(freepars, params, sigma)
        output, err = op.leastsq(modelfunc, params[myfreepars[m]], args=(params, sigma), 
                                 factor=100, ftol=1e-16, xtol=1e-16, gtol=1e-16)
        if (err >= 1) and (err <= 4):
            print("Fit converged without error.")
        else:
            print("WARNING: Fit did not converge. Using last iteration.")

        print('\nMinimizer results:')
        print(output)
        #return

#         # update shared parameters
#         for i in range(len(stepsize)):
#             if stepsize[i] < 0:
#                 params[i] = params[-stepsize[i]-1]
    # ************  **************

        # Create model using initial fit values
        saveext = ''
        if 6 in functype:
            plt.figure(6000+m)
            plt.clf()
        for j in range(numevents):
            saveext     = saveext + fit[j].saveext
            fit[j].fit0 = np.ones(fit[j].nobj)
            k           = 0
            for i in range(cummodels[j],cummodels[j+1]):
                if functype[i] == 6:
                    ipflux, fit[j].binipflux, fit[j].binipstd = myfuncs[i](params[iparams[i]], \
                                                         funcx[i], fit[j].fit0, retbinflux=True, retbinstd=True)
                    fit[j].fit0     *= ipflux
                    fit[j].binipflux = fit[j].binipflux.reshape(fit[j].xygrid[0].shape)
                    fit[j].binipstd  = fit[j].binipstd. reshape(fit[j].xygrid[0].shape)
                    # Define IP confidence level (or statistical significance) per bin
                    # Confidence = signal / (noise / sqrt(sample size))
                    # fit[j].ipconf = np.sqrt(fit[j].numpts)*
                    #                 (fit[j].binipflux-fit[j].binipflux.mean())/fit[j].binipstd
                    fit[j].ipconf    = np.sqrt(fit[j].numpts)*(fit[j].binipflux)/fit[j].binipstd
                    fit[j].ipconf[np.where(np.isnan(fit[j].ipconf))] = 0
                    fit[j].ipconf[np.where(np.isinf(fit[j].ipconf))] = 0
    
                    numsteps = len(event[j].params.sssteps)-1
                    ipmean   = np.zeros(numsteps)
                    ipstd    = np.zeros(numsteps)
                    ipmin    = np.zeros(numsteps)
                    ipmax    = np.zeros(numsteps)
                    for k in range(numsteps):
                        instep = np.where(np.bitwise_and(fit[j].numpts >  event[j].params.sssteps[k  ], \
                                                         fit[j].numpts <= event[j].params.sssteps[k+1]))
                        if len(instep[0]) > 0:
                            ipmean[k] = np.mean(fit[j].ipconf[instep])
                            ipstd [k] = np.std(fit[j].ipconf[instep])
                            ipmin [k] = np.min(fit[j].ipconf[instep])
                            ipmax [k] = np.max(fit[j].ipconf[instep])
                    
                    # PLOT CONFIDENCE vs. NUMBER OF POINTS PER BIN
                    # ERROR BARS ARE MINIMUM AND MAXIMUM, NOT STANDARD DEVIATION
                    plt.figure(6100+m)
                    a = plt.suptitle(str(obj) + ' Mean Statistical Significance of Binned IP Flux', size=16)
                    a = plt.errorbar(event[j].params.sssteps[1:numsteps+1], ipmean, [ipmean-ipmin,ipmax-ipmean], \
                                     fmt='o-', label=event[j].eventname)
                    plt.xscale('log')
                    plt.yscale('log')
                    a = plt.xlabel('Number of Points per Bin', size=14)
                    a = plt.ylabel('Statistical Significance', size=14)
                    a = plt.legend(loc='best')
                    plt.savefig(event[j].modeldir + "/" + str(obj) + 
                                "-fig" + str(6100+m) + "-" + saveext + ".png")
                    # OPEN PROCESS ID TO SAVE allknots
                    fit[j].allknotfile = event[j].modeldir + "/d-" + event[j].eventname + '-allknots-' + \
                                         fit[j].saveext + '.npy'
                    fit[j].allknotpid  = open(fit[j].allknotfile, 'wb')
                else:
                    if domodel[m,i]:
                        fit[j].fit0 *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
                k += 1

        for j in range(numevents):
            fit[j].redchisq = ( sum((fit[j].fit0 - fit[j].flux)**2 / fit[j].sigma**2) /
                                (fit[j].nobj - fit[j].numfreepars[m])                    )
            print("Reduced Chi^2: " + str(fit[j].redchisq), file=printout)
            redchisq[m] = fit[j].redchisq

        #BIN DATA
        for j in range(numevents):
            fit[j].binfluxuc  = np.zeros(nbins[j])
            fit[j].binphase   = np.zeros(nbins[j])
            fit[j].binphaseuc = np.zeros(nbins[j])
            fit[j].binstduc   = np.zeros(nbins[j])
            fit[j].binfit0    = np.zeros(nbins[j])
            fit[j].binbestfit = np.zeros(nbins[j])
            for i in range(nbins[j]):
                startuc              = int(1.* i   *fit[j].nobjuc/nbins[j])
                enduc                = int(1.*(i+1)*fit[j].nobjuc/nbins[j])
                fit[j].binphaseuc[i] = np.median(fit[j].phaseuc[startuc:enduc])
                fit[j].binfluxuc[i]  = sum(fit[j].fluxuc[startuc:enduc]/fit[j].sigmauc[startuc:enduc]**2)/ \
                    			       sum(1/fit[j].sigmauc[startuc:enduc]**2)
                fit[j].binstduc[i]   = np.sqrt(1 / sum(1/fit[j].sigmauc[startuc:enduc]**2))
                start                = int(1.* i   *fit[j].nobj/nbins[j])
                end                  = int(1.*(i+1)*fit[j].nobj/nbins[j])
                fit[j].binphase[i]   = np.median(fit[j].phase[start:end])
                fit[j].binfit0[i]    = np.mean(fit[j].fit0[start:end])

        # Plot binned data with initial fit
        ebfmt   = ['bo',  'go',  'ro',  'co',  'mo',  'yo',  'ko',  'wo' ]
        pltfmt  = ['b-',  'g-',  'r-',  'c-',  'm-',  'y-',  'k-',  'w-' ]
        pltfmt2 = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--', 'w--']
        plt.figure(5100+m)
        plt.clf()
        for j in range(numevents):
            a = plt.plot(fit[j].binphase, fit[j].binfit0,    'r-',     label='Initial Fit ' + str(j))
            a = plt.errorbar(fit[j].binphaseuc, fit[j].binfluxuc, fit[j].binstduc, fmt=ebfmt[j], \
                         ms=4, lw=1, label='Binned Data ' + str(j))
        a = plt.title(str(obj) + ' Binned Eclipse Data First Minimizer')
        a = plt.xlabel('Orbital Phase')
        a = plt.ylabel('Binned Flux')
        a = plt.text(min(fit[0].phaseuc), max(fit[0].binfluxuc), text)
        a = plt.legend(loc='best')
        a.fontsize = 8
        plt.savefig(event[j].modeldir + "/" + str(obj) + "-fig" + str(5100+m) + "-" + saveext + ".png")


        # New least-squares fit
        # FINDME: the idea next, is to calculate chisq in the vicinity.
        bestp = params
        if leastsq == True:
            print("Re-calculating least-squares fit with new errors.")
            output, err = op.leastsq(modelfunc, bestp[myfreepars[m]], args=(bestp, sigma), factor=100,
                                     ftol=1e-16, xtol=1e-16, gtol=1e-16)
            if (err >= 1) and (err <= 4):
                print("Fit converged without error.")
            else:
                print("WARNING: Error with least squares fit!")
            print("Least squares fit best parameters:")
            print(output)
            bestfit.append(bestp)

        #FIX IP MAPPING TO bestmip IF isfixipmap = TRUE
        for j in range(numevents):
            if hasattr(event[j].params, 'isfixipmap') and event[j].params.isfixipmap == True:
                for i in range(cummodels[j],cummodels[j+1]):
                    if functype[i] == 6:
                        foo = funcx  .pop(i)
                        foo = funcxuc.pop(i)
                        funcx  .append([fit[j].bestmip,   fit[j].binipflux,   np.zeros(len(fit[j].wbfipmask  ))])
                        funcxuc.append([fit[j].bestmipuc, fit[j].binipfluxuc, np.zeros(len(fit[j].wbfipmaskuc))])
                        myfuncs[i] = models.fixipmapping

        # Compute separate best eclipse, ramp, and intrapixel models and best parameters
        for j in range(numevents):
            fit[j].allparams  = params[numparams[cummodels[j]]:numparams[cummodels[j+1]]]
            fit[j].meanknots  = 0
            fit[j].medianknots= 0
            fit[j].stdknots   = 0
            fit[j].bestp      = bestp    [numparams[cummodels[j]]:numparams[cummodels[j+1]]]
            fit[j].bestfit    = np.ones(fit[j].nobj)
            k = 0
            for i in range(cummodels[j],cummodels[j+1]):
                if functype[i] == 6:
                    fit[j].allknotpid.close()
                    del(fit[j].allknotpid)
                    fit[j].etc[k]      = fit[j].bestfit
                if domodel[m,i]:
                    fit[j].bestfit *= myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
                #UPDATE BEST-FIT PARAMETER LIST FOR NEXT RUN   # this is the old output file
#                 if fit[j].modelslist[k] != 'ipspline':
#                     parlist[j][k][2][0] = bestp[iparams[i]]
                k += 1

            fit[j].bestecl      = np.ones(fit[j].nobj)
            fit[j].bestramp     = np.ones(fit[j].nobj)
            fit[j].bestip       = np.ones(fit[j].nobj)
            fit[j].bestpos      = np.ones(fit[j].nobj)
            fit[j].bestvs       = np.ones(fit[j].nobj)
            fit[j].bestff       = np.ones(fit[j].nobj)
            fit[j].bestmip      = np.ones(fit[j].nobj)
            fit[j].bestpip      = []
            fit[j].bestfitmip   = []
            fit[j].bestecluc    = np.ones(fit[j].nobjuc)
            fit[j].bestrampuc   = np.ones(fit[j].nobjuc)
            fit[j].bestipuc     = np.ones(fit[j].nobjuc)
            fit[j].bestposuc    = np.ones(fit[j].nobjuc)
            fit[j].bestvsuc     = np.ones(fit[j].nobjuc)
            fit[j].bestffuc     = np.ones(fit[j].nobjuc)
            fit[j].bestmipuc    = np.ones(fit[j].nobjuc)
            fit[j].bestfitmipuc = []
            if   hasattr(fit[j].i, 'flux'):
                fit[j].systemflux   = fit[j].bestp[fit[j].i.flux]
            elif hasattr(fit[j].i, 'flux2'):
                fit[j].systemflux   = fit[j].bestp[fit[j].i.flux2]
            else:
                fit[j].systemflux   = 1
            k = 0
            for i in range(cummodels[j],cummodels[j+1]):
                if   functype[i] == 0:
                    fit[j].bestecl     *= myfuncs[i](bestp[iparams[i]], funcx[i])
                    fit[j].bestecluc   *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
                elif functype[i] == 1:
                    fit[j].bestramp    *= myfuncs[i](bestp[iparams[i]], funcx[i])
                    fit[j].bestrampuc  *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
                elif functype[i] == 2:
                    fit[j].bestip      *= myfuncs[i](bestp[iparams[i]], funcx[i],   fit[j].etc[k])
                    fit[j].bestipuc    *= myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].etc[k])
                    fit[j].bestpip      = bestp[iparams[i]]
                elif functype[i] == 3:
                    fit[j].bestpos     *= myfuncs[i](bestp[iparams[i]], funcx[i],   fit[j].etc[k])
                    fit[j].bestposuc   *= myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].etc[k])
                elif functype[i] == 4:
                    fit[j].bestvs      *= myfuncs[i](bestp[iparams[i]], funcx[i])
                    fit[j].bestvsuc    *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
                    fit[j].bestpvs      = bestp[iparams[i]]
                elif functype[i] == 5:
                    fit[j].bestff      *= myfuncs[i](bestp[iparams[i]], funcx[i])
                    fit[j].bestffuc    *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
                elif functype[i] == 6:
                    fit[j].bestfitmip   = fit[j].bestecl   * fit[j].bestramp   * fit[j].bestip   * \
                                          fit[j].bestpos   * fit[j].bestvs     * fit[j].bestff
                    fit[j].bestfitmipuc = fit[j].bestecluc * fit[j].bestrampuc * fit[j].bestipuc * \
                                          fit[j].bestposuc * fit[j].bestvsuc   * fit[j].bestffuc
                    a, fit[j].binipflux = myfuncs[i](bestp[iparams[i]], funcx[i], \
                                                     fit[j].bestfitmip, retbinflux=True)
                    b,fit[j].binipfluxuc= myfuncs[i](bestp[iparams[i]], funcxuc[i], \
                                                     fit[j].bestfitmipuc, retbinflux=True)
                    fit[j].bestmip     *= a
                    fit[j].bestmipuc   *= b
                    fit[j].bestmipuc[np.where(fit[j].clipmask)] = fit[j].bestmip
                    fit[j].binipflux    = fit[j].binipflux.  reshape(fit[j].xygrid[0].shape)
                    fit[j].binipfluxuc  = fit[j].binipfluxuc.reshape(fit[j].xygrid[0].shape)
                    fit[j].bestpmip     = bestp[iparams[i]]
                    fit[j].systemflux   = np.mean(fit[j].binipflux.flatten()[np.where(fit[j].binfluxmask)])
                k += 1


        #NEED TO RESTART FOR LOOP HERE!
        rampp             = np.copy(bestp)
        text2 = modelnames[np.where(domodel[m,:])]
        rpars = fparams[m,:][iparams[int(rampidx[m])]]
        text3 = ':'
        for p in np.arange(len(rpars)):
            text3 = text3 + ' ' + str(int(rpars[p]))
        for j in range(numevents):
            fit[j].nump = numparams[cummodels[j+1]] - numparams[cummodels[j]]
            k = 0
            for i in range(cummodels[j],cummodels[j+1]):
                if functype[i] == 6:
                    fit[j].etc[k] = fit[j].bestfitmip
                k += 1

            # Compute residuals and reduced chi-square
            fit[j].residuals   = fit[j].flux - fit[j].bestfit
            residuals.append(fit[j].flux - fit[j].bestfit)
            fit[j].redchisq    = sum((fit[j].residuals)**2 / fit[j].sigma**2)/ (fit[j].nobj - fit[j].numfreepars[m])
            redchisq[m]     = np.sum((fit[j].residuals)**2 / fit[j].sigma**2)/ (fit[j].nobj - fit[j].numfreepars[m])
            print("Best Parameters:", file=printout)
            print(fit[j].bestp, file=printout)
            print("Reduced Chi^2: " + str(fit[j].redchisq), file=printout)


            # Compute bayesian information criterion (BIC)
            fit[j].bic = sum((fit[j].residuals/fit[j].sigma)**2) + fit[j].numfreepars[m]*np.log(fit[j].nobj)
            bic[m]  = np.sum((fit[j].residuals/fit[j].sigma)**2) + fit[j].numfreepars[m]*np.log(fit[j].nobj)
            print("BIC = " + str(fit[j].bic), file=printout)

            # Create data without the eclipse
            fit[j].ramp       = np.ones(fit[j].nobj)
            fit[j].rampuc     = np.ones(fit[j].nobjuc)
            if hasattr(fit[j].i, 'depth'):
                rampp[fit[j].i.depth  + numparams[cummodels[j]]] = 0.0
            if hasattr(fit[j].i, 'depth2'):
                rampp[fit[j].i.depth2 + numparams[cummodels[j]]] = 0.0
            k = 0
            for i in range(cummodels[j],cummodels[j+1]):
                if functype[i] == 6:
                    fit[j].ramp   *= fit[j].bestmip
                    fit[j].rampuc *= fit[j].bestmipuc
                elif domodel[m,i]:
                    fit[j].ramp   *= myfuncs[i](rampp[iparams[i]], funcx  [i], fit[j].etc[k])
                    fit[j].rampuc *= myfuncs[i](rampp[iparams[i]], funcxuc[i], fit[j].etc[k])
                k += 1

            fit[j].normflux      = (fit[j].flux       / fit[j].ramp).flatten()
            fit[j].normsigma     = (fit[j].sigma      / fit[j].ramp).flatten()
            fit[j].normbestfit   = (fit[j].bestfit    / fit[j].ramp).flatten()
            fit[j].normresiduals = (fit[j].residuals  / fit[j].ramp).flatten()
            fit[j].normfluxuc    = (fit[j].fluxuc     / fit[j].rampuc).flatten()
            fit[j].normsigmauc   = (fit[j].sigmauc / fit[j].rampuc).flatten()
            sdnr[m] = np.std((fit[j].residuals  / fit[j].ramp).flatten())

            # Sort flux by x, y and radial positions
            yy      = np.sort(fit[j].position[0])
            xx      = np.sort(fit[j].position[1])
            yflux   = (fit[j].flux/fit[j].bestecl/fit[j].bestramp/fit[j].bestpos/fit[j].bestvs/fit[j].bestff) \
                         [np.argsort(fit[j].position[0])]
            xflux   = (fit[j].flux/fit[j].bestecl/fit[j].bestramp/fit[j].bestpos/fit[j].bestvs/fit[j].bestff) \
                         [np.argsort(fit[j].position[1])]
            ybestip = fit[j].bestip[np.argsort(fit[j].position[0])] * fit[j].bestmip[np.argsort(fit[j].position[0])]
            xbestip = fit[j].bestip[np.argsort(fit[j].position[1])] * fit[j].bestmip[np.argsort(fit[j].position[1])]

            # Sort flux by frmvis
            fvsort = np.sort(fit[j].frmvis)
            vsflux = (fit[j].flux / fit[j].bestecl / fit[j].bestramp / fit[j].bestip / \
                      fit[j].bestpos / fit[j].bestff / fit[j].bestmip)[np.argsort(fit[j].frmvis)]

            # Bin data using weighted average
            fit[j].binmeanfit    = np.zeros(nbins[j])
            fit[j].binmedianfit  = np.zeros(nbins[j])
            fit[j].normbinflux   = np.zeros(nbins[j])
            fit[j].normbinsd     = np.zeros(nbins[j])
            fit[j].normbinbest   = np.zeros(nbins[j])
            fit[j].normbinmean   = np.zeros(nbins[j])
            fit[j].normbinmedian = np.zeros(nbins[j])
            fit[j].binramp       = np.zeros(nbins[j])
            fit[j].binxx         = np.zeros(nbins[j])
            fit[j].binyy         = np.zeros(nbins[j])
            fit[j].binxflux      = np.zeros(nbins[j])
            fit[j].binyflux      = np.zeros(nbins[j])
            fit[j].binxflstd     = np.zeros(nbins[j])
            fit[j].binyflstd     = np.zeros(nbins[j])
            fit[j].binvsflux     = np.zeros(nbins[j])
            fit[j].binvsflstd    = np.zeros(nbins[j])
            fit[j].binfrmvis     = np.zeros(nbins[j])
            fit[j].binxbestip    = np.zeros(nbins[j])
            fit[j].binybestip    = np.zeros(nbins[j])
            fit[j].binxbipstd    = np.zeros(nbins[j])
            fit[j].binybipstd    = np.zeros(nbins[j])
            fit[j].binres        = np.zeros(nbins[j])
            for i in range(nbins[j]):
                startuc                 = int(1.0 * i    * fit[j].nobjuc/nbins[j])
                enduc                   = int(1.0 *(i+1) * fit[j].nobjuc/nbins[j])
                fit[j].binstduc[i]      = np.sqrt(1 / sum(1/fit[j].sigmauc[startuc:enduc]**2)) #
                fit[j].normbinflux[i]   = sum(fit[j].normfluxuc[startuc:enduc] / 
                                               fit[j].normsigmauc[startuc:enduc]**2) / \
                                            sum(1/fit[j].normsigmauc[startuc:enduc]**2)
                fit[j].normbinsd[i]     = np.sqrt(1 / sum(1/fit[j].normsigmauc[startuc:enduc]**2))
                start                   = int(1.0* i   *fit[j].nobj/nbins[j])
                end                     = int(1.0*(i+1)*fit[j].nobj/nbins[j])
                fit[j].binbestfit[i]    = np.mean(fit[j].bestfit[start:end])                   #
                fit[j].binramp[i]       = np.mean(fit[j].ramp[start:end])                      #
                fit[j].binres[i]        = np.mean(fit[j].residuals[start:end])
                fit[j].binvsflux[i]     = np.median(vsflux[start:end])
                fit[j].binvsflstd[i]    = np.std(vsflux[start:end]) / np.sqrt(end-start)
                fit[j].binfrmvis[i]     = np.median(fvsort[start:end])

        # Calculate bic values for joint model fits
        if numevents > 1:
            bic  = 0
            nobj = 0
            for j in range(numevents):
                bic  += sum((fit[j].residuals/event[j].fit[0].sigma)**2)
                nobj += fit[j].nobj
            bic += numfreepars[m]*np.log(nobj)
            print("BIC for joint model fit = " + str(bic), file=printout)

        # Plot binned data, best fit
        plot2 = False
        if plot2:
            plt.figure(5200+m)
            plt.clf()
            for j in range(numevents):
                a = plt.errorbar(fit[j].binphaseuc, fit[j].binfluxuc, fit[j].binstduc, fmt=ebfmt[j], 
                             ms=4, linewidth=1, label='Binned Data')
                a = plt.plot(fit[j].binphase, fit[j].binramp, 'k-', label='No Eclipse')
                a = plt.plot(fit[j].binphase, fit[j].binbestfit,    pltfmt[j], label='Best Fit')
            a = plt.title(obj + ' Binned Eclipse Data With Final Results', size=14)
            a = plt.xlabel('Orbital Phase')
            a = plt.ylabel('Flux')
            a = plt.text(min(fit[0].phaseuc), max(fit[0].binfluxuc), text)
            a = plt.legend(loc='best')
            a.fontsize = 8
            plt.savefig(event[j].modeldir + "/" + obj + "-fig" + str(5200+m) + "-" + saveext + ".png")
    
        # Plot normalized binned and best results
        plt.figure(5300+m, figsize=(9,7))
        plt.clf()
        a = plt.axes([0.11,0.33,0.85,0.6])
        plt.setp(a.get_xticklabels(), visible = False)
        for j in range(numevents):
            a = plt.plot(fit[j].phase,    fit[j].normbestfit,    pltfmt[j], label="Best Fit")
            if hasattr(event[j].params, 'interclip'):
                interclip = event[j].params.interclip
                for i in range(len(interclip)):
                    y = np.ones(interclip[i][1]-interclip[i][0])
                    a = plt.plot(fit[j].phaseuc[interclip[i][0]:interclip[i][1]], y, '-w', lw=2)
            
            a = plt.errorbar(fit[j].binphaseuc, fit[j].normbinflux, fit[j].normbinsd, fmt=ebfmt[j], 
                             ms=5, linewidth=1, label='Normalized Data')

        a = plt.title(obj + ' Binned Data, Best Eclipse', size=14)
        a = plt.ylabel('Normalized flux')
        a = plt.legend(loc='lower right')
        a.fontsize = 8
        xmin, xmax      = plt.xlim()
        # Plot Residuals
        plt.axes([0.11,0.08,0.85,0.22])
        plt.xlabel('Orbital Phase')
        plt.ylabel('Residuals',size=14)
        plt.plot([-1,1],[0,0],':', color='0.4')
        plt.plot(fit[j].binphaseuc, fit[j].binres/fit[j].mflux,   'k+', mew=1.2)
        plt.xlim(xmin,xmax)
        # Print info into final plots
        plt.annotate(str(text2)+text3,                                (0.13,0.89), textcoords='figure fraction')
        plt.annotate("BIC    = " + str(fit[j].bic),                   (0.12,0.40), textcoords='figure fraction')
        plt.annotate("SDNR = "   + str(np.std(fit[j].normresiduals)), (0.12,0.34), textcoords='figure fraction')

    # Modify sigma such that reduced chi-square = 1
    # Spitzer over-estimates errors
    minredchisq = np.amin(redchisq)
    newsigma = fit[0].sigma * np.sqrt(minredchisq)
    for m in np.arange(numleasts):
        bic2[m] = np.sum((residuals[m]/newsigma)**2) + fit[j].numfreepars[m]*np.log(fit[j].nobj)
        plt.figure(5300+m)
        plt.annotate("BIC    = " + str(bic2[m])+" (modified)", (0.12,0.365), textcoords='figure fraction')
        plt.savefig(event[j].modeldir + "/" + obj + "-fig" + str(5300+m) + "-" + saveext + ".png")

    # Write to file the best fit values
    for j in np.arange(numevents):
        for m in np.arange(numleasts):
            fout[j].write("*** Model " + str(m) + " ***\n")
            for k in np.arange(len(modelnames)):
                if domodel[m,k]:
                    fout[j].write(modelnames[k]+'\n')
                    fout[j].write(formatpar(fparams[m,:][iparams[k]])+'\n')
                    fout[j].write(formatpar(bestfit[m][iparams[k]])+'\n')
            fout[j].write('\nBIC             BIC2            SDNR\n')
            fout[j].write('%14.9f'%bic[m] + '  %14.9f'%bic2[m] + '  %.11e\n'%sdnr[m])
        fout[j].write('\nminimum reduced chisq: ' + str(minredchisq) + '\n')
        fout[j].close()


    print('\nReturn here\n')
    return
    if False:

        # Plot visit sensitivity and model
        if 4 in functype:
            numvsmodels = 0
            for i in functype:
                if i == 4:
                    numvsmodels += 1

            plt.figure(607+num*numfigs)
            plt.clf()
            a = plt.suptitle(obj + ' Visit # vs. Sensitivity', size=16)
            k = 1
            for j in range(numevents):
                for i in range(cummodels[j],cummodels[j+1]):
                    if functype[i] == 4:
                        a = plt.subplot(numvsmodels,1,k)
                        a = plt.errorbar(fit[j].binfrmvis, fit[j].binvsflux, fit[j].binvsflstd, fmt='go', label='Data')
                        a = plt.plot(fit[j].frmvis, fit[j].bestvs, 'k.', label='Model')
                        for model in fit[j].modelslist:
                            if (model == 'vsspline'):
                                a = plt.plot(event[j].params.vsknots, fit[j].bestpvs, 'bs')
                        a = plt.ylabel('Flux Sensitivity')
                        a = plt.xlabel('Visit Number')
                        a = plt.legend(loc='best')
                        a.fontsize=8

            plt.savefig(event[j].modeldir + "/" + obj + "-fig" + str(num*numfigs+607) + "-" + saveext + ".png")

        #CONTOUR PLOT OF INTERPOLATED INTRA-PIXEL
        if fit[j].isipmapping:
            plt.ioff()
            plt.figure(608+num*numfigs, figsize=(14, 8))
            plt.clf()
            a = plt.suptitle(obj + ' Filled Contour Plots', size=16)
            palette = plt.matplotlib.colors.LinearSegmentedColormap('jet3',plt.cm.datad['jet'],16384)
            palette.set_under(alpha=0.0)
            plt.subplots_adjust(left=0.08,right=0.96,bottom=0.10,top=0.90,hspace=0.2)
            for j in range(numevents):
                for i in range(cummodels[j],cummodels[j+1]):
                    if functype[i] == 6:
                        k   = 2*j+1
                        vmin = fit[j].binipflux[np.where(fit[j].binipflux > 0)].min()
                        vmax = fit[j].binipflux.max()
                        xmin = fit[j].xygrid[0].min()
                        xmax = fit[j].xygrid[0].max()
                        ymin = fit[j].xygrid[1].min()
                        ymax = fit[j].xygrid[1].max()
                        if fit[j].modelslist.__contains__('nnint'):
                            interp = 'nearest'
                        else:
                            interp = 'bilinear'
                        #PLOT FLUX VS POSITION
                        a = plt.subplot(numevents,2,k)
                        a = plt.imshow(fit[j].binipflux, cmap=palette, vmin=vmin, vmax=vmax, origin='lower', 
                               extent=(xmin,xmax,ymin,ymax), aspect='auto', interpolation=interp)
                        a = plt.colorbar(a)
                        if j == 0:
                            a = plt.title('Flux vs. Position')
                        if j == int(numevents/2):
                            a = plt.ylabel('Relative Pixel Position in y')
                        if j == (numevents-1):
                            a = plt.xlabel('Relative Pixel Position in x')
                        if ymin < -0.5:
                            a = plt.hlines(-0.5, xmin, xmax, 'k')
                        if ymax >  0.5:
                            a = plt.hlines( 0.5, xmin, xmax, 'k')
                        if xmin < -0.5:
                            a = plt.vlines(-0.5, ymin, ymax, 'k')
                        if xmax >  0.5:
                            a = plt.vlines( 0.5, ymin, ymax, 'k')
                        a = plt.text(xmin+0.01, ymin+0.01, event[j].eventname, fontsize=10)

                        #PLOT # OF POINTS/BIN VS POSITION
                        lenbinflux = np.zeros(len(fit[j].wbfipmask))
                        for m in range(len(fit[j].wbfipmask)):
                            lenbinflux[m] = len(fit[j].wbfipmask[m])
                        lenbinflux = lenbinflux.reshape(fit[j].xygrid[0].shape)
                        a = plt.subplot(numevents,2,k+1)
                        a = plt.imshow(lenbinflux, cmap=palette, vmin=1, vmax=lenbinflux.max(), origin='lower', 
                                       extent=(xmin,xmax,ymin,ymax), aspect='auto', interpolation=interp)
                        a = plt.colorbar(a)
                        if j == 0:
                            a = plt.title('# of Points/Bin vs. Position')
                        if j == int(numevents/2):
                            a = plt.ylabel('Relative Pixel Position in y')
                        if j == (numevents-1):
                            a = plt.xlabel('Relative Pixel Position in x')
                        if ymin < -0.5:
                            a = plt.hlines(-0.5, xmin, xmax, 'k')
                        if ymax >  0.5:
                            a = plt.hlines( 0.5, xmin, xmax, 'k')
                        if xmin < -0.5:
                            a = plt.vlines(-0.5, ymin, ymax, 'k')
                        if xmax >  0.5:
                            a = plt.vlines( 0.5, ymin, ymax, 'k')
            
            plt.savefig(event[j].modeldir + "/" + obj + "-fig" + str(num*numfigs+608) + "-" + saveext + ".png")
            plt.ion()
            
            '''
                #HISTOGRAM OF POINTS PER BIN
                lenbinflux = np.zeros(len(fit[j].wherebinflux))
                for i in range(len(fit[j].wherebinflux)):
                    lenbinflux[i] = len(fit[j].wherebinflux[i])
                a = plt.subplot(numevents,2,k+1)
                numbins = 59
                a = np.histogram(lenbinflux, bins=[1,numbins+1,lenbinflux.max()])
                gtnumbins = a[0][1]
                a = plt.hist(lenbinflux, bins=numbins, range=(1,numbins+1), align='left')
                y = a[0].max()/2
                a = plt.ylabel('Number of Bins')
                a = plt.xlabel('Number of Points per Bin')
                a = plt.text(numbins/2, y, str(gtnumbins)+' bins > '+str(numbins)+' points/bin')
            '''
        plt.show()
        
        return
