
"""
 MODIFICATION HISTORY:
    Written by:	Kevin Stevenson, UCF  	2008-07-02
                kevin218@knights.ucf.edu
    Finished initial version:       kevin   2008-09-08
    Updated for multi events:       kevin   2009-11-01
    Added ip interpolation:         kevin   2010-06-28
    Added minnumpts & fixipmap:     kevin   2010-08-03
    Least-squares w/ multi events:  kevin   2010-08-19
    Use phase or bjd as time unit   kevin   2010-11-17
    Minor improvements from pipeline meeting
                                    kevin   2010-12-11
    Added Ballard IP method         kevin   2010-12-31
    Added sinusoidal-type models    kevin   2011-01-25
    Added bjdutc and bjdtdb         kevin   2011-03-19
    Write allparams to file         kevin   2011-03-31
    Added autocorr & 2D hist plots  kevin   2011-07-18
    Added param orthogonalization   kevin   2011-07-26
    Added priors                    kevin   2011-10-11
    Non-fixed priors in least-squares       2014-08-21
    Added DEMCz                     kevin   2014-08-21
    Converted to Python3            kbs     2018-11-15
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d
import time, os, sys, importlib
import models_c as mc
import mcmc
import demc
import residperm
import paramedit as pe
import readeventhdf
import smoothing
import correlated_noise as cn
import plots
import nasc

def rundmc(event, num=0, printout=sys.stdout, isinteractive=True):
    """
    This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm
    
    Parameters
    ----------
    event       : list
        List of event objects
    num         : integer
        Iteration number for runs with multiple models
    printout    : object
        File object for directing print statements
    
    Returns
    -------
    None
    
    """
    #When allplots=False: the correlation and histogram figures are not produced.
    #When normflux=Flase: flux is not normalized w.r.t. each position's median value.
    #                     This should only be used when using normflux model.
    numit        = event[0].params.numit.astype(np.int)
    #nbins        = event[0].params.nbins
    chi2flag     = event[0].params.chi2flag
    allplots     = event[0].params.allplots
    leastsq      = event[0].params.leastsq
    #normflux     = event[0].params.normflux
    nbins        = []
    normflux     = []
    obj          = event[0].filename
    numevents    = len(event)
    nummodels    = np.zeros(numevents, dtype=int)
    fit          = []
    params       = []
    pmin         = []
    pmax         = []
    stepsize     = []
    numparams    = [0]
    iparams      = []
    myfuncs      = []
    functype     = []
    parlist      = []
    numfigs      = 20
    pcounter     = 0
    isnoise      = False #FINDME: test, should be False
    if hasattr(event[0].params, 'gamma'):
        gamma    = event[0].params.gamma[num]
    else:
        gamma    = None
    if hasattr(event[0].params, 'ipfactor'):
        ipfactor = event[0].params.ipfactor
    else:
        ipfactor = 1.
    if hasattr(event[0].params, 'isblock'):
        isblock = event[0].params.isblock
    else:
        isblock = False
    if hasattr(event[0].params, 'isresidperm'):
        isresidperm = event[0].params.isresidperm
        nights      = np.zeros(numevents, dtype=int)
        for j in range(numevents):
            nights[j] = event[j].params.night
    else:
        isresidperm = False
    if hasattr(event[0].params, 'night'):
        nights      = np.zeros(numevents, dtype=int)
        for j in range(numevents):
            nights[j] = event[j].params.night
        #nnights     = len(np.unique(nights))
    else:
        nights      = np.zeros(numevents, dtype=int)
    for j in range(numevents):
        if hasattr(event[j].params, 'nbins'):
            pass
        else:
            event[j].params.nbins = int(sum(event[j].good))
    
    for j in range(numevents):
        fit.append(readeventhdf.fits())
        event[j].fit.append(fit[j])
        fit[j].model   = event[j].params.model[num]
        if hasattr(event[j].params, 'blocks'):
            fit[j].blocks = event[j].params.blocks[num]
        
        #UPDATE event.params
        event[j].params.numit    = np.copy(numit)
        event[j].params.chi2flag = chi2flag
        nbins   .append(event[j].params.nbins)
        normflux.append(event[j].params.normflux)
        if hasattr(event[j].params, 'nchains') == False:
            event[j].params.nchains = [20,5]
        if isinstance(event[j].params.nchains, int):
            event[j].params.nchains = [event[j].params.nchains,event[j].params.nchains]
        if hasattr(event[j].params, 'newortho') == False:
            event[j].params.newortho = False
        event[j].isortho = False
        
        #READ IN DATA
        fit[j].mflux = np.mean(event[j].aplev[np.where(event[j].good == 1)])
        #NORMALIZE FLUX AT EACH POSITION IF USING > 1 POSITION
        if ((event[j].npos > 1) and (normflux[j] == True)):
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
        if hasattr(event[j], 'bjdutc') == False:
            event[j].bjdutc = event[j].juldat + event[j].bjdcor/86400.
        if hasattr(event[j], 'bjdtdb') == False:
            event[j].bjdtdb = event[j].bjdutc
            print('****WARNING: BJD_TDB does not exist. Using BJD_UTC.')
        fit[j].bjdutcuc = event[j].bjdutc[np.where(event[j].good == 1)]
        fit[j].bjdtdbuc = event[j].bjdtdb[np.where(event[j].good == 1)]
        fit[j].posuc    = event[j].pos   [np.where(event[j].good == 1)]
        fit[j].frmvisuc = event[j].frmvis[np.where(event[j].good == 1)]
        #if hasattr(event[j].fp, 'P_hat'):
        #    fit[j].P_hatuc  = event[j].fp.P_hat[np.where(event[j].good == 1)]
        if hasattr(event[j], 'apdata'):
            fit[j].apdatauc = event[j].apdata[np.where(event[j].good[0] == 1)]
        
        #SORT DATA BY bjdutcuc IF USING > 1 POSITION
        if event[j].npos > 1:
            isort           = np.argsort(fit[j].bjdutcuc)
            fit[j].phaseuc  = fit[j].phaseuc[isort]
            fit[j].bjdutcuc = fit[j].bjdutcuc[isort]
            fit[j].bjdtdbuc = fit[j].bjdtdbuc[isort]
            fit[j].fluxuc   = fit[j].fluxuc[isort]
            fit[j].sigmauc  = fit[j].sigmauc[isort]
            fit[j].yuc      = fit[j].yuc[isort]
            fit[j].xuc      = fit[j].xuc[isort]
            fit[j].time     = fit[j].time[isort]
            fit[j].posuc    = fit[j].posuc[isort]
            fit[j].frmvisuc = fit[j].frmvisuc[isort]
            if hasattr(event[j], 'apdata'):
                fit[j].apdatauc = fit[j].apdatauc[isort]
            #if hasattr(event[j].fp, 'P_hat'):
            #    fit[j].P_hatuc  = fit[j].P_hatuc[isort]
        
        #OBTAIN INITIAL MODEL PARAMETERS
        nummodels[j] = len(fit[j].model)
        if event[j].params.modelfile == None:
            modelfile = event[j].initvalsfile
        elif event[j].params.modelfile.__class__ == str:
            modelfile = event[j].ancildir + event[j].params.modelfile
        elif len(event[j].params.modelfile) == len(event[j].params.model):
            modelfile = event[j].ancildir + event[j].params.modelfile[num] 
        else:
            modelfile = event[j].ancildir + event[j].params.modelfile[0]
        parlist.append(pe.read(modelfile, fit[j].model, event[j]))
        for i in range(nummodels[j]):
            pars       = parlist[j][i][2]
            if fit[j].model[i] == 'ipspline':
                #READ NUMBER OF KNOTS ALONG x AND y
                numipyk, numipxk = pars[0]
                #TEMPORARILY POPULATE INITIAL INTRA-PIXEL PARAMETERS WITH ones
                temppars = [np.ones(numipxk*numipyk)*1.0,
                            np.ones(numipxk*numipyk)*pars[1][0],
                            np.ones(numipxk*numipyk)*pars[2][0],
                            np.ones(numipxk*numipyk)*pars[3][0]]
                pars     = temppars
            
            params     = np.concatenate((params,   pars[0]),0)
            pmin       = np.concatenate((pmin,     pars[1]),0)
            pmax       = np.concatenate((pmax,     pars[2]),0)
            stepsize   = np.concatenate((stepsize, pars[3]),0)
            numparams  = np.concatenate((numparams, [numparams[-1] + len(pars[0])]),0)
            iparams.append(list(range(pcounter,pcounter+len(pars[0]))))
            
            #CHECK FOR SHARED PARAMETERS, MAKE APPROPRIATE CHANGES TO PARAMETER VALUES AND LIMITS
            for k in range(len(pars[0])):
                if pars[3][k] < 0:
                    params [pcounter+k] = params [int(-pars[3][k]-1)]
                    pmin   [pcounter+k] = pmin   [int(-pars[3][k]-1)]
                    pmax   [pcounter+k] = pmax   [int(-pars[3][k]-1)]
                    iparams[-1][k]      = int(-pars[3][k]-1)
            
            pcounter += len(pars[0])
        
        #SETUP MODELS BY DECLARING FUNCTIONS, FUNCTION TYPES, PARAMETER NAMES, INDICES & SAVE EXTENSIONS
        func, funct, fit[j].parname, fit[j].i, fit[j].saveext = mc.setupmodel(fit[j].model, fit[j].i)
        myfuncs  = np.concatenate(( myfuncs, func), 0)
        for i in range(len(funct)):
            functype.append(funct[i])
        if np.where(funct == 'ipmap') != (nummodels[j] - 1):
            print("ERROR: The interpolation model must be the last model listed.")
            return
        if np.where(funct == 'gp') != (nummodels[j] - 1):
            print("ERROR: Gaussian process model must be the last model listed.")
            return
    
    #NEED TO RESTART FOR LOOP HERE!
    cummodels    = [0]      #Cumulative list of number of models
    funcx        = []
    funcxuc      = []
    ifreepars    = np.array([], dtype=int)
    nonfixedpars = np.array([], dtype=int)
    isprior      = np.array([], dtype=int)
    text         = ""
    maxnumfp     = 0
    numfreepars  = np.array(np.where(stepsize > 0)).flatten().size
    isipmapping  = False
    print("\nCurrent event & model:", file=printout)
    for j in range(numevents):
        print(event[j].eventname, file=printout)
        print(fit[j].model, file=printout)
        cummodels.append(int(cummodels[j] + nummodels[j]))
        fit[j].nump       = numparams[cummodels[j+1]] - numparams[cummodels[j]]
        
        #CHECK FOR PRIORS, ASSIGN INDICES AND priorvals=[MEAN,LOW,HIGH] TO fit
        fit[j].ipriors   = []
        fit[j].priorvals = []
        fit[j].isprior   = np.zeros(fit[j].nump)
        if hasattr(event[j].params, 'priorvars') and len(event[j].params.priorvars) > 0:
            i = 0
            for pvar in event[j].params.priorvars:
                ipvar = getattr(fit[j].i, pvar) + numparams[cummodels[j]]
                fit[j].ipriors.append(ipvar)
                fit[j].isprior[ipvar-numparams[cummodels[j]]] = 1
                params[ipvar] = event[j].params.priorvals[i][0]
                i += 1
            fit[j].priorvals = event[j].params.priorvals
        isprior = np.concatenate((isprior, fit[j].isprior),0)
        
        #SPECIFY INDICES OF FREE AND NONPRIOR PARAMETERS
        fit[j].ifreepars    = np.array(np.where(stepsize[numparams[cummodels[j]]:numparams[cummodels[j+1]]] > 0)).flatten()
        ifreepars = np.concatenate((ifreepars, fit[j].ifreepars + numparams[cummodels[j]]),0)
        fit[j].nonfixedpars = np.array(np.where(stepsize[numparams[cummodels[j]]:numparams[cummodels[j+1]]] != 0)).flatten()
        nonfixedpars = np.concatenate((nonfixedpars, fit[j].nonfixedpars + numparams[cummodels[j]]),0)
        fit[j].numfreepars = fit[j].ifreepars.size
        maxnumfp = np.max((maxnumfp, fit[j].numfreepars))
        
        #DETERMINE THE CENTROID LOCATION RELATIVE TO THE CENTER OF PIXEL
        #RECORD IN WHICH QUADRANT THE CENTER OF LIGHT FALLS
        fit[j].nobjuc   = fit[j].fluxuc.size
        fit[j].quadrant = np.zeros(fit[j].nobjuc)
        fit[j].numq     = np.zeros(4)
        yround          = np.round(np.median(fit[j].yuc))
        xround          = np.round(np.median(fit[j].xuc))
        fit[j].y        = fit[j].yuc - yround
        fit[j].x        = fit[j].xuc - xround
        print("Positions are with respect to (y,x): " + str(int(yround)) + ", " + str(int(xround)), file=printout)
        #DETERMINE PIXEL QUADRANT BASED ON ydiv & xdiv BOUNDARY CONDITIONS
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
        
        #CLIP FIRST AND LAST FEW PTS OF DATA
        if len(event[j].params.preclip) == len(event[j].params.model):
            preclip     = event[j].params.preclip[num]
        else:
            preclip     = event[j].params.preclip[0]
        if len(event[j].params.postclip) == len(event[j].params.model):
            postclip    = fit[j].nobjuc - event[j].params.postclip[num]
        else:
            postclip    = fit[j].nobjuc - event[j].params.postclip[0]
        
        #DEFINE clipmask WHERE 1 IS KEPT AND 0 IS CLIPPED
        fit[j].clipmask = np.zeros(fit[j].nobjuc)
        fit[j].clipmask[preclip:postclip] = 1
        #USE interclip IF IT EXISTS TO CLIP POINTS IN THE MIDDLE
        if hasattr(event[j].params, 'interclip'):
            interclip   = event[j].params.interclip
            for i in range(len(interclip)):
                fit[j].clipmask[interclip[i][0]:interclip[i][1]] = 0

        #DEFINE ipmask WHERE 1 IS USED IN THE INTRAPIXEL INTERPOLATION MODEL AND 0 IS NOT
        fit[j].ipmaskuc = np.ones(fit[j].nobjuc)
        if hasattr(event[j].params, 'ipclip'):
            ipclip = event[j].params.ipclip
            for i in range(len(ipclip)):
                fit[j].ipmaskuc[ipclip[i][0]:ipclip[i][1]] = 0
        
        #DEFINE INTRAPIXEL MASK, (y,x) POSITION & NUMBER OF POINTS AFTER CLIPPING
        fit[j].ipmask   = np.copy(fit[j].ipmaskuc[np.where(fit[j].clipmask)])
        fit[j].position = np.array([fit[j].y[np.where(fit[j].clipmask)], 
                                    fit[j].x[np.where(fit[j].clipmask)], 
                                    fit[j].quadrant[np.where(fit[j].clipmask)]])
        fit[j].nobj     = fit[j].position[0].size
        
        #CALCULATE minnumptsmask = 1 FOR AT LEAST minnumpts IN EACH BIN
        fit[j].minnumptsmask = np.ones(fit[j].nobj, dtype=int)
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'ipmap':
                #Read bin sizes
                if len(event[j].params.ystep) == len(event[j].params.model):
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
                #Number of bins in y,x dimensions
                ysize    = int((ymax-ymin)/ystep + 1)
                xsize    = int((xmax-xmin)/xstep + 1)
                ygrid, ystep  = np.linspace(ymin, ymax, ysize, retstep=True)
                xgrid, xstep  = np.linspace(xmin, xmax, xsize, retstep=True)
                #Minimum number of acceptable points in a bin
                if event[j].params.minnumpts.__class__ == int:
                    minnumpts = event[j].params.minnumpts
                elif len(event[j].params.minnumpts) == len(event[j].params.model):
                    minnumpts = event[j].params.minnumpts[num]
                else:
                    minnumpts = event[j].params.minnumpts[0]
                #CALCULATE BINS FOR 2D BINNING
                #Gets the number of points in each bin.
                for m in range(ysize):
                    wbftemp   = np.where(np.abs(fit[j].position[0]-ygrid[m])-ystep/2. <= 1e-16)[0]
                    for n in range(xsize):
                        wbf       = wbftemp[np.where((np.abs(fit[j].position[1,[wbftemp]]-xgrid[n])-xstep/2. <= 1e-16)[0])]
                        wbfipmask = wbf    [np.where(fit[j].ipmask[wbf] == 1)]
                        if len(wbfipmask) < minnumpts:
                            fit[j].minnumptsmask[wbf] = 0

        #REDEFINE CLIPPED VARIABLES BASED ON minnumpts FOR IP MAPPING
        fit[j].clipmask[np.where(fit[j].clipmask)] *= fit[j].minnumptsmask
        fit[j].isclipmask = np.where(fit[j].clipmask)[0]
        fit[j].phase    = fit[j].phaseuc [fit[j].isclipmask]
        fit[j].bjdutc   = fit[j].bjdutcuc[fit[j].isclipmask]
        fit[j].bjdtdb   = fit[j].bjdtdbuc[fit[j].isclipmask]
        fit[j].flux     = np.copy(fit[j].fluxuc  [fit[j].isclipmask])
        fit[j].sigma    = np.copy(fit[j].sigmauc [fit[j].isclipmask])
        fit[j].pos      = np.copy(fit[j].posuc   [fit[j].isclipmask])
        fit[j].frmvis   = np.copy(fit[j].frmvisuc[fit[j].isclipmask])
        fit[j].ipmask   = np.copy(fit[j].ipmaskuc[fit[j].isclipmask])
        #FINDME: may not need with new pipeline.
        if hasattr(event[j], 'apdata'):
            fit[j].apdata   = np.copy(fit[j].apdatauc[fit[j].isclipmask])
        fit[j].position = np.array([fit[j].y[fit[j].isclipmask], 
                                    fit[j].x[fit[j].isclipmask], 
                                    fit[j].quadrant[fit[j].isclipmask]])
        fit[j].nobj     = fit[j].flux.size
        fit[j].positionuc = np.array([fit[j].y, fit[j].x, fit[j].quadrant])
        
        #DETERMINE BIN LOCATION FOR EACH POSITION
        fit[j].isipmapping = False
        fit[j].numknots = 0
        fit[j].ballardip   = np.ones(fit[j].flux.size)
        fit[j].ballardipuc = np.ones(fit[j].fluxuc.size)
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'ipmap':
                isipmapping           = True
                fit[j].isipmapping    = True    #Using mapping for IP sensitivity
                fit[j].wherebinflux   = []      #List of size = # of bins, def which points fall into each bin
                fit[j].wherebinfluxuc = []      #Un-clipped version of above
                fit[j].wbfipmask      = []      #Intrapixel masked version of wherebinflux
                fit[j].wbfipmaskuc    = []      #Un-clipped version of above
                fit[j].binloc         = np.zeros((2, fit[j].nobj),   dtype=int) - 1
                fit[j].binlocuc       = np.zeros((2, fit[j].nobjuc), dtype=int) - 1
                #Read bin sizes
                if len(event[j].params.ystep) == len(event[j].params.model):
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
                '''
                ymin     = np.round(fit[j].y.min(),3) - 3*ystep
                ymax     = np.round(fit[j].y.max(),3) + 3*ystep
                xmin     = np.round(fit[j].x.min(),3) - 3*xstep
                xmax     = np.round(fit[j].x.max(),3) + 3*xstep
                '''
                #Number of bins in y,x dimensions
                ysize    = int((ymax-ymin)/ystep + 1)
                xsize    = int((xmax-xmin)/xstep + 1)
                ygrid, ystep  = np.linspace(ymin, ymax, ysize, retstep=True)
                xgrid, xstep  = np.linspace(xmin, xmax, xsize, retstep=True)
                #ystep    = np.round(ystep,8)
                #xstep    = np.round(xstep,8)
                fit[j].xygrid = np.meshgrid(xgrid, ygrid)
                fit[j].binfluxmask   = np.zeros((ysize, xsize), dtype=int)
                fit[j].binfluxmaskuc = np.zeros((ysize, xsize), dtype=int)
                fit[j].numpts        = np.zeros((ysize, xsize)) #Number of points per bin
                fit[j].numptsuc      = np.zeros((ysize, xsize))
                #Minimum number of acceptable points in a bin
                if event[j].params.minnumpts.__class__ == int:
                    minnumpts = event[j].params.minnumpts
                elif len(event[j].params.minnumpts) == len(event[j].params.model):
                    minnumpts = event[j].params.minnumpts[num]
                else:
                    minnumpts = event[j].params.minnumpts[0]
                print('Step size in y = ' + str(ystep), file=printout)
                print('Step size in x = ' + str(xstep), file=printout)
                print('Ignoring bins with < ' + str(minnumpts) + ' points.', file=printout)
                print('Computing bin for each position.')
                #ASSIGN BINS FOR 2D BINNING
                for m in range(ysize):
                    wbftemp   = np.where(np.abs(fit[j].position[0]-ygrid[m])-ystep/2. <= 1e-16)[0]
                    wbftempuc = np.where(np.abs(fit[j].y          -ygrid[m])-ystep/2. <= 1e-16)[0]
                    for n in range(xsize):
                        wbf   = wbftemp[np.where((np.abs(fit[j].position[1,[wbftemp]]-xgrid[n])-xstep/2. <= 1e-16)[0])]
                        wbfuc = wbftempuc[np.where(np.abs(fit[j].x[wbftempuc]-xgrid[n])-xstep/2. <= 1e-16)]
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
                        # do the same for unclipped data
                        if len(wbfipmaskuc) > 0:
                            fit[j].binfluxmaskuc[m,n] = 1
                            fit[j].numptsuc     [m,n] = len(wbfipmaskuc)
                            fit[j].binlocuc[0, wbfuc] = m*xsize + n
                            fit[j].wherebinfluxuc.append(wbfuc)
                            fit[j].wbfipmaskuc.append(wbfipmaskuc)
                        else:
                            fit[j].wherebinfluxuc.append([])
                            fit[j].wbfipmaskuc.append([])
                
                #Check for points not associated with a bin
                isnobins    = np.where(fit[j].binloc  [0] == -1)[0]
                isnobinsuc  = np.where(fit[j].binlocuc[0] == -1)[0]
                if len(isnobins) > 0 or len(isnobinsuc) > 0:
                    print('****WARNING: Not all points have been assigned to bins.****')
                    print('Terminating run.')
                    return
                fit[j].numknots = np.sum(fit[j].binfluxmask)
                #Read smoothing paramters
                if len(event[j].params.nx) == len(event[j].params.model):
                    ny = event[j].params.ny[num]
                    nx = event[j].params.nx[num] 
                else:
                    ny = event[j].params.ny[0]
                    nx = event[j].params.nx[0] 
                if len(event[j].params.sx) == len(event[j].params.model):
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
                                                     fit[j].position[0] <= ygrid[m+1]))[0]
                    for n in range(xsize-1):
                        wherexy = wherey[np.where(np.bitwise_and(fit[j].position[1,[wherey]] > xgrid[n  ],
                                                                 fit[j].position[1,[wherey]] <= xgrid[n+1])[0])[0]]
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
                                                     fit[j].y <= ygrid[m+1]))[0]
                    for n in range(xsize-1):
                        wherexy = wherey[np.where(np.bitwise_and(fit[j].x[wherey] > xgrid[n  ],
                                                                 fit[j].x[wherey] <= xgrid[n+1]))[0]]
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
                
                # COMBINE MODEL PARAMETERS INTO ONE LIST               
                fit[j].posflux  = [fit[j].y[fit[j].isclipmask], 
                                   fit[j].x[fit[j].isclipmask], 
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
            
            if functype[i] == 'ballardip':
                print("Computing intrapixel effect")
                ipparams = [params[fit[j].i.sigmay], params[fit[j].i.sigmax], int(params[fit[j].i.nbins])]
                position = [fit[j].position, fit[j].ballardip, fit[j].flux]
                fit[j].ballardip = mc.ballardip(ipparams, position , etc=[fit[j].ballardip])
                #position = [fit[j].positionuc, fit[j].ballardipuc, fit[j].fluxuc]
                #fit[j].ballardipuc = mc.ballardip(ipparams, position, etc=[True])
                '''
                for i in range(nbinstemp):
                    start   = int(1.*i*nobj/nbinstemp)
                    end     = int(1.*(i+1)*nobj/nbinstemp)
                    s       = np.ones(nobj)
                    s[start:end] = 0
                    s[np.where(fit[j].phase >= (params[fit[j].i.midpt] - params[fit[j].i.width]/2.))[0][0]:  \
                      np.where(fit[j].phase >= (params[fit[j].i.midpt] + params[fit[j].i.width]/2.))[0][0]] = 0
                    biny = np.mean(y[start:end])
                    binx = np.mean(x[start:end])
                    weight[start:end] = sum(np.exp(-0.5*((x-binx)/event[j].xprecision)**2) * \
                                            np.exp(-0.5*((y-biny)/event[j].yprecision)**2)*flux*s) / \
                                        sum(np.exp(-0.5*((x-binx)/event[j].xprecision)**2) * \
                                            np.exp(-0.5*((y-biny)/event[j].yprecision)**2)*s)
                    start   = int(1.*i*nobjuc/nbinstemp)
                    end     = int(1.*(i+1)*nobjuc/nbinstemp)
                    s       = np.ones(nobjuc)
                    s[start:end] = 0
                    s[np.where(fit[j].phaseuc >= (params[fit[j].i.midpt] - params[fit[j].i.width]/2.))[0][0]:  \
                      np.where(fit[j].phaseuc >= (params[fit[j].i.midpt] + params[fit[j].i.width]/2.))[0][0]] = 0
                    biny = np.mean(yuc[start:end])
                    binx = np.mean(xuc[start:end])
                    weightuc[start:end] = sum(np.exp(-0.5*((xuc-binx)/event[j].xprecision)**2) * \
                                            np.exp(-0.5*((yuc-biny)/event[j].yprecision)**2)*fluxuc*s) / \
                                        sum(np.exp(-0.5*((xuc-binx)/event[j].xprecision)**2) * \
                                            np.exp(-0.5*((yuc-biny)/event[j].yprecision)**2)*s)
                '''
    
    #RESTART FOR LOOP HERE
    inonprior  = np.where(np.bitwise_and(stepsize > 0,isprior == 0))[0] #Indices of non-prior parameters
    iortholist = []                                                     #List of variable indices to orthogonalize
    opname     = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9']   #Orthogonal parameter names
    #Assigns phase/bjdutc/bjdtdb or x,y-position as independent variable for model fitting
    for j in range(numevents):
        #Use orbits or days as unit of time.
        if hasattr(event[j].params, 'timeunit') and ((event[j].params.timeunit == 'days')
                                                 or  (event[j].params.timeunit == 'days-utc')):
            fit[j].timeunit   = fit[j].bjdutc   - event[j].params.tuoffset
            fit[j].timeunituc = fit[j].bjdutcuc - event[j].params.tuoffset
        elif hasattr(event[j].params, 'timeunit') and (event[j].params.timeunit == 'days-tdb'):
            fit[j].timeunit   = fit[j].bjdtdb   - event[j].params.tuoffset
            fit[j].timeunituc = fit[j].bjdtdbuc - event[j].params.tuoffset
        else:
            fit[j].timeunit   = fit[j].phase
            fit[j].timeunituc = fit[j].phaseuc
            event[j].params.timeunit = 'orbits'
            event[j].params.tuoffset = 0.0
        
        # Initialize ortho variables
        if hasattr(event[j].params, 'ortholist'):
            if len(event[j].params.ortholist) == len(event[j].params.model):
                ortholist = event[j].params.ortholist[num] 
            else:
                ortholist = event[j].params.ortholist[0]
        else:
            ortholist = []
        diminvtrans       = len(ortholist)
        fit[j].invtrans   = np.matrix(np.identity(diminvtrans)) #Inverse transformation matrix, initialized to 1
        fit[j].trans      = np.matrix(np.identity(diminvtrans)) #Transformation matrix, initialized to 1
        fit[j].origin     = np.zeros(diminvtrans)               #Origin of coordinate system, initialized to 0
        fit[j].orthosigma = np.ones(diminvtrans)                #Normalizing uncertainties, initialized to 1
        iortholist.append([])
        
        #ASSIGN INDEPENDENT VARIABLE AND EXTRA PARAMETERS FOR EACH MODEL TYPE
        k = 0
        fit[j].etc = []  
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ecl/tr':
                funcx.  append(fit[j].timeunit)
                funcxuc.append(fit[j].timeunituc)
                if hasattr(event[j].params, 'timeunit') and (event[j].params.timeunit != 'orbits'):
                    fit[j].etc.append([event[j].period])
                else:
                    fit[j].etc.append([1.0])
            elif functype[i] == 'ramp':
                funcx.  append(fit[j].timeunit)
                funcxuc.append(fit[j].timeunituc)
                fit[j].etc.append([])
            elif functype[i] == 'sinusoidal':
                funcx.  append(fit[j].timeunit)
                funcxuc.append(fit[j].timeunituc)
                fit[j].etc.append([])
            elif functype[i] == 'spline':
                funcx.  append([fit[j].timeunit, fit[j].flux])
                funcxuc.append([fit[j].timeunituc, fit[j].fluxuc])
                #fit[j].splineknots  = np.linspace(fit[j].timeunit[2], fit[j].timeunit[-3], event.sp.nknots)
                fit[j].etc.append([])#[fit[j].splineknots])
            #elif functype[i] == 'ippld':
            #    funcx.  append(fit[j].P_hat)
            #    funcxuc.append(fit[j].P_hatuc)
            #    fit[j].etc.append([])
            elif functype[i] == 'ippoly':
                if fit[j].model[k] == 'ipspline':
                    funcx.  append(fit[j].position)
                    funcxuc.append(fit[j].positionuc)
                    #CREATE KNOTS FOR IPSPLINE
                    fit[j].etc.append(np.meshgrid(
                                  np.linspace(fit[j].position[1].min(),fit[j].position[1].max(),numipxk), 
                                  np.linspace(fit[j].position[0].min(),fit[j].position[0].max(),numipyk)))
                elif fit[j].model[k] == 'posfluxlinip' or fit[j].model[k] == 'linip':
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
            elif functype[i] == 'ballardip':
                funcx.  append([fit[j].position, fit[j].ballardip, fit[j].flux])
                funcxuc.append([fit[j].positionuc, fit[j].ballardipuc, fit[j].fluxuc])
                fit[j].etc.append([])
            elif functype[i] == 'posoffset':
                fit[j].wherepos   = []
                fit[j].whereposuc = []
                for i in range(event[j].npos):
                    fit[j].wherepos  .append(np.where(fit[j].pos   == i)[0])
                    fit[j].whereposuc.append(np.where(fit[j].posuc == i)[0])
                funcx.  append([fit[j].nobj,   fit[j].wherepos])
                funcxuc.append([fit[j].nobjuc, fit[j].whereposuc])
                fit[j].etc.append([])
            elif functype[i] == 'vissen':
                funcx.  append([fit[j].frmvis,   event[j].params.vsknots])
                funcxuc.append([fit[j].frmvisuc, event[j].params.vsknots])
                fit[j].etc.append([])
            elif functype[i] == 'flatf':
                funcx.  append(fit[j].apdata)
                funcxuc.append(fit[j].apdatauc)
                fit[j].etc.append([])
            elif functype[i] == 'ipmap':
                funcx.  append(fit[j].posflux)
                #funcxuc.append([fit[j].y, fit[j].x, fit[j].fluxuc, fit[j].whereltraduc])   #medianip
                funcxuc.append(fit[j].posfluxuc)
                fit[j].etc.append([])
            elif functype[i] == 'gp':
                fit[j].gpfile = event[j].modeldir + "/d-" + event[j].eventname + "-gp.dat"
                if hasattr(event[j].params, 'expand') and event[j].params.expand > 1:
                    funcx.  append(fit[j].hirestimeuc)
                    funcxuc.append(fit[j].hirestimeuc)
                else:
                    funcx.  append(fit[j].timeunit)
                    funcxuc.append(fit[j].timeunituc)
                fit[j].etc.append([])
            elif functype[i] == 'noise':
                isnoise = True
                funcx.  append([])
                funcxuc.append([])
                if len(event[j].params.noisewavelet) == len(event[j].params.model):
                    fit[j].etc.append([event[j].params.noisewavelet[num]])
                else:
                    fit[j].etc.append([event[j].params.noisewavelet[0]])
            elif functype[i] == 'ortho':
                event[j].isortho = True
                # Create or load invserse transformation matrix and new coordinate origin
                if event[j].params.newortho == False:
                    for f in os.listdir('.'):
                        if f.endswith(event[j].eventname + '-ortho-' + fit[j].saveext + '.npz'):
                            print('Loading orthogonalization save file:', f, file=printout)
                            fit[j].orthofile = f
                            orthofile = np.load(f)
                            if len(orthofile['origin']) == diminvtrans:
                                fit[j].invtrans   = np.matrix(orthofile['invtrans'])
                                fit[j].trans      = np.matrix(orthofile['trans'])
                                fit[j].origin     = orthofile['origin']
                                fit[j].orthosigma = orthofile['sigma']
                            else:
                                print('***WARNING: Shape of saved inverse transformation matrix does not match number of elements in ortholist. Using identity matrix instead.', file=printout)
                # Load iortholist
                m = 0
                fit[j].iortholist = []
                for orthovar in ortholist:
                    ivar = getattr(fit[j].i, orthovar)
                    fit[j].iortholist.append(ivar)
                    iortholist[j].append(ivar + numparams[cummodels[j]])
                    #params[ivar + numparams[cummodels[j]]] -= fit[j].origin[m]
                    m += 1
                funcx.  append(fit[j].invtrans)
                funcxuc.append(fit[j].invtrans)
                fit[j].etc.append([fit[j].origin,fit[j].orthosigma])
                # Set orthogonal parameter names
                fit[j].opname   = np.copy(fit[j].parname)
                fit[j].opname[fit[j].iortholist] = opname[:diminvtrans]                  
            
            text   += fit[j].model[k] + " "
            k      += 1
    
    for j in range(numevents):
        fit[j].good      = event[j].good
        fit[j].isgood    = np.where(event[j].good)[0]
    for k in np.unique(nights):
        tonight   = np.where(nights == k)[0]
        # Construct white light curve from weighted spectroscopic lights curves
        # Fot Spitzer, this will default to ones.
        nasc.constructWhiteLC(event, fit, tonight, k)
        # Construct reference star spectroscopic light curves
        # Fot Spitzer, this will default to no reference star.
        nasc.refStarSpecLC2(event, fit, tonight)
    
    #RESTART FOR LOOP HERE
    flux      = []
    phase     = []
    bjdutc    = []
    bjdtdb    = []
    sigma     = []
    for j in range(numevents):
        flux.  append(fit[j].flux)
        phase. append(fit[j].phase)
        bjdutc.append(fit[j].bjdutc)
        bjdtdb.append(fit[j].bjdtdb)
        sigma. append(fit[j].sigma)
    
    #****LEAST-SQUARES FIT*****
    if leastsq == True:
        import scipy.optimize as op
        def modelfunc(freepars, params, sigma):
            #params[inonprior] = freepars
            params[ifreepars] = freepars
            for i in range(len(stepsize)):
                if stepsize[i] < 0:
                    params[i] = params[int(-stepsize[i]-1)]
            params[np.where(params < pmin)] = pmin[np.where(params < pmin)]
            params[np.where(params > pmax)] = pmax[np.where(params > pmax)]
            residuals = []
            #pedit     = np.copy(params)
            for j in range(numevents):
                fit0 = np.ones(fit[j].nobj)
                k    = 0
                for i in range(cummodels[j],cummodels[j+1]):
                    if   functype[i] == 'ortho':
                        pass
                    if   functype[i] == 'noise':
                        pass
                    elif (functype[i] == 'ipmap') or (functype[i] == 'spline'):
                        fit[j].etc[k] = fit0
                        fit0 *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
                    elif functype[i] == 'gp':
                        fit[j].etc[k] = [fit[j].flux/fit0, sigma[j],True,fit[j].gpfile,None]
                        gpmodel,params[iparams[i]] = myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
                        fit0 *= gpmodel
                    elif functype[i] == 'ballardip':
                        fit[j].ballardip = myfuncs[i](params[iparams[i]], funcx[i] , etc=[fit0])
                        fit0 *= fit[j].ballardip
                    else:
                        subparams = params[iparams[i]]
                        fit0 *= myfuncs[i](subparams, funcx[i], fit[j].etc[k])
                        params[iparams[i]] = subparams
                    k    += 1
                residuals = np.concatenate((residuals,(fit0 - flux[j])/sigma[j]),axis=0)
                #Apply priors, if they exist
                if len(fit[j].ipriors) > 0:
                    pbar   = fit[j].priorvals[:,0]  #prior mean
                    psigma = np.zeros(len(pbar))    #prior standard deviation
                    # Determine psigma based on which side of asymmetric Gaussian nextp is on
                    for i in range(len(fit[j].ipriors)):
                        if params[fit[j].ipriors[i]] < pbar[i]:
                            psigma[i] = fit[j].priorvals[i,1]
                        else:
                            psigma[i] = fit[j].priorvals[i,2]
                        #print(params[fit[j].ipriors[i]],pbar[i],psigma[i])
                        #print(sum(residuals**2),(np.sqrt(fit[j].nobj)*(params[fit[j].ipriors[i]] - pbar[i])/psigma[i])**2)
                        #plt.pause(1)
                        residuals = np.concatenate((residuals,[np.sqrt(fit[j].nobj)*(params[fit[j].ipriors[i]] - pbar[i])/psigma[i]]),axis=0)
            return residuals
        
        #Function for minimizer
        def fminfunc(freepars):
            params[inonprior] = freepars
            for i in range(len(stepsize)):
                if stepsize[i] < 0:
                    params[i] = params[int(-stepsize[i]-1)]
            #params[np.where(params < pmin)] = pmin[np.where(params < pmin)]
            #params[np.where(params > pmax)] = pmax[np.where(params > pmax)]
            logL = 0
            #pedit     = np.copy(params)
            for j in range(numevents):
                fit0 = np.ones(fit[j].nobj)
                k    = 0
                for i in range(cummodels[j],cummodels[j+1]):
                    if   functype[i] == 'ortho':
                        pass
                    elif (functype[i] == 'ipmap') or (functype[i] == 'spline'):
                        fit[j].etc[k] = fit0
                        fit0 *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
                    elif functype[i] == 'ballardip':
                        fit[j].ballardip = myfuncs[i](params[iparams[i]], funcx[i] , etc=[fit0])
                        fit0 *= fit[j].ballardip
                    elif hasattr(fit[j], 'timebins') and (functype[i] == 'ecl/tr' 
                                                      or  functype[i] == 'ramp' 
                                                      or  functype[i] == 'sinusoidal'):
                        # Average over high-resolution model
                        hiresmodel = myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
                        for n in range(len(fit[j].timebins)):
                            fit0[n] *= np.mean(hiresmodel[fit[j].timebins[n]])
                    elif functype[i] == 'noise':
                        etc       = fit[j].etc[k]
                        noisefunc = myfuncs[i]
                        noisepars = params[iparams[i]]
                    else:
                        subparams = params[iparams[i]]
                        fit0 *= myfuncs[i](subparams, funcx[i], fit[j].etc[k])
                        params[iparams[i]] = subparams
                    k    += 1
                logL += noisefunc(noisepars, fit0 - flux[j], etc)
                #residuals = np.concatenate((residuals,(fit0 - flux[j])),axis=0)
            return logL
        
        def constraint(x):
            return 1
        
        if isnoise == True:
            # Perform minimization using a truncated Newton algorithm
            print("Minimizing log-likelihood while considering correlated noise.")
            bounds = np.concatenate((pmin[:,np.newaxis][inonprior], pmax[:,np.newaxis][inonprior]), axis=1)
            #output, nfeval, rc = op.fmin_tnc(fminfunc, x0=params[inonprior], approx_grad=True, bounds=bounds, maxfun=max(500,50*len(params[inonprior])))#, scale=stepsize[inonprior].tolist())
            #print(output)
            output = op.fmin_cobyla(fminfunc, params[inonprior], constraint, rhobeg=stepsize[inonprior], rhoend=stepsize[inonprior]/1000.)
            params[inonprior] = output
        else:        
            print("Calculating least-squares fit.")
            output, cov_x, infodict, mesg, err = op.leastsq(modelfunc, params[ifreepars], args=(params, sigma), \
                 factor=100, ftol=1e-16, xtol=1e-16, gtol=1e-16, diag=1./stepsize[ifreepars], full_output=True)
            print(mesg)
            print(output)
        '''
        if (err >= 1) and (err <= 4):
            print("Fit converged without error.")
            print(output)
        else:
            print("WARNING: Fit did not converge. Using last iteration.")
        '''
        # update shared parameters
        '''
        for i in range(len(stepsize)):
            if stepsize[i] < 0:
                params[i] = params[-stepsize[i]-1]
        '''
    #**************************
    # Setup is finished.
    #**************************
    
    #CREATE MODEL USING INITIAL VALUES
    saveext = ''
    #if 'ipmap' in functype:
    #    plt.figure(600+num*numfigs)
    #    plt.clf()
    for j in range(numevents):
        saveext     = saveext + fit[j].saveext
        fit[j].fit0 = np.ones(fit[j].nobj)
        k           = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'ipmap':
                ipflux, fit[j].binipflux, fit[j].binipstd = myfuncs[i](params[iparams[i]], \
                                                     funcx[i], fit[j].fit0, retbinflux=True, retbinstd=True)
                fit[j].fit0     *= ipflux
                fit[j].binipflux = fit[j].binipflux.reshape(fit[j].xygrid[0].shape)
                fit[j].binipstd  = fit[j].binipstd. reshape(fit[j].xygrid[0].shape)
                #Define IP confidence level (or statistical significance) per bin
                #Confidence = signal / (noise / sqrt(sample size))
                #fit[j].ipconf    = np.sqrt(fit[j].numpts)*(fit[j].binipflux-fit[j].binipflux.mean())/fit[j].binipstd
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
                
                #PLOT CONFIDENCE vs. NUMBER OF POINTS PER BIN
                #ERROR BARS ARE MINIMUM AND MAXIMUM, NOT STANDARD DEVIATION
                '''
                plt.figure(600+num*numfigs)
                a = plt.suptitle(str(obj) + ' Mean Statistical Significance of Binned IP Flux', size=16)
                a = plt.errorbar(event[j].params.sssteps[1:numsteps+1], ipmean, [ipmean-ipmin,ipmax-ipmean], \
                                 fmt='o-', label=event[j].eventname)
                plt.xscale('log')
                plt.yscale('log')
                a = plt.xlabel('Number of Points per Bin', size=14)
                a = plt.ylabel('Statistical Significance', size=14)
                a = plt.legend(loc='best')
                #plt.errorbar(event[j].params.sssteps[1:numsteps+1],ipmean, ipstd, fmt='o-')
                #plt.plot(event[j].params.sssteps[1:numsteps+1], ipmean,'o-')
                #plt.plot(event[j].params.sssteps[1:numsteps+1], ipmin, '+')
                #plt.plot(event[j].params.sssteps[1:numsteps+1], ipmax, '+')
                plt.savefig(event[j].modeldir + "/" + str(obj) + "-fig" + str(num*numfigs+600) + "-" + saveext + ".png")
                '''
                #np.var(event[0].fit[0].binipstd[np.where(event[0].fit[0].numpts == 2)])
                #OPEN PROCESS ID TO SAVE allknots
                fit[j].allknotfile = event[j].modeldir + "/d-" + event[j].eventname + '-allknots-' + \
                                     fit[j].saveext + '.npy'
                fit[j].allknotpid  = open(fit[j].allknotfile, 'wb')
            elif functype[i] == 'spline':
                fit[j].etc[k] = fit[j].fit0
                fit[j].fit0 *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
            elif functype[i] == 'gp':
                fit[j].etc[k] = [fit[j].flux/fit[j].fit0, fit[j].sigma, False, None, None]
                fit[j].fit0  *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
            else:
                fit[j].fit0 *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
            k += 1
    
    # Transform parameters in ortholist to orthogonal parameters
    orthop = np.copy(params)
    for j in range(numevents):
        if event[j].isortho == True and event[j].params.newortho == False:
            orthop[iortholist[j]] = mc.orthoTrans(params[iortholist[j]], fit[j].trans, \
                                                 [fit[j].origin, fit[j].orthosigma])
    
    # Build starting parameters for multiple chains
    totnump = len(params)
    nchains = event[0].params.nchains[0]
    if nchains > 1:
        orthops = np.zeros((nchains,totnump)) + orthop
        for p in ifreepars:
            startpts = orthops[1:,p] + np.random.normal(0,stepsize[p]/ipfactor,nchains-1)
            '''
            outofrange = np.where(startpts < pmin[p])[0]
            if len(outofrange) > 0:
                startpts[outofrange] = pmin[p] + np.abs(np.random.normal(0,stepsize[p])/ipfactor/10.)
                #print(p,pmin[p],outofrange)
            outofrange = np.where(startpts > pmax[p])[0]
            if len(outofrange) > 0:
                startpts[outofrange] = pmax[p] - np.abs(np.random.normal(0,stepsize[p])/ipfactor/10.)
                #print(p,pmax[p],outofrange)
            orthops[1:,p] = startpts
            '''
            outofrange = np.where(startpts < pmin[p])[0]
            counter = 0
            while (len(outofrange) > 0) and (counter < 100):
                startpts[outofrange] = pmin[p] + np.abs(np.random.normal(0,stepsize[p]/ipfactor/10.,1))
                outofrange = np.where(startpts < pmin[p])[0]
                counter += 1
            if counter == 100:
                print("***WARNING: Parameter limits are too small.")
            outofrange = np.where(startpts > pmax[p])[0]
            counter = 0
            while (len(outofrange) > 0) and (counter < 100):
                startpts[outofrange] = pmax[p] - np.abs(np.random.normal(0,stepsize[p]/ipfactor/10.,1))
                outofrange = np.where(startpts < pmin[p])[0]
                counter += 1
            if counter == 100:
                print("***WARNING: Parameter limits are too small.")
            orthops[1:,p] = startpts
    else:
        orthops = orthop
    
    #RUN BURN IN OF MCMC MODEL FITTING ROUTINE
    print('MARK: ' + time.ctime() + ' : Start burn-in of mcmc model', file=printout)
    print('Number of iterations: ' + str(int(numit[0])), file=printout)
    if nchains == 1:
        allparams, bestop, numaccept, numit[0] = mcmc.mcmc9(flux, orthops, pmin, pmax, stepsize, numit[0], \
                         sigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, fit, nchains=0)
    elif isblock == True:
        allparams, bestop, numaccept, numit[0] = demc.demc_block(flux, orthops, pmin, pmax, stepsize, numit[0], \
                         sigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, fit, isGR=False, gamma=gamma)
    else:
        #allparams, bestop, numaccept, numit[0] = demc.demcz(flux, orthops, stepsize/ipfactor, pmin, pmax, stepsize, \
        #               numit[0], sigma, numparams, cummodels, functype, myfuncs, funcx, \
        #               iortholist, nights, fit, isGR=False, gamma=gamma)
        allparams, bestop, numaccept, numit[0] = demc.demc(flux, orthops, pmin, pmax, stepsize, numit[0], \
                 sigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, nights, fit, isGR=False, gamma=gamma)
    print('MARK: ' + time.ctime() + ' : End burn-in of mcmc model', file=printout)
    
    bestp = np.copy(bestop)
    for j in range(numevents):
        # Transform parameters in ortholist to original parameters
        if event[j].isortho == True and event[j].params.newortho == False:
            bestp[iortholist[j]] = mc.orthoInvTrans(bestop[iortholist[j]], fit[j].invtrans, \
                                                   [fit[j].origin, fit[j].orthosigma])
        # Record results for each fit
        fit[j].allparams  = allparams[numparams[cummodels[j]]:numparams[cummodels[j+1]]]
        fit[j].bestp      = bestp    [numparams[cummodels[j]]:numparams[cummodels[j+1]]]
        fit[j].bestfit   = np.ones(fit[j].nobj)
        fit[j].bestfituc = np.ones(fit[j].nobjuc)
        k = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'spline':
                fit[j].bestspl   = myfuncs[i](bestp[iparams[i]], funcx  [i], fit[j].bestfit)
                fit[j].bestspluc = myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].bestfituc)
                fit[j].bestfit  *= fit[j].bestspl
                fit[j].bestfituc*= fit[j].bestspluc
            elif functype[i] == 'ipmap':
                fit[j].bestmip,   fit[j].binipflux   = myfuncs[i](bestp[iparams[i]], funcx  [i], 
                                                                  fit[j].bestfit,   retbinflux=True)
                fit[j].bestmipuc, fit[j].binipfluxuc = myfuncs[i](bestp[iparams[i]], funcxuc[i], 
                                                                  fit[j].bestfituc, retbinflux=True)
                fit[j].bestmipuc[fit[j].isclipmask] = fit[j].bestmip
                fit[j].bestfit  *= fit[j].bestmip
                fit[j].bestfituc*= fit[j].bestmipuc
                fit[j].allknotpid.close()
                fit[j].allknotpid  = open(fit[j].allknotfile, 'w')
            elif functype[i] == 'gp':
                fit[j].etc[k]    = [fit[j].flux/fit[j].bestfit, fit[j].sigma, False, None, None]
                fit[j].bestgp    = myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
                fit[j].bestgpuc  = myfuncs[i](params[iparams[i]], funcxuc[i], [fit[j].fluxuc/fit[j].bestfituc, fit[j].sigmauc])
                fit[j].bestgpuc[fit[j].isclipmask] = fit[j].bestgp
                fit[j].bestfit  *= fit[j].bestgp
                fit[j].bestfituc*= fit[j].bestgpuc
            elif functype[i] == 'ballardip':
                print("Re-computing intrapixel effect")
                fit[j].ballardip   = myfuncs[i](bestp[iparams[i]], funcx[i],   etc=[fit[j].bestfit  ])
                fit[j].ballardipuc = myfuncs[i](bestp[iparams[i]], funcxuc[i], etc=[fit[j].bestfituc])
                fit[j].ballardipuc[fit[j].isclipmask] = fit[j].ballardip
                fit[j].bestfit    *= fit[j].ballardip
                fit[j].bestfituc  *= fit[j].ballardipuc
            else:
                fit[j].bestfit    *= myfuncs[i](bestp[iparams[i]], funcx  [i], fit[j].etc[k])
                fit[j].bestfituc  *= myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].etc[k])
            k += 1
    
    for j in range(numevents):
        fit[j].redchisq = sum((fit[j].bestfit - fit[j].flux)**2 / fit[j].sigma**2) / \
                              (fit[j].nobj - fit[j].numfreepars)
        print("Reduced Chi^2: " + str(fit[j].redchisq), file=printout)
        fit[j].arate    = 100.0 * numaccept / numit[0]
    
    print("Acceptance rate = " + str(fit[j].arate) + "%")
    print("Acceptance rate = " + str(fit[j].arate) + "%", file=printout)
    print("Best Parameters:", file=printout)
    print(bestp, file=printout)
    
    for j in range(numevents):
        fignum   = 6000+num*numfigs+j*100
        savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
        plots.trace(event[j], fit[j], fignum, savefile)
        if not isinteractive:
            plt.close(fignum)
        else:
            plt.pause(0.01)
    
    #BIN DATA
    for j in range(numevents):
        fit[j].binfluxuc  = np.zeros(nbins[j])
        fit[j].binphase   = np.zeros(nbins[j])
        fit[j].binphaseuc = np.zeros(nbins[j])
        fit[j].binbjdutc  = np.zeros(nbins[j])
        fit[j].binbjdutcuc= np.zeros(nbins[j])
        fit[j].binbjdtdb  = np.zeros(nbins[j])
        fit[j].binbjdtdbuc= np.zeros(nbins[j])
        fit[j].binstduc   = np.zeros(nbins[j])
        fit[j].binfit0    = np.zeros(nbins[j])
        fit[j].binbestfit = np.zeros(nbins[j])
        for i in range(nbins[j]):
            startuc              = int(1.*i*fit[j].nobjuc/nbins[j])
            enduc                = int(1.*(i+1)*fit[j].nobjuc/nbins[j])
            fit[j].binphaseuc[i] = np.median(fit[j].phaseuc [startuc:enduc])
            fit[j].binbjdutcuc[i]= np.median(fit[j].bjdutcuc[startuc:enduc])
            fit[j].binbjdtdbuc[i]= np.median(fit[j].bjdtdbuc[startuc:enduc])
            fit[j].binfluxuc[i]  = sum(fit[j].fluxuc[startuc:enduc]/fit[j].sigmauc[startuc:enduc]**2)/ \
                			       sum(1/fit[j].sigmauc[startuc:enduc]**2)
            fit[j].binstduc[i]   = np.sqrt(1 / sum(1/fit[j].sigmauc[startuc:enduc]**2))
            start                = int(1.*i*fit[j].nobj/nbins[j])
            end                  = int(1.*(i+1)*fit[j].nobj/nbins[j])
            fit[j].binphase[i]   = np.median(fit[j].phase [start:end])
            fit[j].binbjdutc[i]  = np.median(fit[j].bjdutc[start:end])
            fit[j].binbjdtdb[i]  = np.median(fit[j].bjdtdb[start:end])
            fit[j].binfit0[i]    = np.mean(fit[j].fit0[start:end])
            fit[j].binbestfit[i] = np.mean(fit[j].bestfit[start:end])
    
    #MODIFY sigma SUCH THAT REDUCED CHI-SQUARED = 1
    #SPITZER OVER-ESTIMATES ERRORS
    newsigma = []
    for j in range(numevents):
        if chi2flag == 1:
            fit[j].newsigma   = fit[j].sigma   * np.sqrt(fit[j].redchisq)
            fit[j].newsigmauc = fit[j].sigmauc * np.sqrt(fit[j].redchisq)
        else:
            fit[j].newsigma   = np.copy(fit[j].sigma)
            fit[j].newsigmauc = np.copy(fit[j].sigmauc)
        newsigma.append(fit[j].newsigma)
    
    #NEW LEAST-SQUARES FIT USING MODIFIED sigma VALUES
    if leastsq == True:
        if isnoise == True:
            # Perform minimization using a truncated Newton algorithm
            '''
            http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_tnc.html#scipy.optimize.fmin_tnc
            Wright S., Nocedal J. (2006), 'Numerical Optimization'
            Nocedal, J, and S J Wright. 2006. Numerical Optimization. Springer New York.
            '''
            print("Minimizing log-likelihood while considering correlated noise and constraints.")
            bounds = np.concatenate((pmin[:,np.newaxis][inonprior], pmax[:,np.newaxis][inonprior]), axis=1)
            #output, nfeval, rc = op.fmin_tnc(fminfunc, x0=bestp[inonprior], approx_grad=True, bounds=bounds, maxfun=max(500,50*len(params[inonprior])))#, scale=stepsize[inonprior].tolist())
            #print(output)
            output = op.fmin_cobyla(fminfunc, params[inonprior], constraint, rhobeg=stepsize[inonprior], rhoend=stepsize[inonprior]/1000.)
            bestp[inonprior] = output
        else:
            print("Re-calculating least-squares fit with new errors.")
            # ccampo added 8-24-2010 to allow new errors
            # ccampo changed params to bestp
            output, err = op.leastsq(modelfunc, bestp[ifreepars], args=(bestp, newsigma), factor=100, \
                                         ftol=1e-16, xtol=1e-16, gtol=1e-16, diag=1./stepsize[ifreepars])
            if (err >= 1) and (err <= 4):
                print("Fit converged without error.")
            else:
                print("WARNING: Error with least squares fit!")
            print("Least squares fit best parameters:")
            print(output)
        # ccampo: update best parameters w/ shared values
        for i in range(len(stepsize)):
            if stepsize[i] < 0:
                bestp[i] = bestp[int(-stepsize[i]-1)]

    #FIX IP MAPPING TO bestmip IF isfixipmap = TRUE
    for j in range(numevents):
        if hasattr(event[j].params, 'isfixipmap') and event[j].params.isfixipmap == True:
            for i in range(cummodels[j],cummodels[j+1]):
                if functype[i] == 'ipmap':
                    foo = funcx  .pop(i)
                    foo = funcxuc.pop(i)
                    funcx  .append([fit[j].bestmip,   fit[j].binipflux,   np.zeros(len(fit[j].wbfipmask  ))])
                    funcxuc.append([fit[j].bestmipuc, fit[j].binipfluxuc, np.zeros(len(fit[j].wbfipmaskuc))])
                    myfuncs[i] = mc.fixipmapping
    
    #Set isoptimize to False in Gaussian Process
    for j in range(numevents):
        k = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'gp':
                fit[j].etc[k][2:5] = [False, None, None]
            k += 1
    
    # Transfor
    # Transform parameters in ortholist to orthogonal parameters
    bestop = np.copy(bestp)
    for j in range(numevents):
        if event[j].isortho == True and event[j].params.newortho == False:
            bestop[iortholist[j]] = mc.orthoTrans(bestp[iortholist[j]], fit[j].trans, \
                                                 [fit[j].origin, fit[j].orthosigma])
    
    # Build starting parameters for multiple chains
    nchains = event[0].params.nchains[1]
    if nchains > 1:
        bestops     = np.zeros((nchains,totnump)) + bestop
        stdpburnin  = np.std(allparams,axis=1)
        for p in ifreepars:
            startpts = bestops[1:,p] + np.random.normal(0,stdpburnin[p],nchains-1)
            '''
            outofrange = np.where(startpts < pmin[p])[0]
            if len(outofrange) > 0:
                startpts[outofrange] = pmin[p] + np.abs(np.random.normal(0,stdpburnin[p]/10.,1))
            outofrange = np.where(startpts > pmax[p])[0]
            if len(outofrange) > 0:
                startpts[outofrange] = pmax[p] - np.abs(np.random.normal(0,stdpburnin[p]/10.,1))
            bestops[1:,p] = startpts
            '''
            outofrange = np.where(startpts < pmin[p])[0]
            counter = 0
            while (len(outofrange) > 0) and (counter < 100):
                startpts[outofrange] = pmin[p] + np.abs(np.random.normal(0,np.std(allparams[p])/10.,1))
                outofrange = np.where(startpts < pmin[p])[0]
                counter += 1
            if counter == 100:
                print("***WARNING: Parameter limits are too small.")
            outofrange = np.where(startpts > pmax[p])[0]
            counter = 0
            while (len(outofrange) > 0) and (counter < 100):
                startpts[outofrange] = pmax[p] - np.abs(np.random.normal(0,np.std(allparams[p])/10.,1))
                outofrange = np.where(startpts < pmin[p])[0]
                counter += 1
            if counter == 100:
                print("***WARNING: Parameter limits are too small.")
            bestops[1:,p] = startpts
            #bestops[1:,p] += np.random.normal(0,stepsize[p]/1.,nchains-1)
    else:
        bestops = bestop
    
    #RUN MCMC FITTING ROUTINE FOR ANALYSIS
    print('MARK: ' + time.ctime() + ' : Start full mcmc model', file=printout)
    print('Max number of iterations: ' + str(int(numit[1])), file=printout)
    if nchains == 1 and isresidperm == False:
        allorthop, bestop, numaccept, numit[1] = mcmc.mcmc9(flux, bestops, pmin, pmax, stepsize, numit[1], \
                       newsigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, fit, nchains=4)
    elif isblock == True:
        allorthop, bestop, numaccept, numit[1] = demc.demc_block(flux, bestops, pmin, pmax, stepsize, numit[1], \
                       newsigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, fit, isGR=True, gamma=gamma)
    elif isresidperm == True:
        print("Caution: Priors become fixed with residual permutation.")
        allorthop = residperm.permute2(event, fit, bestop, pmin, pmax, flux, cummodels, functype, myfuncs, funcx, funcxuc, iparams, stepsize, inonprior, nights, numit[1])
        numaccept = numit[1]
    else:
        allorthop, bestop, numaccept, numit[1] = demc.demcz(flux, bestops, stdpburnin, pmin, pmax, stepsize, \
                       numit[1], newsigma, numparams, cummodels, functype, myfuncs, funcx, \
                       iortholist, nights, fit, isGR=True, gamma=gamma)
        #allorthop, bestop, numaccept, numit[1] = demc.demc(flux, bestops, pmin, pmax, stepsize, numit[1], \
        #         newsigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, nights, fit, isGR=True, gamma=gamma)
    print('MARK: ' + time.ctime() + ' : End full mcmc model', file=printout)
    print('Completed ' + str(int(numit[1])) + ' iterations.', file=printout)
    
    bestp     = np.copy(bestop)
    allparams = np.copy(allorthop)
    for j in range(numevents):
        # Record trace of decorrelated parameters, if any, otherwise same as fit[j].allparams
        fit[j].allorthop  = np.copy(allorthop[numparams[cummodels[j]]:numparams[cummodels[j+1]]])
        fit[j].bestop     = np.copy(bestp    [numparams[cummodels[j]]:numparams[cummodels[j+1]]])
        # Transform parameters in ortholist to original parameters
        if event[j].isortho == True and event[j].params.newortho == False:
            bestp[iortholist[j]]     = mc.orthoInvTrans(bestop[iortholist[j]], fit[j].invtrans, \
                                                   [fit[j].origin, fit[j].orthosigma])
            allparams[iortholist[j]] = mc.orthoInvTrans(allorthop[iortholist[j]], fit[j].invtrans, \
                                                   [fit[j].origin, fit[j].orthosigma])
    
    # ccampo: checks best parameters against minimizer output to make sure they are the same.
    if leastsq == True:
        if np.all(np.abs(output - bestp[ifreepars]) <= 1e-8) == True:
            print("Least-squares minimizer values are the best fitting paramters!", file=printout)
        else:
            print("WARNING: minimizer values DO NOT MATCH best MCMC values!!!", file=printout)
            print("MINIMIZER VALUES:", file=printout)
            print(output, file=printout)
            print("MCMC VALUES:", file=printout)
            print(bestp[ifreepars], file=printout)
            print("DIFFERENCES:", file=printout)
            print(bestp[ifreepars] - output, file=printout)
            
            if isnoise == True:
                # Perform minimization using a truncated Newton algorithm
                print("Doing a last minimizing with final MCMC values.")
                bounds = np.concatenate((pmin[:,np.newaxis][inonprior], pmax[:,np.newaxis][inonprior]), axis=1)
                #output, nfeval, rc = op.fmin_tnc(fminfunc, x0=params[inonprior], approx_grad=True, bounds=bounds, maxfun=max(500,50*len(params[inonprior])))#, scale=stepsize[inonprior].tolist())
                output = op.fmin_cobyla(fminfunc, params[inonprior], constraint, rhobeg=stepsize[inonprior], rhoend=stepsize[inonprior]/1000.)
                print("Minimizer's best parameters:")
                print(output)
                bestp[inonprior] = output
            else:
                print("Doing a last least-squares fit with final MCMC values...")                 
                # ccampo added 8-24-2010 to allow new errors          
                # ccampo changed params to bestp
                output, err = op.leastsq(modelfunc, bestp[ifreepars], args=(bestp, newsigma), factor=100, \
                                         ftol=1e-16, xtol=1e-16, gtol=1e-16, diag=1./stepsize[ifreepars])
                if (err >= 1) and (err <= 4):
                    print("Fit converged without error.")
                else:
                    print("WARNING: Error with least squares fit!")
                print("Least squares fit best parameters:")
                print(output)
            # ccampo: update best parameters w/ shared values
            for i in range(len(stepsize)):
                if stepsize[i] < 0:
                    bestp[i] = bestp[int(-stepsize[i]-1)]
        # ccampo: adding fit.minimizerp variable
        for j in range(numevents):
            fit[j].minimizerp = output
    
    #RECORD BEST MODELS FOR EACH MODEL TYPE
    for j in range(numevents):
        fit[j].allparams  = allparams[numparams[cummodels[j]]:numparams[cummodels[j+1]]]
        fit[j].meanknots  = 0
        fit[j].medianknots= 0
        fit[j].stdknots   = 0
        fit[j].bestp      = bestp    [numparams[cummodels[j]]:numparams[cummodels[j+1]]]
        fit[j].bestfit    = np.ones(fit[j].nobj)
        k = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'spline':
                fit[j].etc[k]   = fit[j].bestfit
                fit[j].bestfit *= myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
            elif functype[i] == 'ipmap':
                fit[j].allknotpid.close()
                del(fit[j].allknotpid)
                fit[j].etc[k]   = fit[j].bestfit
                fit[j].bestfit *= myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
            elif functype[i] == 'gp':
                fit[j].bestfitgp    = np.copy(fit[j].bestfit)
                fit[j].etc[k]       = [fit[j].flux/fit[j].bestfit, fit[j].newsigma, False, None, None]
                #foo = myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
                fit[j].bestgp       = myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
                fit[j].bestfit     *= fit[j].bestgp
            else:
                fit[j].bestfit *= myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
            #UPDATE BEST-FIT PARAMETER LIST FOR NEXT RUN
            if fit[j].model[k] != 'ipspline':
                parlist[j][k][2][0] = bestp[iparams[i]]
            k += 1
        
        #COMPUTE SEPARATE BEST ECLIPSE, RAMP & INTRA-PIXEL MODELS AND BEST PARAMETERS
        fit[j].bestecl      = np.ones(fit[j].nobj)
        fit[j].bestramp     = np.ones(fit[j].nobj)
        fit[j].bestsin      = np.ones(fit[j].nobj)
        fit[j].bestip       = np.ones(fit[j].nobj)
        fit[j].bestpos      = np.ones(fit[j].nobj)
        fit[j].bestvs       = np.ones(fit[j].nobj)
        fit[j].bestff       = np.ones(fit[j].nobj)
        fit[j].bestmip      = np.ones(fit[j].nobj)
        fit[j].bestspl      = np.ones(fit[j].nobj)
        fit[j].bestpip      = []
        fit[j].bestfitmip   = []
        fit[j].bestecluc    = np.ones(fit[j].nobjuc)
        fit[j].bestrampuc   = np.ones(fit[j].nobjuc)
        fit[j].bestsinuc    = np.ones(fit[j].nobjuc)
        fit[j].bestipuc     = np.ones(fit[j].nobjuc)
        fit[j].bestposuc    = np.ones(fit[j].nobjuc)
        fit[j].bestvsuc     = np.ones(fit[j].nobjuc)
        fit[j].bestffuc     = np.ones(fit[j].nobjuc)
        fit[j].bestmipuc    = np.ones(fit[j].nobjuc)
        fit[j].bestspluc    = np.ones(fit[j].nobjuc)
        fit[j].bestgpuc     = np.ones(fit[j].nobjuc)
        fit[j].bestfitmipuc = []
        #COMPUTE SYSTEM FLUX
        fit[j].systemflux       = 1
        if hasattr(fit[j].i, 'flux'):
            fit[j].systemflux  *= fit[j].bestp[fit[j].i.flux]
        if hasattr(fit[j].i, 'flux2'):
            fit[j].systemflux  *= fit[j].bestp[fit[j].i.flux2]
        if hasattr(fit[j].i, 'flux3'):
            fit[j].systemflux  *= fit[j].bestp[fit[j].i.flux3]
        if hasattr(fit[j].i, 'trflux'):
            fit[j].systemflux  *= fit[j].bestp[fit[j].i.trflux]
        if hasattr(fit[j].i, 'trflux2'):
            fit[j].systemflux  *= fit[j].bestp[fit[j].i.trflux2]
        if hasattr(fit[j].i, 'trspf'):
            fit[j].systemflux  *= fit[j].bestp[fit[j].i.trspf]
        if hasattr(fit[j].i, 'trspf2'):
            fit[j].systemflux  *= fit[j].bestp[fit[j].i.trspf2]
        #fit[j].binipflux    = np.ones(fit[j].xygrid[0].shape)
        k = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ecl/tr':
                fit[j].bestecl     *= myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
                fit[j].bestecluc   *= myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].etc[k])
            elif functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'ramp':
                fit[j].bestramp    *= myfuncs[i](bestp[iparams[i]], funcx[i])
                fit[j].bestrampuc  *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
            elif functype[i] == 'sinusoidal':
                fit[j].bestsin     *= myfuncs[i](bestp[iparams[i]], funcx[i])
                fit[j].bestsinuc   *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
            elif functype[i] == 'ippoly':
                fit[j].bestip      *= myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].etc[k])
                fit[j].bestipuc    *= myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].etc[k])
                fit[j].bestpip      = bestp[iparams[i]]
            elif functype[i] == 'posoffset':
                fit[j].bestpos     *= myfuncs[i](bestp[iparams[i]], funcx[i],   fit[j].etc[k])
                fit[j].bestposuc   *= myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].etc[k])
            elif functype[i] == 'vissen':
                fit[j].bestvs      *= myfuncs[i](bestp[iparams[i]], funcx[i])
                fit[j].bestvsuc    *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
                fit[j].bestpvs      = bestp[iparams[i]]
            elif functype[i] == 'flatf':
                fit[j].bestff      *= myfuncs[i](bestp[iparams[i]], funcx[i])
                fit[j].bestffuc    *= myfuncs[i](bestp[iparams[i]], funcxuc[i])
            elif functype[i] == 'spline':
                fit[j].bestfitspl   = fit[j].bestecl   * fit[j].bestramp   * fit[j].bestip   * \
                                      fit[j].bestpos   * fit[j].bestvs     * fit[j].bestff   * \
                                      fit[j].bestsin
                fit[j].bestfitspluc = fit[j].bestecluc * fit[j].bestrampuc * fit[j].bestipuc * \
                                      fit[j].bestposuc * fit[j].bestvsuc   * fit[j].bestffuc * \
                                      fit[j].bestsinuc
                fit[j].bestspl     *= myfuncs[i](bestp[iparams[i]], funcx[i], fit[j].bestfitspl)
                fit[j].bestspluc   *= myfuncs[i](bestp[iparams[i]], funcxuc[i], fit[j].bestfitspluc)
            elif functype[i] == 'ipmap':
                fit[j].bestfitmip   = fit[j].bestecl   * fit[j].bestramp   * fit[j].bestip   * \
                                      fit[j].bestpos   * fit[j].bestvs     * fit[j].bestff   * \
                                      fit[j].bestsin
                fit[j].bestfitmipuc = fit[j].bestecluc * fit[j].bestrampuc * fit[j].bestipuc * \
                                      fit[j].bestposuc * fit[j].bestvsuc   * fit[j].bestffuc * \
                                      fit[j].bestsinuc
                fit[j].bestmip, fit[j].binipflux = myfuncs[i](bestp[iparams[i]], funcx[i], \
                                                 fit[j].bestfitmip, retbinflux=True)
                fit[j].bestmipuc ,fit[j].binipfluxuc= myfuncs[i](bestp[iparams[i]], funcxuc[i], \
                                                 fit[j].bestfitmipuc, retbinflux=True)
                #fit[j].bestmip     *= a
                #fit[j].bestmipuc   *= b
                fit[j].bestmipuc *= np.mean(fit[j].bestmip)/np.mean(fit[j].bestmipuc[fit[j].isclipmask])
                fit[j].bestmipuc[fit[j].isclipmask] = fit[j].bestmip
                fit[j].binipflux    = fit[j].binipflux.  reshape(fit[j].xygrid[0].shape)
                fit[j].binipfluxuc  = fit[j].binipfluxuc.reshape(fit[j].xygrid[0].shape)
                fit[j].bestpmip     = bestp[iparams[i]]
                #fit[j].systemflux   = np.mean(fit[j].binipflux.flatten()[np.where(fit[j].binfluxmask)])
            elif functype[i] == 'gp':
                #fit[j].bestfitgp    = fit[j].bestecl   * fit[j].bestramp   * fit[j].bestip   * \
                #                      fit[j].bestpos   * fit[j].bestvs     * fit[j].bestff   * \
                #                      fit[j].bestsin
                fit[j].bestfitgpuc  = fit[j].bestecluc * fit[j].bestrampuc * fit[j].bestipuc * \
                                      fit[j].bestposuc * fit[j].bestvsuc   * fit[j].bestffuc * \
                                      fit[j].bestsinuc
                #fit[j].bestgp      *= myfuncs[i](bestp[iparams[i]], funcx[i], [fit[j].flux/fit[j].bestfitgp, fit[j].newsigma])
                fit[j].bestgpuc    *= myfuncs[i](bestp[iparams[i]], funcxuc[i], [fit[j].fluxuc/fit[j].bestfitgpuc, fit[j].newsigmauc, False, None, None])
                fit[j].bestpgp      = bestp[iparams[i]]
            k += 1
        
        #UPDATE INITIAL PARAMETERS
        if event[j].params.modelfile == None:
            modelfile = event[j].initvalsfile
        elif event[j].params.modelfile.__class__ == str:
            modelfile = event[j].ancildir + event[j].params.modelfile
        elif len(event[j].params.modelfile) == len(event[j].params.model):
            modelfile = event[j].ancildir + event[j].params.modelfile[num] 
        else:
            modelfile = event[j].ancildir + event[j].params.modelfile[0]
        pe.write(modelfile, parlist[j])
        ier = os.system("cp %s %s/." % (modelfile, event[j].modeldir))
        
        # Calculate acceptance rate
        fit[j].arate      = 100.0 * numaccept / numit[1]
        print("Acceptance rate = " + str(fit[j].arate) + "%")
        print("Acceptance rate = " + str(fit[j].arate) + "%", file=printout)
    
    print("Calculating autocorrelation of free parameters.")
    maxstepsize = 1
    for j in range(numevents):
        acsize = np.int(np.min((numit[1],2e3)))
        fit[j].autocorr = np.zeros((len(fit[j].nonfixedpars), acsize))
        fit[j].ess      = np.zeros(len(fit[j].nonfixedpars)) + numit[1]
        k = 0
        for i in fit[j].nonfixedpars:
            # Calculate autocorrelation
            meanapi  = np.mean(fit[j].allparams[i])
            autocorr = np.correlate(fit[j].allparams[i,:acsize]-meanapi, 
                                    fit[j].allparams[i,:acsize]-meanapi, mode='full')
            fit[j].autocorr[k] = (autocorr[int(autocorr.size/2):] / autocorr.max())
            # Don't calculate effective sample size for fixed priors during residual permutation calculation
            if np.std(fit[j].allparams[i]) > 1e-12:
                try:
                    # Calculate effective sample size (ESS, see Kass et al.; 1998)
                    cutoff = np.where(fit[j].autocorr[k] < 0.01)[0][0]  #First instance where autocorr < 0.01
                    fit[j].ess[k] = numit[1] / (1 + 2*np.sum(fit[j].autocorr[k,:cutoff]))
                except:
                    print('Could not calculate ESS; autocorrelation never reaches < 0.01.')
                    fit[j].ess[k] = 2
            k += 1
        miness = np.min(fit[j].ess)
        print(event[j].eventname, "effective sample size:", np.round(miness), file=printout)
        fit[j].stepsize          = np.round(numit[1]/miness)
        if hasattr(event[j].params, 'stepsize') == False or event[j].params.stepsize == None:
            event[j].params.stepsize = fit[j].stepsize
            print(event[j].eventname, "new stepsize:", event[j].params.stepsize, file=printout)
        maxstepsize = int(np.max((maxstepsize, event[j].params.stepsize)))
    
    # Compute transformation matrix and plot orthogonal parameter correlations
    for j in range(numevents):
        if event[j].params.newortho == True:
            print("Computing transformation matrix for " +event[j].eventname+ ".")
            steps = int(np.ceil(numit[0]/2000.))
            #Orthogonalization using Principal Component Analysis (PCA) class
            orthodata = fit[j].allparams[fit[j].iortholist,::steps]
            fit[j].ortho = plt.mlab.PCA(orthodata.T)      #ortho object
            fit[j].origin = fit[j].ortho.mu             #means of orthodata
            # Compute inverse transformation matrix
            # Y = trans * (data-mu).T
            # (data-mu).T = invtrans * Y
            # Note: np.matrix type, not np.array
            fit[j].trans      = np.matrix(fit[j].ortho.Wt)
            fit[j].invtrans   = np.matrix(np.linalg.inv(fit[j].trans))
            fit[j].orthosigma = fit[j].ortho.sigma
            print(event[j].eventname + " inverse transformation matrix:",file=printout)
            print(fit[j].invtrans,file=printout)
            
            # Update best-fit orthogonal parameters
            #fit[j].origin = np.copy(neworigin)
            fit[j].bestop[fit[j].iortholist] = mc.orthoTrans(fit[j].bestp[fit[j].iortholist], fit[j].trans, \
                                                            [fit[j].origin, fit[j].orthosigma])
            bestop[iortholist[j]] = mc.orthoTrans(bestp[iortholist[j]], fit[j].trans, \
                                                            [fit[j].origin, fit[j].orthosigma])
            print(event[j].eventname + " best parameters after orthogonalization:",file=printout)
            print(fit[j].bestop,file=printout)
            
            # Correlation plot of original parameters
            fignum   = 6000+num*numfigs+j*100
            savefile = event[j].modeldir +"/"+ event[j].eventname +"-fig"+ str(fignum) +"-"+ fit[j].saveext +".png"
            plots.hist2d(event[j], fit[j], fignum, savefile, allparams=orthodata, iparams=fit[j].iortholist)
            
            # Correlation plot of transformed parameters
            fignum   = 6010+num*numfigs+j*100
            savefile = event[j].modeldir +"/"+ event[j].eventname +"-fig"+ str(fignum) +"-"+ fit[j].saveext +".png"
            plots.hist2d(event[j], fit[j], fignum, savefile, allparams=fit[j].ortho.Y.T, parname=fit[j].opname, iparams=fit[j].iortholist)
            '''
            # Update arrays for MCMC
            k = 0
            for i in range(cummodels[j],cummodels[j+1]):
                if functype[i] == 'ortho':
                    funcx[k]      = fit[j].invtrans
                    funcxuc[k]    = fit[j].invtrans
                    fit[j].etc[k] = fit[j].origin
                k += 1
            '''
    
    #Calculate mean, median and standard deviation of output parameters
    totnump           = len(params)
    medianp           = np.zeros((totnump, 2))
    meanp             = np.zeros((totnump, 2))
    medianp[:, 0]     = np.median(allparams[:,::maxstepsize],axis=1)
    medianp[:, 1]     = np.std   (allparams[:,::maxstepsize],axis=1)
    meanp[:, 0]       = np.mean  (allparams[:,::maxstepsize],axis=1)
    meanp[:, 1]       = np.copy  (medianp[:, 1])
    
    #NEED TO RESTART FOR LOOP HERE!
    rampp             = np.copy(bestp)
    for j in range(numevents):
        fit[j].medianp    = medianp[numparams[cummodels[j]]:numparams[cummodels[j+1]]]
        fit[j].meanp      = meanp  [numparams[cummodels[j]]:numparams[cummodels[j+1]]]
        fit[j].medianfit  = np.ones(fit[j].nobj)
        fit[j].meanfit    = np.ones(fit[j].nobj)
        k = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'spline':
                fit[j].etc[k] = fit[j].bestfitspl
                fit[j].medianfit *= myfuncs[i](medianp[iparams[i],0], funcx[i], fit[j].etc[k])
                fit[j].meanfit   *= myfuncs[i](meanp  [iparams[i],0], funcx[i], fit[j].etc[k])
            elif functype[i] == 'ipmap':
                fit[j].etc[k] = fit[j].bestfitmip
                fit[j].medianfit *= myfuncs[i](medianp[iparams[i],0], funcx[i], fit[j].etc[k])
                fit[j].meanfit   *= myfuncs[i](meanp  [iparams[i],0], funcx[i], fit[j].etc[k])
            elif functype[i] == 'gp':
                fit[j].etc[k] = [fit[j].flux/fit[j].bestfitgp, fit[j].newsigma, False, None, None]
                fit[j].medianfit *= myfuncs[i](medianp[iparams[i],0], funcx[i], fit[j].etc[k])
                fit[j].meanfit   *= myfuncs[i](meanp  [iparams[i],0], funcx[i], fit[j].etc[k])
            else:
                fit[j].medianfit *= myfuncs[i](medianp[iparams[i],0], funcx[i], fit[j].etc[k])
                fit[j].meanfit   *= myfuncs[i](meanp  [iparams[i],0], funcx[i], fit[j].etc[k])
            k += 1
        
        #COMPUTE RESIDUALS AND REDUCED CHI^2
        fit[j].residuals   =      fit[j].flux - fit[j].bestfit
        fit[j].redchisq    = sum((fit[j].residuals)**2 / fit[j].newsigma**2) /   \
                                 (fit[j].nobj - fit[j].numfreepars)
        fit[j].oldredchisq = sum((fit[j].residuals)**2 / fit[j].sigma**2) /      \
                                 (fit[j].nobj - fit[j].numfreepars)
        print("Median Parameters With Standard Deviations:", file=printout)
        print(fit[j].medianp, file=printout)
        print("Best Parameters:", file=printout)
        print(fit[j].bestp, file=printout)
        print("Reduced Chi^2: " + str(fit[j].redchisq), file=printout)
        
        #COMPUTE BAYESIAN INFORMATION CRITERION (BIC)
        #fit[j].bic = fit[j].nobj*np.log(sum(fit[j].residuals**2)/fit[j].nobj) + \
        #             fit[j].numfreepars*np.log(fit[j].nobj) #OLD VERSION - DO NOT USE
        
        #Residuals are nomalized by dividing by the variance of the residuals from:
        #the previous model fit from the same channel, or
        #from itself if this is the first model of this channel.
        #This results in the same normalizing factor for all models of the same channel.
        #bicloc = j
        #for i in range(0,j):
        #    if event[j].photchan == event[i].photchan:
        #        bicloc = i
        #fit[j].bic = sum(fit[j].residuals**2)/(np.std(event[bicloc].fit[0].residuals))**2 + \
        #                 fit[j].numfreepars*np.log(fit[j].nobj) #+ \
                         #np.sum(np.log(fit[j].numpts[np.where(fit[j].numpts > 0)]))
        try:
            fit[j].aic = sum((fit[j].residuals/event[j].fit[0].newsigma)**2) + 2*fit[j].numfreepars
            fit[j].bic = sum((fit[j].residuals/event[j].fit[0].newsigma)**2) +   fit[j].numfreepars*np.log(fit[j].nobj)
            print("AIC = " + str(fit[j].aic), file=printout)
            print("BIC = " + str(fit[j].bic), file=printout)
        except:
            fit[j].aic = 0
            fit[j].bic = 0
            print("Error computing AIC/BIC. Likely cause is a shape mismatch due to different number of knots.")
        
        #Create unclipped ramp with eclipse
        fit[j].rampuc2    = np.ones(fit[j].nobjuc)
        k = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'ipmap':
                fit[j].rampuc2 *= fit[j].bestmipuc
            elif functype[i] == 'spline':
                fit[j].rampuc2 *= fit[j].bestspluc
            else:
                fit[j].rampuc2 *= myfuncs[i](rampp[iparams[i]], funcxuc[i], fit[j].etc[k])
            k += 1
        
    #NEED TO RESTART FOR LOOP HERE!
    for j in range(numevents):
        #Create data without the eclipse
        #inotecl      = np.array(np.where(functype > 0))
        fit[j].noecl      = np.ones(fit[j].nobj)
        fit[j].ramp       = np.ones(fit[j].nobj)
        fit[j].rampuc     = np.ones(fit[j].nobjuc)
        if hasattr(fit[j].i, 'depth'):
            rampp[fit[j].i.depth + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'depth2'):
            rampp[fit[j].i.depth2 + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'depth3'):
            rampp[fit[j].i.depth3 + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'trqrprs'):
            rampp[fit[j].i.trqrprs + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'trq2rprs'):
            rampp[fit[j].i.trq2rprs + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'trrprs'):
            rampp[fit[j].i.trrprs + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'trrprs2'):
            rampp[fit[j].i.trrprs2 + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'rprs'):
            rampp[fit[j].i.rprs + numparams[cummodels[j]]] = 0.0
        if hasattr(fit[j].i, 'rprs2'):
            rampp[fit[j].i.rprs2 + numparams[cummodels[j]]] = 0.0
        k = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'ipmap':
                fit[j].noecl  *= fit[j].bestmip
                fit[j].ramp   *= fit[j].bestmip
                fit[j].rampuc *= fit[j].bestmipuc
            elif functype[i] == 'spline':
                fit[j].noecl  *= fit[j].bestspl
                fit[j].ramp   *= fit[j].bestspl
                fit[j].rampuc *= fit[j].bestspluc
            elif functype[i] == 'gp':
                fit[j].noecl  *= fit[j].bestgp
                fit[j].ramp   *= fit[j].bestgp
                fit[j].rampuc *= fit[j].bestgpuc
            #DO NOT REMOVE SINUSOIDAL FUNCTIONS FROM NORMALIZED RESULTS (FOR PHASE CURVES)
            elif functype[i] == 'sinusoidal':
                fit[j].noecl   *= myfuncs[i](rampp[iparams[i]], funcx  [i], fit[j].etc[k])
            else:
                fit[j].noecl   *= myfuncs[i](rampp[iparams[i]], funcx  [i], fit[j].etc[k])
                fit[j].ramp    *= myfuncs[i](rampp[iparams[i]], funcx  [i], fit[j].etc[k])
                fit[j].rampuc  *= myfuncs[i](rampp[iparams[i]], funcxuc[i], fit[j].etc[k])
            k += 1
        
        # normalize data
        fit[j].normflux      = (fit[j].flux       / fit[j].ramp).flatten()
        fit[j].normsigma     = (fit[j].newsigma   / fit[j].ramp).flatten()
        fit[j].normmeanfit   = (fit[j].meanfit    / fit[j].ramp).flatten()
        fit[j].normmedianfit = (fit[j].medianfit  / fit[j].ramp).flatten()
        fit[j].normbestfit   = (fit[j].bestfit    / fit[j].ramp).flatten()
        fit[j].normfluxuc    = (fit[j].fluxuc     / fit[j].rampuc).flatten()
        fit[j].normsigmauc   = (fit[j].newsigmauc / fit[j].rampuc).flatten()
        fit[j].normresuc     = (fit[j].fluxuc     / fit[j].rampuc2).flatten()
        fit[j].normresiduals = (fit[j].residuals  / fit[j].ramp.mean()).flatten()
        
        #COMPUTE SDNR
        fit[j].sdnr = np.std(fit[j].normresiduals)
        
        #SORT flux BY x, y AND radial POSITIONS
        yy      = np.sort(fit[j].position[0])
        xx      = np.sort(fit[j].position[1])
        yflux   = (fit[j].flux/fit[j].bestecl/fit[j].bestramp/fit[j].bestsin/fit[j].bestpos/ \
                   fit[j].bestvs/fit[j].bestff)[np.argsort(fit[j].position[0])]
        xflux   = (fit[j].flux/fit[j].bestecl/fit[j].bestramp/fit[j].bestsin/fit[j].bestpos/ \
                   fit[j].bestvs/fit[j].bestff)[np.argsort(fit[j].position[1])]
        ybestip = fit[j].bestip   [np.argsort(fit[j].position[0])] * \
                  fit[j].bestmip  [np.argsort(fit[j].position[0])] * \
                  fit[j].ballardip[np.argsort(fit[j].position[0])]
        xbestip = fit[j].bestip   [np.argsort(fit[j].position[1])] * \
                  fit[j].bestmip  [np.argsort(fit[j].position[1])] * \
                  fit[j].ballardip[np.argsort(fit[j].position[1])]
        
        #SORT flux BY frmvis
        fvsort = np.sort(fit[j].frmvis)
        vsflux = (fit[j].flux / fit[j].bestecl / fit[j].bestramp / fit[j].bestsin / fit[j].bestip / \
                  fit[j].bestpos / fit[j].bestff / fit[j].bestmip)[np.argsort(fit[j].frmvis)]
        
        #BIN DATA USING WEIGHTED AVERAGE
        fit[j].binmeanfit    = np.zeros(nbins[j])
        fit[j].binmedianfit  = np.zeros(nbins[j])
        fit[j].normbinfluxuc = np.zeros(nbins[j])
        fit[j].normbinflux   = np.zeros(nbins[j])
        fit[j].normbinsduc   = np.zeros(nbins[j])
        fit[j].normbinsd     = np.zeros(nbins[j])
        fit[j].normbinresuc  = np.zeros(nbins[j])
        fit[j].normbinbest   = np.zeros(nbins[j])
        fit[j].normbinmean   = np.zeros(nbins[j])
        fit[j].normbinmedian = np.zeros(nbins[j])
        fit[j].binnoecl      = np.zeros(nbins[j])
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
        fit[j].binresstd     = np.zeros(nbins[j])
        for i in range(nbins[j]):
            startuc                 = int(1.*i*fit[j].nobjuc/nbins[j])
            enduc                   = int(1.*(i+1)*fit[j].nobjuc/nbins[j])
            fit[j].normbinresuc[i]  = np.mean(fit[j].normresuc[startuc:enduc])
            fit[j].binstduc[i]      = np.sqrt(1 / sum(1/fit[j].newsigmauc[startuc:enduc]**2))
            fit[j].normbinfluxuc[i] = sum(fit[j].normfluxuc[startuc:enduc] /      \
                                          fit[j].normsigmauc[startuc:enduc]**2) / \
                                      sum(1/fit[j].normsigmauc[startuc:enduc]**2)
            fit[j].normbinsduc[i]   = np.sqrt(1 / sum(1/fit[j].normsigmauc[startuc:enduc]**2))
            start                   = int(1.*i*fit[j].nobj/nbins[j])
            end                     = int(1.*(i+1)*fit[j].nobj/nbins[j])
            fit[j].binbestfit[i]    = np.mean(fit[j].bestfit[start:end])
            fit[j].binmeanfit[i]    = np.mean(fit[j].meanfit[start:end])
            fit[j].binmedianfit[i]  = np.mean(fit[j].medianfit[start:end])
            fit[j].binnoecl[i]      = np.mean(fit[j].noecl[start:end])
            fit[j].normbinflux[i]   = sum(fit[j].normflux[start:end] /      \
                                          fit[j].normsigma[start:end]**2) / \
                                      sum(1/fit[j].normsigma[start:end]**2)
            fit[j].normbinbest[i]   = np.mean(fit[j].normbestfit[start:end])
            fit[j].normbinmean[i]   = np.mean(fit[j].normmeanfit[start:end])
            fit[j].normbinmedian[i] = np.mean(fit[j].normmedianfit[start:end])
            fit[j].normbinsd[i]     = np.sqrt(1 / sum(1/fit[j].normsigmauc[start:end]**2))
            fit[j].binvsflux[i]     = np.median(vsflux[start:end])
            fit[j].binvsflstd[i]    = np.std(vsflux[start:end]) / np.sqrt(end-start)
            fit[j].binfrmvis[i]     = np.median(fvsort[start:end])
            #fit[j].binxx[i]         = np.mean(xx[start:end])
            #fit[j].binyy[i]         = np.mean(yy[start:end])
            #fit[j].binxflux[i]      = np.median(xflux[start:end])
            #fit[j].binyflux[i]      = np.median(yflux[start:end])
            #fit[j].binxflstd[i]     = np.std(xflux[start:end]) / np.sqrt(end-start)
            #fit[j].binyflstd[i]     = np.std(yflux[start:end]) / np.sqrt(end-start)
            xxrange                =      xx[np.where(np.bitwise_and((xx >= xx[0] + 1.*i    *(xx[-1]-xx[0])/nbins[j]),\
                                                                     (xx <= xx[0] + 1.*(i+1)*(xx[-1]-xx[0])/nbins[j])))]
            xfluxrange             =   xflux[np.where(np.bitwise_and((xx >= xx[0] + 1.*i    *(xx[-1]-xx[0])/nbins[j]),\
                                                                     (xx <= xx[0] + 1.*(i+1)*(xx[-1]-xx[0])/nbins[j])))]
            xbestiprng             = xbestip[np.where(np.bitwise_and((xx >= xx[0] + 1.*i    *(xx[-1]-xx[0])/nbins[j]),\
                                                                     (xx <= xx[0] + 1.*(i+1)*(xx[-1]-xx[0])/nbins[j])))]
            yyrange                =      yy[np.where(np.bitwise_and((yy >= yy[0] + 1.*i    *(yy[-1]-yy[0])/nbins[j]),\
                                                                     (yy <= yy[0] + 1.*(i+1)*(yy[-1]-yy[0])/nbins[j])))]
            yfluxrange             =   yflux[np.where(np.bitwise_and((yy >= yy[0] + 1.*i    *(yy[-1]-yy[0])/nbins[j]),\
                                                                     (yy <= yy[0] + 1.*(i+1)*(yy[-1]-yy[0])/nbins[j])))]
            ybestiprng             = ybestip[np.where(np.bitwise_and((yy >= yy[0] + 1.*i    *(yy[-1]-yy[0])/nbins[j]),\
                                                                     (yy <= yy[0] + 1.*(i+1)*(yy[-1]-yy[0])/nbins[j])))]
            fit[j].binxx[i]        = np.mean(xxrange)
            fit[j].binyy[i]        = np.mean(yyrange)
            fit[j].binxflux[i]     = np.mean(xfluxrange)
            fit[j].binyflux[i]     = np.mean(yfluxrange)
            fit[j].binxflstd[i]    = np.std(xfluxrange) / np.sqrt(xfluxrange.size)
            fit[j].binyflstd[i]    = np.std(yfluxrange) / np.sqrt(yfluxrange.size)
            fit[j].binxbestip[i]   = np.mean(xbestiprng)
            fit[j].binybestip[i]   = np.mean(ybestiprng)
            fit[j].binxbipstd[i]   = np.std(xbestiprng) / np.sqrt(xbestiprng.size)
            fit[j].binybipstd[i]   = np.std(ybestiprng) / np.sqrt(ybestiprng.size)
            fit[j].binres[i]       = np.mean(fit[j].residuals[start:end])
            fit[j].binresstd[i]    = np.sqrt(1 / sum(1/fit[j].residuals[start:end]**2))
        
        fit[j].normbinres   = fit[j].normbinflux/fit[j].normbinbest - 1.
        #fit[j].normbinresuc = fit[j].binrampuc2/fit[j].normbinbest - 1.
        '''
        #Create savefile for ipPlotting.py
        savefile = event[j].modeldir + "/d-" + event[j].eventname + "-ip-" + fit[j].saveext + ".npz"
        np.savez(savefile, x=fit[j].position[1], y=fit[j].position[0], 
                 flux=fit[j].flux/fit[j].bestecl/fit[j].bestramp/fit[j].bestsin/fit[j].bestpos/fit[j].bestvs, 
                 bestpip=fit[j].bestpip, ydiv=event[j].params.ydiv, xdiv=event[j].params.xdiv)
        '''

    #CALCULATE AIC/BIC VALUE FOR JOINT MODEL FITS
    if numevents > 1:
        aic  = 0
        bic  = 0
        nobj = 0
        for j in range(numevents):
            aic  += sum((fit[j].residuals/event[j].fit[0].newsigma)**2)
            bic  += sum((fit[j].residuals/event[j].fit[0].newsigma)**2)
            nobj += fit[j].nobj
        aic += 2*numfreepars
        bic +=   numfreepars*np.log(nobj)
        print("AIC for joint model fit = " + str(aic), file=printout)
        print("BIC for joint model fit = " + str(bic), file=printout)
    
    #Assign abscissa time unit (orbits or days)
    for j in range(numevents):
        if event[j].params.timeunit == 'orbits':
            fit[j].tuall      = event[j].phase.flatten()    #Include all frames
            fit[j].timeunituc = fit[j].phaseuc              #Use good, unclipped frames
            fit[j].timeunit   = fit[j].phase                #Use only clipped frames
            fit[j].abscissa   = fit[j].binphase             #Binned clipped frames
            fit[j].abscissauc = fit[j].binphaseuc           #Binned unclipped
            fit[j].xlabel     = 'Orbital Phase'
        elif event[j].params.timeunit == 'days-utc':
            fit[j].tuall      = event[j].bjdutc.flatten() - event[j].params.tuoffset
            fit[j].timeunituc = fit[j].bjdutcuc    - event[j].params.tuoffset
            fit[j].timeunit   = fit[j].bjdutc      - event[j].params.tuoffset
            fit[j].abscissa   = fit[j].binbjdutc   - event[j].params.tuoffset
            fit[j].abscissauc = fit[j].binbjdutcuc - event[j].params.tuoffset
            fit[j].xlabel     = 'BJD_UTC - ' + str(event[j].params.tuoffset)
        elif event[j].params.timeunit == 'days-tdb':
            fit[j].tuall      = event[j].bjdtdb.flatten() - event[j].params.tuoffset
            fit[j].timeunituc = fit[j].bjdtdbuc    - event[j].params.tuoffset
            fit[j].timeunit   = fit[j].bjdtdb      - event[j].params.tuoffset
            fit[j].abscissa   = fit[j].binbjdtdb   - event[j].params.tuoffset
            fit[j].abscissauc = fit[j].binbjdtdbuc - event[j].params.tuoffset
            fit[j].xlabel     = 'BJD_TDB - ' + str(event[j].params.tuoffset)
        else:
            fit[j].tuall      = event[j].bjdutc.flatten() - event[j].params.tuoffset
            fit[j].timeunituc = fit[j].bjdutcuc    - event[j].params.tuoffset
            fit[j].timeunit   = fit[j].bjdutc      - event[j].params.tuoffset
            fit[j].abscissa   = fit[j].binbjdutc   - event[j].params.tuoffset
            fit[j].abscissauc = fit[j].binbjdutcuc - event[j].params.tuoffset
            fit[j].xlabel     = 'BJD - ' + str(event[j].params.tuoffset)
    
    print('Producing figures.')
    #PLOT BINNED DATA AND BEST FIT
    for j in range(numevents):
        fignum   = 6001+num*numfigs+j*100
        savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
        plots.binlc(event[j], fit[j], fignum, savefile, j=j)

    #PLOT NORMALIZED BINNED DATA AND BEST FIT
    for j in range(numevents):
        fignum   = 6002+num*numfigs+j*100
        savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
        if hasattr(event[j].params, 'interclip'):
            interclip = event[j].params.interclip
        else:
            interclip = None
        plots.normlc(event[j], fit[j], fignum, savefile, j=j, interclip=interclip)
    
    if allplots == True:
        #allparams TRACE PARAMETER VALUES FOR ALL STEPS
        for j in range(numevents):
            fignum   = 6003+num*numfigs+j*100
            savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
            plots.trace(event[j], fit[j], fignum, savefile)
        
        #allorthop TRACE PARAMETER VALUES FOR ALL STEPS
        for j in range(numevents):
            if event[j].isortho == True and event[j].params.newortho == False:
                fignum   = 6013+num*numfigs+j*100
                savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
                plots.trace(event[j], fit[j], fignum, savefile=savefile, allparams=fit[j].allorthop, parname=fit[j].opname)
        
        #allparams AUTOCORRELATION PLOT
        for j in range(numevents):
            fignum   = 6004+num*numfigs+j*100
            savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
            plots.autocorr(event[j], fit[j], fignum, savefile)
        
        #allorthop AUTOCORRELATION PLOT
        for j in range(numevents):
            if event[j].isortho == True and event[j].params.newortho == False:
                fignum   = 6014+num*numfigs+j*100
                savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
                plots.autocorr(event[j], fit[j], fignum, savefile=savefile, allparams=fit[j].allorthop, parname=fit[j].opname)
        
        #allparams CORRELATION PLOTS WITH 2D HISTOGRAMS
        for j in range(numevents):
            fignum   = 6005+num*numfigs+j*100
            savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
            try:
                plots.hist2d(event[j], fit[j], fignum, savefile=savefile)
            except:
                print("Failed to create Figure " + str(fignum))
        
        #allorthop CORRELATION PLOTS WITH 2D HISTOGRAMS
        for j in range(numevents):
            if event[j].isortho == True and event[j].params.newortho == False:
                fignum   = 6015+num*numfigs+j*100
                savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
                plots.hist2d(event[j], fit[j], fignum, savefile=savefile, allparams=fit[j].allorthop, parname=fit[j].opname)
        
        #allparams 1D HISTOGRAMS
        for j in range(numevents):
            fignum   = 6006+num*numfigs+j*100
            savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
            plots.histograms(event[j], fit[j], fignum, savefile)
        
        #allorthop 1D HISTOGRAMS
        for j in range(numevents):
            if event[j].isortho == True and event[j].params.newortho == False:
                fignum   = 6016+num*numfigs+j*100
                savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
                plots.histograms(event[j], fit[j], fignum, savefile=savefile, allparams=fit[j].allorthop, parname=fit[j].opname)
        
        #Plot projections of position sensitivity along x and y
        for j in range(numevents):
            fignum   = 6007+num*numfigs+j*100
            savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
            plots.ipprojections(event[j], fit[j], fignum, savefile=savefile)
        
        #BLISS MAP
        for j in range(numevents):
            if fit[j].isipmapping:
                # Minimum number of acceptable points in a bin
                if event[j].params.minnumpts.__class__ == int:
                    minnumpts = event[j].params.minnumpts
                elif len(event[j].params.minnumpts) == len(event[j].params.model):
                    minnumpts = event[j].params.minnumpts[num]
                else:
                    minnumpts = event[j].params.minnumpts[0]
                # Plot BLISS Map
                fignum   = 6008+num*numfigs+j*100
                savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
                plots.blissmap(event[j], fit[j], fignum, savefile=savefile, minnumpts=minnumpts)
                # Plot Pointing Histogram
                fignum   = 6009+num*numfigs+j*100
                savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
                plots.pointingHist(event[j], fit[j], fignum, savefile=savefile, minnumpts=minnumpts)
        
        #PLOT RMS vs. BIN SIZE
        for j in range(numevents):
            if hasattr(event[j].params, 'rmsbins'):
                fit[j].rms, fit[j].stderr, fit[j].binsz = cn.computeRMS(fit[j].normresiduals, binstep=1,
                                                                        maxnbins = event[j].params.rmsbins)
            else:
                fit[j].rms, fit[j].stderr, fit[j].binsz = cn.computeRMS(fit[j].normresiduals, binstep=1)
            # Compute standard error for noisy data, instead of denoised data
            if hasattr(event[j].params, 'noisysdnr') and event[j].params.noisysdnr != None:
                fit[j].stderr = cn.computeStdErr(event[j].params.noisysdnr, fit[j].normresiduals.size, fit[j].binsz)
            fignum   = 6011+num*numfigs+j*100
            savefile = event[j].modeldir+"/"+event[j].eventname+"-fig"+str(fignum)+"-"+ fit[j].saveext+".png"
            try:
                plots.rmsplot(event[j], fit[j], fignum, savefile=savefile)
            except:
                print("Failed to create Figure " + str(fignum))
        
        #PLOT VISIT SENSITIVITY AND MODEL
        if 'vissen' in functype:
            numvsmodels = 0
            for i in functype:
                if i == 'vissen':
                    numvsmodels += 1
            
            plt.figure(612+num*numfigs)
            plt.clf()
            a = plt.suptitle(obj + ' Visit # vs. Sensitivity', size=16)
            k = 1
            for j in range(numevents):
                for i in range(cummodels[j],cummodels[j+1]):
                    if functype[i] == 'vissen':
                        a = plt.subplot(numvsmodels,1,k)
                        a = plt.errorbar(fit[j].binfrmvis, fit[j].binvsflux, fit[j].binvsflstd, fmt='go', label='Data')
                        a = plt.plot(fit[j].frmvis, fit[j].bestvs, 'k.', label='Model')
                        for model in fit[j].model:
                            if (model == 'vsspline'):
                                a = plt.plot(event[j].params.vsknots, fit[j].bestpvs, 'bs')
                        a = plt.ylabel('Flux Sensitivity')
                        a = plt.xlabel('Visit Number')
                        a = plt.legend(loc='best')
                        a.fontsize=8
            
            plt.savefig(event[j].modeldir + "/" + obj + "-fig" + str(num*numfigs+1608) + "-" + saveext + ".png")
    
    if isinteractive == False:
        plt.close('all')
    else:
        plt.pause(0.01)     #Executes draw command for figures
    
    if hasattr(event[0].params, 'savedata') and event[0].params.savedata == False:
        pass
    else:
        print('Writing save files.')
        for j in range(numevents):
            #WRITE fit[j].allparams TO FILE, DELETE fit[j].allparams
            fit[j].allparamsfile = event[j].modeldir + "/d-" + event[j].eventname + '-allparams-' + \
                                   fit[j].saveext + '.npy'
            pid  = open(fit[j].allparamsfile, 'w')
            np.save(pid, fit[j].allparams)
            pid.close()
            del fit[j].allparams, fit[j].allorthop
            
            if event[j].params.newortho == True:
                # Write ortho save file
                fit[j].orthofile = "d-" + event[j].eventname + '-ortho-' + fit[j].saveext + '.npz'
                pid  = open(fit[j].orthofile, 'w')
                np.savez(pid, invtrans=fit[j].invtrans, trans=fit[j].trans, origin=fit[j].origin, sigma=fit[j].orthosigma)
                pid.close()
    
    return

