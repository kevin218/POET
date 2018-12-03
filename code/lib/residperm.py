
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import sys, nasc

def modelfunc(freepars, event, fit, params, pmin, pmax, flux, cummodels, functype, myfuncs, funcx, iparams, inonprior, stepsize, funcxuc, nights):
    '''
    Compute normalized residuals for least-squares fit
    '''
    params[inonprior] = freepars
    for i in range(len(stepsize)):
        if stepsize[i] < 0:
            params[i] = params[-stepsize[i]-1]
    params[np.where(params < pmin)] = pmin[np.where(params < pmin)]
    params[np.where(params > pmax)] = pmax[np.where(params > pmax)]
    
    for k in np.unique(nights):
        tonight   = np.where(nights == k)[0]
        # Construct non-analytic systematic model
        nasc.nasmodel(event, fit, params, iparams, cummodels, functype, myfuncs, funcxuc, tonight)
            
    residuals = []
    bestfit   = []
    #pedit     = np.copy(params)
    numevents = len(fit)
    for j in range(numevents):
        fit0    = np.ones(fit[j].nobj)
        k       = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'ipmap':
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
            #elif functype[i] == 'noise':
            #    pass
            else:
                subparams = params[iparams[i]]
                fit0 *= myfuncs[i](subparams, funcx[i], fit[j].etc[k])
                params[iparams[i]] = subparams
            k    += 1
        #Multiply analytic model by non-analytic systematics model
        fit0 *= fit[j].bestsys
        residuals = np.concatenate((residuals,(fit0 - flux[j])/fit[j].sigma),axis=0)
    return residuals

def modelfit(event, fit, params, flux, cummodels, functype, myfuncs, funcx, iparams, funcxuc, nights):
    '''
    Compute model fit and residuals using supplied paramter values.
    '''
    #params[inonprior] = freepars
    #params[np.where(params < pmin)] = pmin[np.where(params < pmin)]
    #params[np.where(params > pmax)] = pmax[np.where(params > pmax)]
    for k in np.unique(nights):
        tonight   = np.where(nights == k)[0]
        # Construct non-analytic systematic model
        nasc.nasmodel(event, fit, params, iparams, cummodels, functype, myfuncs, funcxuc, tonight)
            
    residuals = []
    bestfit   = []
    #pedit     = np.copy(params)
    numevents = len(fit)
    for j in range(numevents):
        fit0    = np.ones(fit[j].nobj)
        k       = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pass
            elif functype[i] == 'noise':
                pass
            elif functype[i] == 'ipmap':
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
            #elif functype[i] == 'noise':
            #    pass
            else:
                subparams = params[iparams[i]]
                fit0 *= myfuncs[i](subparams, funcx[i], fit[j].etc[k])
                params[iparams[i]] = subparams
            k    += 1
        #Multiply analytic model by non-analytic systematics model
        fit0 *= fit[j].bestsys
        residuals.append(fit0 - flux[j])
        bestfit.append(fit0)
    return bestfit, residuals

'''
bestfit = np.ones(100)
residuals = np.random.normal(0,0.3,100)

'''

def permute2(event, fit, params, pmin, pmax, flux, cummodels, functype, myfuncs, funcx, funcxuc, iparams, stepsize, inonprior, nights, maxnumit=10000):
    '''
    
    '''
    print('Performing error estimation using the residual permutation technique.')
    bestfit, residuals = modelfit(event, fit, params, flux, cummodels, functype, myfuncs, funcx, iparams, funcxuc, nights)
    
    numevents    = len(fit)
    newflux      = [[] for j in range(numevents)]
    nump         = len(params)
    nnights      = len(np.unique(nights))
    numit        = np.zeros(nnights)+maxnumit
    
    # Determine numit for each night
    for j in range(len(residuals)):
        numit[nights[j]] = np.min((numit[nights[j]],len(residuals[j])))
    
    # Form indices
    totnumit = np.product(numit)
    if totnumit <= 1e9:
        ind      = np.zeros((totnumit,nnights),dtype=np.int)
        counter  = 1
        for j in range(nnights):
            ind[:,j] = (np.arange(totnumit)/counter)%numit[j]
            counter *= numit[j]
        # Permute indices and crop to maxnumit
        ind = np.random.permutation(ind)[:maxnumit]
        totnumit =  int(np.min((totnumit,maxnumit)))
    else:
        totnumit = maxnumit
        ind      = np.zeros((totnumit,nnights),dtype=np.int)
        for j in range(nnights):
            ind[:,j] = np.random.randint(0,numit[j],totnumit)
    
    allparams    = np.zeros((nump,totnumit))
    for i in range(totnumit):
        sys.stdout.write('\r'+str(i+1)+'/'+str(totnumit)+'      ')
        sys.stdout.flush()
        for j in range(numevents):
            # Create new realization of data 
            k          = nights[j]
            newflux[j] = bestfit[j] + np.concatenate((residuals[j][ind[i,k]:],residuals[j][:ind[i,k]]))
        '''
        plt.figure(2)
        plt.plot(newflux[0],'o')
        plt.pause(3)
        '''
        # Calculate least-squares best fit
        output = op.leastsq(modelfunc, params[inonprior], \
            args=(event, fit, params, pmin, pmax, newflux, cummodels, functype, myfuncs, funcx, iparams, inonprior, stepsize, funcxuc, nights), \
            factor=100, ftol=1e-8, xtol=1e-8, gtol=1e-8, diag=1./stepsize[inonprior], full_output=False)
        params[inonprior] = output[0]
        allparams[:,i]    = params
    
    # Update parameters
    for i in range(len(stepsize)):
        if stepsize[i] < 0:
            allparams[i] = allparams[-stepsize[i]-1]
    
    return allparams

#########################################################################
#                               OLD                                     #
#########################################################################

#Residual permutation
def permute(fit, params, pmin, pmax, bestfit, residuals, newflux, cummodels, functype, myfuncs, funcx, iparams, stepsize, inonprior, ispermute, maxnumit=10000, nk=1, nnights=1, ni=1, nnumit=1):
    '''
    
    '''
    # Set number of iterations to length of shortest dataset
    numit          = maxnumit
    for i in range(len(residuals)):
        numit = np.min((numit,len(residuals[i])))
    nump           = len(params)
    numevents      = len(fit)
    allparams      = np.zeros((nump,numit))
    allparams[:,0] = params
    for i in range(1,numit):
        for j in range(numevents):
            if ispermute[j]:
                # Permute residuals to create new realization of data
                newflux[j] = bestfit[j] + np.concatenate((residuals[j][i:],residuals[j][:i]))
            #else:
            #    newflux.append(bestfit[j])
        #print(i,numit)
        sys.stdout.write('\r'+str(nk+1)+'/'+str(nnights)+' '+str(ni+1)+'/'+str(nnumit)+' '+str(i+1)+'/'+str(numit)+'      ')
        sys.stdout.flush()
        '''
        plt.figure(2)
        plt.plot(newflux[0],'o')
        plt.pause(3)
        '''
        # Calculate least-squares best fit
        output = op.leastsq(modelfunc, params[inonprior], \
            args=(event, fit, params, pmin, pmax, newflux, cummodels, functype, myfuncs, funcx, iparams, inonprior, stepsize, funcxuc, nights), \
            factor=100, ftol=1e-8, xtol=1e-8, gtol=1e-8, diag=1./stepsize[inonprior], full_output=False)
        allparams[inonprior,i] = output[0]
    
    return allparams

def driver(fit, params, pmin, pmax, flux, cummodels, functype, myfuncs, funcx, iparams, stepsize, inonprior, nights, maxnumit=10000):
    '''
    
    '''
    print('Performing error estimation using the residual permutation technique.')
    bestfit, residuals = modelfit(fit, params, flux, cummodels, functype, myfuncs, funcx, iparams)
    '''
    plt.figure(1)
    plt.clf()
    for i in range(len(residuals)):
        plt.plot(bestfit[i] + residuals[i],'o')
    plt.pause(3)
    '''
    numevents    = len(fit)
    newflux      = [[] for j in range(numevents)]
    nump         = len(params)
    nnights      = len(np.unique(nights))
    if nnights == 1:
        ispermute = np.ones(numevents, dtype=bool)
        allparams = permute(fit, params, pmin, pmax, bestfit, residuals, newflux, cummodels, functype, myfuncs, funcx, iparams, stepsize, inonprior, ispermute, maxnumit)
    else:
        allparams = np.empty(shape=(nump,0))
        for k in range(nnights-1):
            tonight   = (np.where(nights == k)[0])#.tolist()
            ispermute = np.ones(numevents, dtype=bool)
            ispermute[tonight] = False
            # Set number of iterations to length of shortest dataset
            numit          = maxnumit
            for j in range(len(residuals)):
                if ispermute[j] == False:
                    numit = np.min((numit,len(residuals[j])))
            # Permute residuals for the ith night then fix while permuting all other nights
            for i in range(numit):
                for j in tonight:
                    newflux[j] = bestfit[j] + np.concatenate((residuals[j][i:],residuals[j][:i]))
                allp = permute(fit, params, pmin, pmax, bestfit, residuals, newflux, cummodels, functype, myfuncs, funcx, iparams, stepsize, inonprior, ispermute, maxnumit, k, nnights-1, i, numit)
                allparams = np.concatenate((allparams,allp),axis=1)
    
    #allparams = np.reshape(allparams,(nump, (m+1)*nchains))
    for i in range(len(stepsize)):
        if stepsize[i] < 0:
            allparams[i] = allparams[-stepsize[i]-1]
    return allparams

