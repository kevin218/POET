

import scipy.optimize as op

def modelfunc(freepars, params, pmin, pmax, sigma):
    params[ifreepars] = freepars
    params[np.where(params < pmin)] = pmin[np.where(params < pmin)]
    params[np.where(params > pmax)] = pmax[np.where(params > pmax)]
    residuals = []
    for j in range(numevents):
        fit0 = np.ones(fit[j].nobj)
        k    = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'ipmap':
                fit[j].etc[k] = fit0
            fit0 *= myfuncs[i](params[iparams[i]], funcx[i], fit[j].etc[k])
            k    += 1
        residuals = np.concatenate((residuals,(fit0 - flux[j])/sigma[j]),axis=0)
    return residuals

def minimize(params, pmin, pmax, sigma):
    #modelfunc = lambda freepars, params, sigma: callmodelfunc(freepars, params, sigma)
    return op.leastsq(modelfunc, params[ifreepars], args=(params, pmin, pmax, sigma), factor=100, \
                             ftol=1e-16, xtol=1e-16, gtol=1e-16, diag=1./stepsize[ifreepars])

#params, pmin, pmax, sigma, ifreepars, numevents, fit[j].nobj, cummodels, functype, fit[j].etc, myfuncs, iparams, funcx, flux
