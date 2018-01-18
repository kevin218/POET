'''
This package includes functions for calculating fractional flux at time/phase points for a non-linear limb darkening model given specific input:
    trnlld()
'''

import numpy as np
from trlduniform import trlduniform

def trnlld(params, t, etc):
    '''
    This function computes the primary transit shape using non-linear limb-darkening equations for a
    planet granted a non-linear function, as provided by Mandel & Agol (2002).

    Parameters
    ----------
    midpt:  Center of eclipse
    rprs:   Planet radius / stellar radius
    cosi:   Cosine of the inclination
    ars:    Semi-major axis / stellar radius
    flux:   Flux offset from 0
    c#:     Limb-darkening coefficients
    t:	    Array of phase/time points
    p:      Period in same units as t

    Returns
    -------
    This function returns the flux for each point in t.

    References
    ----------
    
    Mandel & Agol (2002)
    /home/esp01/doc/Mandel+Agol-2002_eq8.pdf
    /home/esp01/code/MandelAgol/occultnl.pro
    
    Revisions
    ---------
    2011-04-19  M. Oliver Bowman  
                bowman@knights.ucf.edu
                Converted to Python
    '''
    
    #DEFINE PARAMETERS
    midpt, rprs, cosi, ars, flux, p, c1, c2, c3, c4 = params
    
    #COMPUTE z(t) FOR TRANSIT ONLY (NOT ECLIPSE) AND Sigma*4
    #NOTE: z(t) ASSUMES A CIRCULAR ORBIT
    z = ars * np.sqrt(np.sin(2 * np.pi * (t - midpt)/p)**2 + (cosi * np.cos(2 * np.pi * (t - midpt)/p))**2)
    z[np.where(np.bitwise_and((t - midpt)%p > p/4.,(t - midpt) % p < p * 3./4))] = ars
    
    # make grid in radius:
    # Call magnification of uniform source:
    mulimb0 = trlduniform(midpt, rprs, cosi, ars, flux, p, t, etc)
    fac = max(abs(mulimb0 - 1))
    
    # Calculate omega
    omega = 4. * ((1. - c1 - c2 - c3 - c4) / 4. + c1 / 5. + c2 / 6. + c3 / 7. + c4 / 8.)
    
    # Initialize indexes and flux/time array
    nb = len(z)
    indx = np.where(mulimb0 != 1)[0]
    mulimb = mulimb0[indx]
    mulimbf = np.zeros((nb, 5))
    mulimbf[:, 0] += 1.
    mulimbf[:, 1] += 0.8
    mulimbf[:, 2] += 2. / 3.
    mulimbf[:, 3] += 4. / 7.
    mulimbf[:, 4] += 0.5
    nr = 2
    dmumax = 1
    
    # Loop for all points t
    while (dmumax > (fac * .001)):
        mulimbp = mulimb
        nr *= 2
        dx  = 0.5 * np.pi / nr
        x   = dx * np.arange(nr + 1)
        th  = x + 0.5 * dx
        r   = np.sin(x)
        sig = np.sqrt(np.cos(th[nr - 1]))
        mulimbhalf  = sig ** 3 * mulimb0[indx] / (1 - r[nr - 1])
        mulimb1     = sig ** 4 * mulimb0[indx] / (1 - r[nr - 1])
        mulimb3half = sig ** 5 * mulimb0[indx] / (1 - r[nr - 1])
        mulimb2     = sig ** 6 * mulimb0[indx] / (1 - r[nr - 1])
        for i in range(1, nr):
            # Calculate uniform magnification at intermediate radii:
            mu = trlduniform(midpt, rprs/r[i], cosi, ars/r[i], flux, p, t[indx], etc=[])
            sig1 = np.sqrt(np.cos(th[i - 1]))
            sig2 = np.sqrt(np.cos(th[i]))
            mulimbhalf  = mulimbhalf  + r[i] ** 2 * mu * (sig1 ** 3 / (r[i] - r[i-1]) - sig2 ** 3 / (r[i+1] - r[i]))
            mulimb1     = mulimb1     + r[i] ** 2 * mu * (sig1 ** 4 / (r[i] - r[i-1]) - sig2 ** 4 / (r[i+1] - r[i]))
            mulimb3half = mulimb3half + r[i] ** 2 * mu * (sig1 ** 5 / (r[i] - r[i-1]) - sig2 ** 5 / (r[i+1] - r[i]))
            mulimb2     = mulimb2     + r[i] ** 2 * mu * (sig1 ** 6 / (r[i] - r[i-1]) - sig2 ** 6 / (r[i+1] - r[i]))
        
        mulimb = ((1 - c1 - c2 - c3 - c4) * mulimb0[indx] + c1 * mulimbhalf * dx + c2 * mulimb1 * dx + c3 * mulimb3half * dx + c4 * mulimb2 * dx) / omega
        ix1 = np.where((mulimb + mulimbp) != 0)[0]
        dmumax = max(abs(mulimb[ix1] - mulimbp[ix1]) / (mulimb[ix1] + mulimbp[ix1]))
    
    # Adjust mulimbf if needed for later use 
    mulimbf[indx, 0] = mulimb0[indx]
    mulimbf[indx, 1] = mulimbhalf * dx
    mulimbf[indx, 2] = mulimb1 * dx
    mulimbf[indx, 3] = mulimb3half * dx
    mulimbf[indx, 4] = mulimb2 * dx
    mulimb0[indx] = mulimb
    return mulimb0

