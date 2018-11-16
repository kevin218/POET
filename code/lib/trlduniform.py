'''
This package includes functions for calculating fractional flux at time/phase points for a uniform limb darkening model given specific input:
    trlduniform()
'''

import numpy as np

def trlduniform(midpt, rprs, cosi, ars, flux, p, t, etc):
    '''
    This function computes the primary transit shape using a uniform limb-darkening equation, as provided by Mandel & Agol (2002).

    Parameters
    ----------
    midpt:  Center of eclipse
    rprs:   Planet radius / stellar radius
    cosi:   Cosine of the inclination
    ars:    Semi-major axis / stellar radius
    flux:   Flux offset from 0
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
    
    # Set exception case for avoiding invalid boundreis
    if(abs(rprs - .5) < .001):
        rprs = .5
    
    #COMPUTE z(t) FOR TRANSIT ONLY (NOT ECLIPSE) AND Sigma*4
    #NOTE: z(t) ASSUMES A CIRCULAR ORBIT
    z = ars * np.sqrt(np.sin(2 * np.pi * (t - midpt)/p)**2 + (cosi * np.cos(2 * np.pi * (t - midpt)/p))**2)
    z[np.where(np.bitwise_and((t - midpt)%p > p/4.,(t - midpt)%p < p * 3./4))] = ars
    
    #CALCULATE TRANSIT SHAPE WITH LIMB-DARKENING
    y = np.ones(len(t))
    if rprs == 0:
        return y*flux
    
    # Where planet is not in transit
    miss        = np.where(z > (1+rprs))
    y[miss]     = 1
    
    # Where planet is in ingress/egress
    iingress    = np.where( np.bitwise_and( np.bitwise_and( (1-rprs) < z, z <= (1+rprs) ), z > (rprs- 1)) )[0]
    kap1        = np.arccos((1 - rprs**2 + z[iingress]**2) / 2 / z[iingress])
    kap0        = np.arccos((rprs**2 + z[iingress]**2 - 1) / 2 / rprs / z[iingress])
    lambdae     = (rprs**2 * kap0 + kap1 - 0.5 * np.sqrt(4 * z[iingress]**2 - (1 + z[iingress]**2 - rprs**2)**2)) / np.pi
    y[iingress] = 1 - lambdae
    
    # Where planet is in transit
    itrans      = np.where(z <= (1-rprs))
    y[itrans]   = 1 - rprs**2
    
    # If planet covers source
    cover       = np.where(z <= rprs - 1)
    y[cover]    = 0
    
    return y*flux
    
