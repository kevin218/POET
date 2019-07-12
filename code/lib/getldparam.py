'''
This package includes functions for calculating the four parameters of a 
non-linear limb darkening model given specific input:
    getlimbparam()
'''

import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as intr
from integrate import integrate
from numpy.linalg import lstsq

def getlimbparam(kuruczfile, chan, Nparam=4, lamfact=1000., spit=None, returnmu=False):
    """
    Reads in Kurucz grid data and Spitzer specral responses to return four 
    parameters for a non-linear limb darkening model.
    
    This function requires a Kurucz grid of flux responses for a star over 
    effective temperature, specific gravity, metalicity and angle. The return 
    parameters are to fit the following model,
    
        I(mu) / I(1) = 1 - a1(1 - mu^0.5) - a2(1 - mu) 
                         - a3(1 - mu^1.5) - a4(1 - mu^2)
    
    where mu is cos(angle) from star's center, I(mu) is flux/intensity given 
    angle, and a# are the four parameters.
    
    Paramaters
    ----------
    kuruczfile : str
        String giving the location of the Kurucz grid (Only for the target 
        star's TEFF, log(g), M/H (or Fe/H))
        ex: kuruczfile = "/home/esp01/ancil/limbdarkening/WASP-11"
    chan : int
        The Spitzer channel to be used for model fitting
        ex: chan       = 1
    
    Returns
    -------
    coeff : 1D array
        Contains four parameters for the non-linear limb darkening model as 
        follows:

        I(mu) / I(1) = 1 - a1(1 - mu^[1/2]) - a2(1 - mu) 
                         - a3(1 - mu^[3/2]) - a4(1 - mu^2)
        
        coeff array is in format
        array([ a1, a2, a3, a4])
    
    Examples
    --------
    >>> from getlimbparam import getlimbparam
    >>> kuruczfile = "/home/esp01/ancil/limbdarkening/WASP-11"
    >>> chan = 2
    >>> paramschan2 = getlimbparam(kuruczfile, chan)
    >>> print(paramschan2)
    [ 0.64940497 -0.85966835  0.85679169 -0.31649174]
    >>> chan = 3
    >>> paramschan3 = getlimbparam(kuruczfile, chan)
    >>> print(paramschan3)
    [ 0.60509989 -0.86795052  0.86195088 -0.31567865]
    
    Revisions
    ---------
    2011-02-22 0.0 bowman@knights.ucf.edu Inital version.
    2011-03-22 0.1 bowman@knights.ucf.edu Added documentation, cleaned up 
                                          code style.
    2011-05-11 0.2 bowman@knights.ucf.edu Further cleaned up style.
    """
    
    # Define conversion factors
    kurconv = 100000.
    
    # Detect proper spitzer file
    if   chan == 1:
        spitfile = "/home/kevin/Documents/esp01/ancil/filter/irac_tr1_2004-08-09.dat"
    elif chan == 2:
        spitfile = "/home/kevin/Documents/esp01/ancil/filter/irac_tr2_2004-08-09.dat"
    elif chan == 3:
        spitfile = "/home/kevin/Documents/esp01/ancil/filter/irac_tr3_2004-08-09.dat"
    elif chan == 4:
        spitfile = "/home/kevin/Documents/esp01/ancil/filter/irac_tr4_2004-08-09.dat"
    
    # Retrieve data from Kurucz grid and Spitzer spectral response table for 
    # the appropriate channel
    
    # First column of kurucz is the wavelength
    kurucz  = np.loadtxt(kuruczfile, skiprows = 3, unpack=True)
    isfilter = False
    if spit == None:
        isfilter = True
        spit    = np.loadtxt(spitfile, unpack=True)
        spit    = spit[0:2] # Use second column for spectral response
    
    # Convert kurucz grid flux to standard units
    for i in range(np.shape(kurucz)[1]):
        for j in range(np.shape(kurucz)[0] - 2):
            kurucz[j + 2, i] *= kurucz[1, i] / kurconv
    
    # Find where kurucz is approx. equal to spit's min and max
    lowspit  = spit[0].min() * lamfact
    highspit = spit[0].max() * lamfact
    
    #lowspitin  = np.where(spit[0] == spit[0].min())[0]
    #highspitin = np.where(spit[0] == spit[0].max())[0]
    
    lowkur  = np.abs(kurucz[0] -  lowspit).argmin()
    highkur = np.abs(kurucz[0] - highspit).argmin()
    
    # Cut kurucz grid to match wavelength range with spitzer
    kurucznew = kurucz[:, lowkur:highkur]
    kurucznew[0] /= lamfact
    
    # Spline for spitzer data fit at exact kurucz wavelengths
    points  = np.size(kurucznew[0])
    
    spitfit = np.zeros((2, points))
    spitfit[0] = kurucznew[0]
    
    for i in range(points):
        tck = intr.splrep(spit[0], spit[1])
        spitfit[1, i] = intr.splev(spitfit[0, i], tck)
    
    # Make array for angle vs flux
    angvflux  = np.zeros((2, kurucznew.shape[0] - 1))
    angvflux[0] = np.loadtxt(kuruczfile, skiprows=2, usecols=range(0,17))[0]
    
    # Integrate spitzer / kurucz data over wavelength
    intbottom = integrate(spitfit[0], spitfit[1])
    #print("Assuming kurucznew is actual spectrum (not filter) that needs to be normalized. See Line 130.")
    for i in range(kurucznew.shape[0] - 1):
        if isfilter == True:
            #Use filter
            inttop    = integrate(spitfit[0], spitfit[1] * kurucznew[i + 1])
        else:
            #Use actual spectrum
            inttop    = integrate(spitfit[0], spitfit[1] * kurucznew[i + 1]/kurucznew[0])
        #inttop    = integrate(kurucznew[0] * spitfit[0], spitfit[1] * kurucznew[i + 1])
        angvflux[1,i] = inttop / intbottom
    
    # I(mu) / I(1) = 1 - a1(1 - mu^[1/2]) - a2(1 - mu) 
    #                  - a3(1 - mu^[3/2]) - a4(1 - mu^2)
    # Therefore:
    # 1 - I(mu) / I(1) =   a1(1 - mu^[1/2]) + a2(1 - mu) 
    #                    + a3(1 - mu^[3/2]) + a4(1 - mu^2)
    
    
    # Use least sqaures to find parameters (a1, etc.) for best 
    # model fit using the above equation
    #mu = angvflux[0]
    #y  = 1. - angvflux[1] / angvflux[1,0]
    #Trim mu < 0.1
    imu     = np.where(angvflux[0] >= 0.1)[0]
    mu      = angvflux[0,imu]
    y       = 1. - angvflux[1,imu] / angvflux[1,0]
    
    A = np.zeros((len(mu), Nparam))
    if Nparam == 4:
        A[:, 0] = (1 - mu**0.5)
        A[:, 1] = (1 - mu     )
        A[:, 2] = (1 - mu**1.5)
        A[:, 3] = (1 - mu**2. )
    elif Nparam == 2:
        A[:, 0] = (1. - mu)
        A[:, 1] = (1. - mu)**2
    elif Nparam == 1:
        A[:, 0] = (1. - mu)
    
    coeff = lstsq(A, y)[0]
        
    # If wanted, plot results using 
    '''
    plt.ion()
    plt.plot(angvflux[0], angvflux[1], '.')
    plt.plot(angvflux[0], (1 - coeff[0] * (1 - angvflux[0]**0.5) - coeff[1] \
         * (1 - angvflux[0]) - coeff[2] * (1 - angvflux[0]**1.5) - coeff[3] \
         * (1 - angvflux[0]**2.)) * angvflux[1,0])
    plt.ylabel("Intensity in $ergs * cm^-2 * s^-1 * hz ^-1 * ster^-1$ ")
    
    plt.xlabel('cosine of Angle (center = 1)')
    plt.title('Angle vs Intensity') # '''
    if returnmu:
        return (coeff, mu, 1-y)
    else:
        return(coeff)
