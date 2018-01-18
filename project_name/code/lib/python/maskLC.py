
import smooth
import numpy as np

def boxCarMedianMask(data, window_len, maxstd, mask=None, unct=None):
    '''
    Mask bad datapoints using a box car median filter
    
    INPUTS
    ------
    data        : Light curve array
    window_len  : Length of box car filter
    maxstd      : Maximum number of standard deviations above which poitns are masked
    mask        : Mask array (optional)
    unct        : Uncertainty array (optional)
    
    OUTPUT
    ------
    mask
    
    HISTORY
    -------
    Written by Kevin Stevenson      June 2015
    '''
    npts    = len(data)
    if window_len % 2 == 0:
        print("Median filter length ("+str(window_len)+") must be odd. Adding 1.")
        window_len += 1
    if mask == None:
        mask = np.ones(npts)
    #Median smooth light curve
    smdata      = smooth.medfilt(data, window_len)
    residuals   = data - smdata
    #Estimate typical uncertainty
    if unct == None:
        unct = np.std(residuals[np.where(mask)])
    #Count number of standard deviations from median
    stdev   = np.abs(residuals/unct)
    #Mask bad points
    mask[np.where(stdev > maxstd)] = 0
    
    return mask
