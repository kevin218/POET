import numpy as np
import models
from scipy.optimize import leastsq

# chisquared - Rising exp
def chisqre(p, x, y, sigma):
    return (models.risingexp(p, x) - y)/sigma

# chisquared - Rising exp + linear
def chisqrel(p, x, y, sigma):
    return (models.relramp(p, x) - y)/sigma

# chisquared - Rising exp + quadratic
def chisqreq(p, x, y, sigma):
    return (models.reqramp(p, x) - y)/sigma

# chisquared - Rising exp + rising exp
def chisqre2(p, x, y, sigma):
    return (models.re2ramp(p, x) - y)/sigma

def corrfit(guess, x, y, yerr):
    args = (x, y, yerr)
    # fit y as a function of x
    fitpar, err = leastsq(chisqre, guess, args=args, factor=100, ftol=1e-16, xtol=1e-16, gtol=1e-16)
    # return fitted parameters
    return fitpar
