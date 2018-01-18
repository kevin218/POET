import numpy as np
import models
from scipy.optimize import leastsq

# chisquared
def chisq(p, x, y, sigma):
    return (models.risingexp(p, x) - y)/sigma

def corrfit(guess, x, y, yerr):
    """
    guess is an array in the format:
    guess[0] = goal
    guess[1] = curvature
    guess[2] = offset

    These are parameters of the rising exponential function,
    f(x) = goal*(1 - exp(-curvature*(x - offset)))
    """
    args = (x, y, yerr)
    # fit y as a function of x
    fitpar, err = leastsq(chisq, guess, args=args, factor=100, ftol=1e-16, xtol=1e-16, gtol=1e-16)
    # return fitted parameters
    return fitpar
