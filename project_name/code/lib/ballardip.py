
import numpy as np

def ballardip (params, pos, etc):
    sigmay, sigmax, nbins = params
    y, x, q = position
    flux = etc[0]
    nobj = y.size
    
    weight  = np.ones(nobj)
    for i in range(nbins):
        start   = int(1.*i*nobj/nbins)
        end     = int(1.*(i+1)*nobj/nbins)
        s       = np.ones(nobj)
        s[start:end] = 0
        biny = np.mean(y[start:end])
        binx = np.mean(x[start:end])
        weight[start:end] = sum(np.exp(-0.5*((x-binx)/sigmax)**2) * \
                                np.exp(-0.5*((y-biny)/sigmay)**2)*flux*s) / \
                            sum(np.exp(-0.5*((x-binx)/sigmax)**2) * \
                                np.exp(-0.5*((y-biny)/sigmay)**2)*s)
    
    return weight
