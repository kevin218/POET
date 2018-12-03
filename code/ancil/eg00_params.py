
import numpy as np
import sys

# p is an event.params object
def modelparams(p):
    p.planetname = 'Example-0b'
    
    #Name of file(s) containing initial parameter values
    p.modelfile  = ['egOO-initvals.txt']

    #SPECIFY OUTPUT FILE FOR print STATEMENTS
    #You can no longer use sys.stdout for standard output (to screen)
    p.printout   =  "results.txt"
    
    #CHOOSE 'orbits', 'days-utc' OR 'days-tdb' FOR MODELING AND PLOTING TIME UNIT
    p.timeunit = 'orbits'
    p.tuoffset = 0.0

    #Model names and descriptions listed in /home/esp01/doc/models.txt.
    p.model = [['mandelecl', 'linramp', 'bilinint'],
               ['ortho', 'mandelecl', 'seramp', 'nnint']]
    
    #Model Constants
    p.numit      = np.array([1e3, 1e4])       #Number of iterations, [burn-in,final]
    p.nbins      = 60                         #Number of bins
    p.nchains    = 4                          #Number of chains in Gelman-Rubin test
    p.newortho   = False                      #Set True to orthogonalize parameters in ortholist
    p.ortholist  = [['flux', 'ser0', 'ser1'],
                    ['ser0', 'ser1']]         #List(s) of correlated parameters
    p.chi2flag   = 1                          #1: Adjust sigma s.t. reduced Chi^2 = 1
    p.numcalc    = 50000                      #Number of temperature calculations
    p.stepsize   = 100                        #Stepsize for histograms and correlation plots
    #p.rmsbins    = 5000                      #Max number of bins in RMS calc, default is 0.1*length(data)
    p.allplots   = 5                          #Disply no plots (0), some plots (1-4), or all plots (5)
    p.leastsq    = True                       #Perform least-squares fit before MCMC
    p.normflux   = True                       #Normalize flux at each position (set False for posflux)
    p.noisysdnr  = None                       #SDNR of noisy data set (None unless using denoised data)
    p.noisewavelet = ['None']                 #List of wavelet names (e.g. 'db19', 'haar', etc)
    p.savedata   = False                      #Set False to skip saving data
    p.isresidperm  = False                    #Set True to implement residual permutation technicque for estimating uncertainties
    p.night      = 0                          #Set to 0 unless using Divide White technique
    
    p.preclip    = [0]                        #Remove first 'preclip' points
    p.postclip   = [0]                        #Remove last 'postclip' points
    #p.interclip  = [[lower,upper]]           #Remove points within specified range from model fits,
                                              #Can handle multiple ranges
    
    #Priors
    #p.priorvars  = ['cosi', 'ars']           #List of variables with priors, applies to each model
    #p.priorvals  = np.array(([[0.0849, 0.002, 0.002],[9.1027, 0.067, 0.060]]))  
                                              #Include [mode, lower sigma, upper sigma] to describe
                                              #each asymmetric Gaussian prior
    
    #Bin sizes for interpolative ip models
    p.xstep      = [0.01]                     #Bin size in x
    p.ystep      = [0.01]                     #Bin size in y
    p.minnumpts  = [1]                        #Minimum # of points in a bin
    p.isfixipmap = False                      #Fix intrapixel mapping to best-fit values after burn-in
    p.sssteps    = [1,2,3,4,6,8,10,30,100,300,1000]
    #p.ipclip     = [[lower,upper]]           #Remove points within specified range from intrapixel mapping
    
    #Smoothing parameters
    p.issmoothing= False                     #Turn smoothing on/off (True/False)
    p.nx         = [3.]                      #Kernel halfwidth in x
    p.ny         = [3.]                      #Kernel halfwidth in y
    p.sx         = [1.5]                     #Gaussian kernel std dev in x
    p.sy         = [1.5]                     #Gaussian kernel std dev in y
    
    #Used with quadip4, otherwise set to 0
    p.xdiv       = 0.0                       #Divide into quadrants here along x axis
    p.ydiv       = 0.0                       #Divide into quadrants here along y axis
    return
