import numpy as np
import time, timer
import gelmanrubin as gr

def mcmc(y, x, params, pmin, pmax, stepsize, numit, sigma, numparams, nummodels, myfuncs, funcx, numaccept):
   """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, x, params, pmin, pmax stepsize, numit, sigma, numparams, nummodels, myfuncs, funcx, numaccept)

 INPUTS:
	y:	   Array containing dependent data
	x:	   Array containing independent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Nummodels: Number of models used
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	Numaccept: Number of accepted steps

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-05-02
			kbstevenson@gmail.com
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
			kbstevenson@gmail.com
   """

   complete   = 0
   nobj       = len(y)
   nump       = len(params)
   nextp      = np.copy(params)
   bestp      = np.copy(params)
   allparams  = np.zeros((nump, numit))
   inotfixed  = np.where(stepsize > 0)
   iequal     = np.array(np.where(stepsize < 0))
   outside    = np.zeros(nump)
   
   #Calc chi-squared for model type using current params
   ymodel     = np.ones(nobj)
   for i in range(nummodels):
      ymodel *= myfuncs[i](params[numparams[i]:numparams[i+1]], funcx[i])
   
   currchisq  = sum((ymodel - y)**2 / sigma**2)
   bestchisq  = currchisq
   
   #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
   for j in range(0, numit):
      	#Take step in random direction for adjustable parameters
        nextp[inotfixed] = np.random.normal(params[inotfixed],stepsize[inotfixed])
      	#CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside         = np.array(np.where((nextp < pmin) | (nextp > pmax)))
        if (ioutside.size > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
	#UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]]))]
	#COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        ymodel     = np.ones(nobj)
        for i in range(nummodels):
            ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i])
	
        nextchisq  = sum((ymodel - y)**2 / sigma**2)
        accept = np.exp(0.5 * (currchisq - nextchisq))

        if (accept >= 1) or (np.random.uniform(0, 1) <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        if (j % (numit / 5) == 0): 
            print(str(complete * 20) + "% complete at " + time.ctime())
            print("Number of times new step is outside its boundary:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            complete  += 1

   return allparams, bestp, numaccept

def mcmc2(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, myfuncs, funcx, fit):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, nummodels, myfuncs, funcx, numaccept)

 INPUTS:
	y:	   Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	#Numaccept: Number of accepted steps

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
			kevin218@knights.ucf.edu
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
			kevin218@knights.ucf.edu
    """
    
    numaccept  = 0
    complete   = 0
    nump       = len(params)
    nextp      = np.copy(params)
    bestp      = np.copy(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize > 0)
    iequal     = np.array(np.where(stepsize < 0))
    outside    = np.zeros(nump)
    numevents  = len(fit)
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    #Calc chi-squared for model type using current params
    currchisq  = 0
    for j in range(numevents):
        ymodel     = np.ones(fit[j].nobj)
        for i in range(cummodels[j],cummodels[j+1]):
            ymodel *= myfuncs[i](params[numparams[i]:numparams[i+1]], funcx[i], ymodel)
        
        currchisq  += sum((ymodel - y[j])**2 / sigma[j]**2)
    
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numnotfixed = len(inotfixed[0])
    step = np.random.normal(0, stepsize[inotfixed], (numit,numnotfixed))
    unif = np.random.uniform(0, 1, numit)

    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[inotfixed] = params[inotfixed] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside         = np.array(np.where((nextp < pmin) | (nextp > pmax)))
        if (ioutside.size > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        nextchisq   = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            for i in range(cummodels[k],cummodels[k+1]):
                ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], ymodel)
            
            nextchisq  += sum((ymodel - y[k])**2 / sigma[k]**2)
        
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        if (j % (numit / 5) == 0): 
            print(str(complete * 20) + "% complete at " + time.ctime())
            print("Number of times new step is outside its boundary:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            complete  += 1
    
    return allparams, bestp, numaccept

def mcmc3(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit)

 INPUTS:
	y:	   Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
    Functype:  Define function type (eclipse, ramp, ip, etc), see models.py
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	#Numaccept: Number of accepted steps

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
			kevin218@knights.ucf.edu
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
			kevin218@knights.ucf.edu
    Updated for ipspline, nnint & bilinint:	
			Kevin Stevenson, UCF  	2010-06-09
			kevin218@knights.ucf.edu
    """
    
    numaccept  = 0
    complete   = 0
    nump       = len(params)
    nextp      = np.copy(params)
    bestp      = np.copy(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize > 0)
    iequal     = np.array(np.where(stepsize < 0))
    outside    = np.zeros(nump)
    numevents  = len(fit)
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    #Calc chi-squared for model type using current params
    currchisq  = 0
    for j in range(numevents):
        ymodel     = np.ones(fit[j].nobj)
        m          = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'ipmap':
                fit[j].etc[m] = ymodel
            ymodel *= myfuncs[i](params[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            m      += 1
        
        currchisq  += sum((ymodel - y[j])**2 / sigma[j]**2)
    
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numnotfixed = len(inotfixed[0])
    step = np.random.normal(0, stepsize[inotfixed], (numit,numnotfixed))
    unif = np.random.uniform(0, 1, numit)

    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[inotfixed] = params[inotfixed] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside         = np.array(np.where((nextp < pmin) | (nextp > pmax)))
        if (ioutside.size > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        nextchisq   = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            m          = 0
            for i in range(cummodels[k],cummodels[k+1]):
                if functype[i] == 'ipmap':
                    fit[k].etc[m] = ymodel
                ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                m      += 1
            
            nextchisq  += sum((ymodel - y[k])**2 / sigma[k]**2)
        
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        if (j % (numit / 5) == 0): 
            print(str(complete * 20) + "% complete at " + time.ctime())
            print("Number of times new step is outside its boundary:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            complete  += 1
    
    return allparams, bestp, numaccept


def mcmc4(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit)

 INPUTS:
	y:	   Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
    Functype:  Define function type (eclipse, ramp, ip, etc), see models.py
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	#Numaccept: Number of accepted steps

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
			kevin218@knights.ucf.edu
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
			kevin218@knights.ucf.edu
    Updated for ipspline, nnint & bilinint:	
			Kevin Stevenson, UCF  	2010-06-09
			kevin218@knights.ucf.edu
    """
    
    numaccept  = 0
    complete   = 0
    nump       = len(params)
    nextp      = np.copy(params)
    bestp      = np.copy(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize > 0)
    iequal     = np.array(np.where(stepsize < 0))
    outside    = np.zeros(nump)
    numevents  = len(fit)
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    #Calc chi-squared for model type using current params
    currchisq  = 0
    for j in range(numevents):
        ymodel          = np.ones(fit[j].nobj)
        m               = 0
        #allknots = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'ipmap':
                #fit[j].allknots = np.zeros((numit, fit[j].xygrid[0].shape[0]*fit[j].xygrid[0].shape[1]))
                fit[j].etc[m]   = ymodel
            ymodel *= myfuncs[i](params[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            m      += 1
        
        currchisq  += sum(((ymodel - y[j]) / sigma[j])**2)
    
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numnotfixed = len(inotfixed[0])
    step = np.random.normal(0, stepsize[inotfixed], (numit,numnotfixed))
    unif = np.random.uniform(0, 1, numit)

    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[inotfixed] = params[inotfixed] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside         = np.array(np.where((nextp < pmin) | (nextp > pmax)))
        if (ioutside.size > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        nextchisq   = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            m          = 0
            for i in range(cummodels[k],cummodels[k+1]):
                if functype[i] == 'ipmap':
                    ipmodel, allknots = myfuncs[i](params[numparams[i]:numparams[i+1]], funcx[i], \
                                                      ymodel, retbinflux=True)
                    ymodel              *= ipmodel
                    #COMPRESS allknots BY SUBTRACTING 1, MULTIPLYING BY 30000 THEN CONVERTING TO INT16
                    np.save(fit[k].allknotpid, np.int16((allknots-1.)*30000))
                else:
                    ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                m      += 1
            
            nextchisq  += sum(((ymodel - y[k]) / sigma[k])**2)
            del(ymodel, m)
        
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        if (j % (numit / 5) == 0): 
            print(str(complete * 20) + "% complete at " + time.ctime())
            print("Number of times new step is outside its boundary:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            complete  += 1
            for k in range(numevents):
                for i in range(cummodels[k],cummodels[k+1]):
                    if functype[i] == 'ipmap':
                        fit[k].allknotpid.flush()
    
    return allparams, bestp, numaccept

def mcmc5(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit)

 INPUTS:
	y:	   Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
    Functype:  Define function type (eclipse, ramp, ip, etc), see models.py
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	fit:       List of fit objects

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
			kevin218@knights.ucf.edu
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
			kevin218@knights.ucf.edu
    Updated for ipspline, nnint & bilinint:	
			Kevin Stevenson, UCF  	2010-06-09
			kevin218@knights.ucf.edu
    """
    import models_c
    numaccept  = 0
    complete   = 0
    nump       = len(params)
    nextp      = np.copy(params)
    bestp      = np.copy(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize > 0)
    iequal     = np.array(np.where(stepsize < 0))
    outside    = np.zeros(nump)
    numevents  = len(fit)
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    #Calc chi-squared for model type using current params
    currchisq  = 0
    for j in range(numevents):
        ymodel          = np.ones(fit[j].nobj)
        m               = 0
        #allknots = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'ipmap':
                #fit[j].allknots = np.zeros((numit, fit[j].xygrid[0].shape[0]*fit[j].xygrid[0].shape[1]))
                fit[j].etc[m]   = ymodel
            ymodel *= myfuncs[i](params[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            m      += 1
        
        currchisq  += models_c.chisq(ymodel, y[j], sigma[j])
        #currchisq  += sum(((ymodel - y[j]) / sigma[j])**2)
    
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numnotfixed = len(inotfixed[0])
    step = np.random.normal(0, stepsize[inotfixed], (numit,numnotfixed))
    unif = np.random.uniform(0, 1, numit)

    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[inotfixed] = params[inotfixed] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside         = np.array(np.where((nextp < pmin) | (nextp > pmax)))
        if (ioutside.size > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        nextchisq   = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            m          = 0
            for i in range(cummodels[k],cummodels[k+1]):
                if functype[i] == 'ipmap':
                    ipmodel, allknots = myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], \
                                                      ymodel, retbinflux=True)
                    ymodel              *= ipmodel
                    #COMPRESS allknots BY SUBTRACTING 1, MULTIPLYING BY 300000 THEN CONVERTING TO INT16
                    np.save(fit[k].allknotpid, np.int16((allknots[np.where(allknots > 0)]-1.)*300000))
                else:
                    ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                m      += 1
            
            nextchisq  += models_c.chisq(ymodel, y[k], sigma[k])
            #nextchisq  += sum(((ymodel - y[k]) / sigma[k])**2)
            del(ymodel, m)
        
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        if (j % (numit / 5) == 0): 
            print(str(complete * 20) + "% complete at " + time.ctime())
            print("Number of times new step is outside its boundary:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            complete  += 1
            for k in range(numevents):
                for i in range(cummodels[k],cummodels[k+1]):
                    if functype[i] == 'ipmap':
                        fit[k].allknotpid.flush()
    
    return allparams, bestp, numaccept

def mcmc6(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit, nchains=4):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit)

 INPUTS:
	y:	   Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
    Functype:  Define function type (eclipse, ramp, ip, etc), see models.py
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	fit:       List of fit objects

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
			kevin218@knights.ucf.edu
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
			kevin218@knights.ucf.edu
    Updated for ipspline, nnint & bilinint:	
			Kevin Stevenson, UCF  	2010-06-09
			kevin218@knights.ucf.edu
    Updated for Gelman-Rubin statistic:	
			Kevin Stevenson, UCF  	2011-07-06
			kevin218@knights.ucf.edu
    """
    import models_c
    numaccept  = 0
    nump       = len(params)
    nextp      = np.copy(params)
    bestp      = np.copy(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize != 0)
    iequal     = np.array(np.where(stepsize < 0))
    ifree      = np.where(stepsize > 0)
    outside    = np.zeros(nump)
    numevents  = len(fit)
    intsteps   = np.min((numit/5,1e5))
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    
    #Calc chi-squared for model type using current params
    currchisq  = 0
    for j in range(numevents):
        ymodel          = np.ones(fit[j].nobj)
        m               = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if functype[i] == 'ipmap':
                fit[j].etc[m]   = ymodel
            ymodel *= myfuncs[i](params[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            m      += 1
        
        currchisq  += models_c.chisq(ymodel, y[j], sigma[j])
    
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numfree     = len(ifree[0])
    numnotfixed = len(inotfixed[0])
    step = np.random.normal(0, stepsize[ifree], (numit,numfree))
    unif = np.random.uniform(0, 1, numit)
    
    #START TIMER
    clock = timer.Timer(numit,progress = np.arange(0,1.01,0.05))
    
    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[ifree] = params[ifree] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside         = np.array(np.where((nextp < pmin) | (nextp > pmax)))
        if (ioutside.size > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        nextchisq   = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            m          = 0
            for i in range(cummodels[k],cummodels[k+1]):
                if functype[i] == 'ipmap':
                    ipmodel, allknots = myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], \
                                                      ymodel, retbinflux=True)
                    ymodel              *= ipmodel
                    #COMPRESS allknots BY SUBTRACTING 1, MULTIPLYING BY 300000 THEN CONVERTING TO INT16
                    np.save(fit[k].allknotpid, np.int16((allknots[np.where(allknots > 0)]-1.)*300000))
                else:
                    ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                m      += 1
            
            nextchisq  += models_c.chisq(ymodel, y[k], sigma[k])
            #nextchisq  += sum(((ymodel - y[k]) / sigma[k])**2)
            del(ymodel, m)
        
        #CALCULATE ACCEPTANCE PROBABILITY
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        #PRINT INTERMEDIATE INFO
        if ((j+1) % intsteps == 0) and (j > 0): 
            print("\n" + time.ctime())
            print("Number of times parameter tries to step outside its prior:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            #Flush allknots to file
            for k in range(numevents):
                for i in range(cummodels[k],cummodels[k+1]):
                    if functype[i] == 'ipmap':
                        fit[k].allknotpid.flush()
            
            #Apply Gelman-Rubin statistic
            psrf, meanpsrf = gr.convergetest(allparams[inotfixed[0],0:j+1],nchains)
            numconv = np.sum(psrf < 1.01)
            print("Gelman-Rubin statistic for free parameters:")
            print(psrf)
            #print(str(numconv) +"/"+ str(numnotfixed) +" parameters have converged.")
            if numconv == numnotfixed:
                print("All parameters have converged to within 1% of unity. Halting MCMC.")
                allparams = allparams[:,0:j+1]
                break
            #else:
                #print("Worst Gelman-Rubin statistic: " + str(np.max(psrf)))
        
        clock.check(j+1)
    
    return allparams, bestp, numaccept, j+1

def mcmc7(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, fit, nchains=4):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit)

 INPUTS:
	y:	   Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
    Functype:  Define function type (eclipse, ramp, ip, etc), see models.py
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	fit:       List of fit objects

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 MODIFICATION HISTORY:
 	Written by:	
 	        Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
    Updated for ipspline, nnint & bilinint:	
			Kevin Stevenson, UCF  	2010-06-09
    Updated for Gelman-Rubin statistic:	
			Kevin Stevenson, UCF  	2011-07-06
    Added principal component analysis:	
			Kevin Stevenson, UCF  	2011-07-22
    """
    import models_c as mc
    numaccept  = 0
    nump       = len(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize != 0)
    iequal     = np.array(np.where(stepsize < 0))
    ifree      = np.where(stepsize > 0)
    outside    = np.zeros(nump)
    numevents  = len(fit)
    intsteps   = np.min((numit/5,1e5))
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    
    #Calc chi-squared for model type using current params
    nextp      = np.copy(params)    #Proposed parameters
    bestp      = np.copy(params)    #Best-fit parameters
    pedit      = np.copy(params)    #Editable parameters
    currchisq  = 0
    for j in range(numevents):
        ymodel          = np.ones(fit[j].nobj)
        m               = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pedit[iortholist[j]] = myfuncs[i](pedit[iortholist[j]], funcx[i], fit[j].etc[m])
            elif functype[i] == 'ipmap':
                fit[j].etc[m]   = ymodel
                ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            elif functype[i] == 'posoffset':
                # Record change in Position 0 => cannot orthogonalize position parameters
                ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            else:
                ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            m      += 1
        
        currchisq  += mc.chisq(ymodel, y[j], sigma[j])
    
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numfree     = len(ifree[0])
    numnotfixed = len(inotfixed[0])
    step = np.random.normal(0, stepsize[ifree], (numit,numfree))
    unif = np.random.uniform(0, 1, numit)
    
    #START TIMER
    clock = timer.Timer(numit,progress = np.arange(0.05,1.01,0.05))
    
    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[ifree] = params[ifree] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside     = np.where(np.bitwise_or(nextp < pmin, nextp > pmax))[0]
        if (len(ioutside) > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        pedit        = np.copy(nextp)
        nextchisq    = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            m          = 0
            for i in range(cummodels[k],cummodels[k+1]):
                if   functype[i] == 'ortho':
                    #MODIFY COPY OF nextp ONLY
                    pedit[iortholist[k]] = myfuncs[i](pedit[iortholist[k]], funcx[i], fit[k].etc[m])
                elif functype[i] == 'ipmap':
                    ipmodel, allknots = myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], \
                                                      ymodel, retbinflux=True)
                    ymodel *= ipmodel
                    #COMPRESS allknots BY SUBTRACTING 1, MULTIPLYING BY 300000 THEN CONVERTING TO INT16
                    np.save(fit[k].allknotpid, np.int16((allknots[np.where(allknots > 0)]-1.)*300000))
                elif functype[i] == 'posoffset':
                    # Record change in Position 0 => cannot orthogonalize position parameters
                    ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                else:
                    ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                m      += 1
            
            nextchisq  += mc.chisq(ymodel, y[k], sigma[k])
            del(ymodel, m)
        
        #CALCULATE ACCEPTANCE PROBABILITY
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        #PRINT INTERMEDIATE INFO
        if ((j+1) % intsteps == 0) and (j > 0): 
            print("\n" + time.ctime())
            print("Number of times parameter tries to step outside its prior:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            #Flush allknots to file
            for k in range(numevents):
                for i in range(cummodels[k],cummodels[k+1]):
                    if functype[i] == 'ipmap':
                        fit[k].allknotpid.flush()
            
            #Apply Gelman-Rubin statistic
            if nchains > 0:
                psrf, meanpsrf = gr.convergetest(allparams[inotfixed[0],0:j+1],nchains)
                numconv = np.sum(psrf < 1.01)
                print("Gelman-Rubin statistic for free parameters:")
                print(psrf)
                if numconv == numnotfixed:
                    print("All parameters have converged to within 1% of unity. Halting MCMC.")
                    allparams = allparams[:,0:j+1]
                    break
            
        
        clock.check(j+1)
    
    return allparams, bestp, numaccept, j+1

def mcmc8(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, fit, nchains=4):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit)

 INPUTS:
	y:	   Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
    Functype:  Define function type (eclipse, ramp, ip, etc), see models.py
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	fit:       List of fit objects

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 MODIFICATION HISTORY:
 	Written by:	
 	        Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
    Updated for ipspline, nnint & bilinint:	
			Kevin Stevenson, UCF  	2010-06-09
    Updated for Gelman-Rubin statistic:	
			Kevin Stevenson, UCF  	2011-07-06
    Added principal component analysis:	
			Kevin Stevenson, UCF  	2011-07-22
	Added priors:
	        Kevin Stevenson, UCF  	2011-10-11
    """
    import models_c as mc
    numaccept  = 0
    nump       = len(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize != 0)[0]
    iequal     = np.array(np.where(stepsize < 0))
    ifree      = np.where(stepsize > 0)[0]
    outside    = np.zeros(nump)
    numevents  = len(fit)
    intsteps   = np.min((numit/5,1e5))
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    
    #Calc chi-squared for model type using current params
    nextp      = np.copy(params)    #Proposed parameters
    bestp      = np.copy(params)    #Best-fit parameters
    pedit      = np.copy(params)    #Editable parameters
    currchisq  = 0
    for j in range(numevents):
        ymodel          = np.ones(fit[j].nobj)
        m               = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pedit[iortholist[j]] = myfuncs[i](pedit[iortholist[j]], funcx[i], fit[j].etc[m])
            elif functype[i] == 'ipmap':
                fit[j].etc[m]   = ymodel
                ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            elif functype[i] == 'posoffset':
                # Record change in Position 0 => cannot orthogonalize position parameters
                ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            else:
                ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            m      += 1
        # Calculate chi^2
        currchisq  += mc.chisq(ymodel, y[j], sigma[j])
        # Apply prior, if one exists
        if len(fit[j].ipriors) > 0:
        #for i in range(len(fit[j].ipriors)):
            pbar   = fit[j].priorvals[:,0]  #prior mean
            psigma = np.zeros(len(pbar))    #prior standard deviation
            # Determine psigma based on which side of asymmetric Gaussian nextp is on
            for i in range(len(fit[j].ipriors)):
                if fit[j].ipriors[i] < pbar[i]:
                    psigma[i] = fit[j].priorvals[i,1]
                else:
                    psigma[i] = fit[j].priorvals[i,2]
                currchisq += ((nextp[fit[j].ipriors[i]] - pbar[i])/psigma[i])**2
        
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numfree     = len(ifree)
    numnotfixed = len(inotfixed)
    step = np.random.normal(0, stepsize[ifree], (numit,numfree))
    unif = np.random.uniform(0, 1, numit)
    
    #START TIMER
    clock = timer.Timer(numit,progress = np.arange(0.05,1.01,0.05))
    
    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[ifree] = params[ifree] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside     = np.where(np.bitwise_or(nextp < pmin, nextp > pmax))[0]
        if (len(ioutside) > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        pedit        = np.copy(nextp)
        nextchisq    = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            m          = 0
            for i in range(cummodels[k],cummodels[k+1]):
                if   functype[i] == 'ortho':
                    #MODIFY COPY OF nextp ONLY
                    pedit[iortholist[k]] = myfuncs[i](pedit[iortholist[k]], funcx[i], fit[k].etc[m])
                elif functype[i] == 'ipmap':
                    ipmodel, allknots = myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], \
                                                      ymodel, retbinflux=True)
                    ymodel *= ipmodel
                    #COMPRESS allknots BY SUBTRACTING 1, MULTIPLYING BY 300000 THEN CONVERTING TO INT16
                    np.save(fit[k].allknotpid, np.int16((allknots[np.where(allknots > 0)]-1.)*300000))
                elif functype[i] == 'posoffset':
                    # Record change in Position 0 => cannot orthogonalize position parameters
                    ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                else:
                    ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                m      += 1
            # Calculate chi^2
            nextchisq  += mc.chisq(ymodel, y[k], sigma[k])
            # Apply prior, if one exists
            if len(fit[k].ipriors) > 0:
                pbar   = fit[k].priorvals[:,0]  #prior mean
                psigma = np.zeros(len(pbar))    #prior standard deviation
                # Determine psigma based on which side of asymmetric Gaussian nextp is on
                for i in range(len(fit[k].ipriors)):
                    if fit[k].ipriors[i] < pbar[i]:
                        psigma[i] = fit[k].priorvals[i,1]
                    else:
                        psigma[i] = fit[k].priorvals[i,2]
                    nextchisq += ((nextp[fit[k].ipriors[i]] - pbar[i])/psigma[i])**2
                del(pbar,psigma)
            del(ymodel, m)
        
        #CALCULATE ACCEPTANCE PROBABILITY
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        #PRINT INTERMEDIATE INFO
        if ((j+1) % intsteps == 0) and (j > 0): 
            print("\n" + time.ctime())
            print("Number of times parameter tries to step outside its prior:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            #Flush allknots to file
            for k in range(numevents):
                for i in range(cummodels[k],cummodels[k+1]):
                    if functype[i] == 'ipmap':
                        fit[k].allknotpid.flush()
            
            #Apply Gelman-Rubin statistic
            if nchains > 0:
                psrf, meanpsrf = gr.convergetest(allparams[inotfixed,0:j+1],nchains)
                numconv = np.sum(psrf < 1.01)
                print("Gelman-Rubin statistic for free parameters:")
                print(psrf)
                if numconv == numnotfixed:
                    print("All parameters have converged to within 1% of unity. Halting MCMC.")
                    allparams = allparams[:,0:j+1]
                    break
            
        
        clock.check(j+1)
    
    return allparams, bestp, numaccept, j+1


def mcmc9(y, params, pmin, pmax, stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, iortholist, fit, nchains=4):
    """
 NAME:
	MCMC

 PURPOSE:
	This function runs Markov Chain Monte Carlo model fitting using the Metropolis-Hastings algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	allparams, bestp, numaccept = MCMC(y, params, pmin, pmax stepsize, numit, sigma, numparams, cummodels, functype, myfuncs, funcx, fit)

 INPUTS:
	y:         Array containing dependent data
	Params:    Array of initial guess for parameters
	Pmin:      Array of parameter minimum values
	Pmax:      Array of parameter maximum values
	stepsize:  Array of 1-sigma change in parameter per iteration
	Numit:	   Number of iterations to perform
	Sigma:	   Standard deviation of data noise in y
	Numparams: Number of parameters for each model
	Cummodels: Cumulative number of models used
    Functype:  Define function type (eclipse, ramp, ip, etc), see models.py
	Myfuncs:   Pointers to model functions
	Funcx:	   Array of x-axis values for myfuncs
	fit:       List of fit objects
	nchains:   Number of chains for G-R statistic

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 MODIFICATION HISTORY:
 	Written by:	
 	        Kevin Stevenson, UCF  	2008-05-02
			kevin218@knights.ucf.edu
 	Finished updating:	
			Kevin Stevenson, UCF  	2008-06-21
    Updated for multi events:	
			Kevin Stevenson, UCF  	2009-11-01
    Updated for ipspline, nnint & bilinint:	
			Kevin Stevenson, UCF  	2010-06-09
    Updated for Gelman-Rubin statistic:	
			Kevin Stevenson, UCF  	2011-07-06
    Added principal component analysis:	
			Kevin Stevenson, UCF  	2011-07-22
	Added priors:
	        Kevin Stevenson, UCF  	2011-10-11
	Added increased resolution models for long exposure times:
	        Kevin Stevenson, UCF  	2012-07-23
    """
    import models_c as mc
    numaccept  = 0
    nump       = len(params)
    allparams  = np.zeros((nump, numit))
    inotfixed  = np.where(stepsize != 0)[0]
    iequal     = np.array(np.where(stepsize < 0))
    ifree      = np.where(stepsize > 0)[0]
    outside    = np.zeros(nump)
    numevents  = len(fit)
    intsteps   = np.min((numit/5,1e5))
    isrednoise = False
    
    #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
    if (iequal.size > 0):
        for i in range(iequal.size):
            params[iequal[0][i]] = params[int(abs(stepsize[iequal[0][i]])-1)]
    
    #Calc chi-squared for model type using current params
    nextp      = np.copy(params)    #Proposed parameters
    bestp      = np.copy(params)    #Best-fit parameters
    pedit      = np.copy(params)    #Editable parameters
    currchisq  = 0
    for j in range(numevents):
        ymodel          = np.ones(fit[j].nobj)
        m               = 0
        for i in range(cummodels[j],cummodels[j+1]):
            if   functype[i] == 'ortho':
                pedit[iortholist[j]] = myfuncs[i](pedit[iortholist[j]], funcx[i], fit[j].etc[m])
            elif (functype[i] == 'ipmap') or (functype[i] == 'spline'):
                fit[j].etc[m]   = ymodel
                ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            elif functype[i] == 'posoffset':
                # Record change in Position 0 => cannot orthogonalize position parameters
                ymodel *= myfuncs[i](nextp[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            elif hasattr(fit[j], 'timebins') and (functype[i] == 'ecl/tr' 
                                              or  functype[i] == 'ramp' 
                                              or  functype[i] == 'sinusoidal'):
                # Average over high-resolution model
                hiresmodel = myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
                if len(fit[j].timebins) == fit[j].nobj:
                    for k in range(len(fit[j].timebins)):
                        ymodel[k] *= np.mean(hiresmodel[fit[j].timebins[k]])
                else:
                    for k in range(len(fit[j].timebinsuc)):
                        ymodel[k] *= np.mean(hiresmodel[fit[j].timebinsuc[k]])
            elif functype[i] == 'noise':
                isrednoise   = True
                wavelet      = fit[j].etc[m]
                noisepars    = pedit[numparams[i]:numparams[i+1]]
                noisefunc    = myfuncs[i]
            else:
                ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[j].etc[m])
            m      += 1
        # Calculate chi^2
        if isrednoise == False:
            currchisq  += mc.chisq(ymodel, y[j], sigma[j])
        else:
            currchisq  += noisefunc(noisepars, ymodel-y[j], wavelet)
        # Apply prior, if one exists
        if len(fit[j].ipriors) > 0:
        #for i in range(len(fit[j].ipriors)):
            pbar   = fit[j].priorvals[:,0]  #prior mean
            psigma = np.zeros(len(pbar))    #prior standard deviation
            # Determine psigma based on which side of asymmetric Gaussian nextp is on
            for i in range(len(fit[j].ipriors)):
                if nextp[fit[j].ipriors[i]] < pbar[i]:
                    psigma[i] = fit[j].priorvals[i,1]
                else:
                    psigma[i] = fit[j].priorvals[i,2]
                currchisq += ((nextp[fit[j].ipriors[i]] - pbar[i])/psigma[i])**2
        
    bestchisq  = currchisq
    
    #GENERATE RANDOM NUMBERS FOR MCMC
    numfree     = len(ifree)
    numnotfixed = len(inotfixed)
    step = np.random.normal(0, stepsize[ifree], (numit,numfree))
    unif = np.random.uniform(0, 1, numit)
    
    #START TIMER
    clock = timer.Timer(numit,progress = np.arange(0.05,1.01,0.05))
    
    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(0, numit):
        #Take step in random direction for adjustable parameters
        nextp[ifree] = params[ifree] + step[j]
        #CHECK FOR NEW STEPS OUTSIDE BOUNDARIES
        ioutside     = np.where(np.bitwise_or(nextp < pmin, nextp > pmax))[0]
        if (len(ioutside) > 0):
            nextp[ioutside]    = np.copy(params[ioutside])
            outside[ioutside] += 1
        #UPDATE PARAMTER(S) EQUAL TO OTHER PARAMETER(S)
        if (iequal.size > 0):
            for i in range(iequal.size):
                nextp[iequal[0][i]] = nextp[int(abs(stepsize[iequal[0][i]])-1)]
        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        pedit        = np.copy(nextp)
        nextchisq    = 0
        for k in range(numevents):
            ymodel     = np.ones(fit[k].nobj)
            m          = 0
            for i in range(cummodels[k],cummodels[k+1]):
                if   functype[i] == 'ortho':
                    #MODIFY COPY OF nextp ONLY
                    pedit[iortholist[k]] = myfuncs[i](pedit[iortholist[k]], funcx[i], fit[k].etc[m])
                elif functype[i] == 'ipmap':
                    ipmodel, allknots = myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], \
                                                      ymodel, retbinflux=True)
                    ymodel *= ipmodel
                    #COMPRESS allknots BY SUBTRACTING 1, MULTIPLYING BY 300000 THEN CONVERTING TO INT16
                    np.save(fit[k].allknotpid, np.int16((allknots[np.where(allknots > 0)]-1.)*300000))
                elif functype[i] == 'spline':
                    ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], ymodel)
                elif functype[i] == 'posoffset':
                    # Record change in Position 0 => cannot orthogonalize position parameters
                    ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                elif hasattr(fit[k], 'timebins') and (functype[i] == 'ecl/tr' 
                                                  or  functype[i] == 'ramp' 
                                                  or  functype[i] == 'sinusoidal'):
                    # Average over high-resolution model
                    hiresmodel = myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                    if len(fit[k].timebins) == fit[k].nobj:
                        for n in range(len(fit[k].timebins)):
                            ymodel[n] *= np.mean(hiresmodel[fit[k].timebins[n]])
                    else:
                        for n in range(len(fit[k].timebinsuc)):
                            ymodel[n] *= np.mean(hiresmodel[fit[k].timebinsuc[n]])
                elif functype[i] == 'noise':
                    noisepars = pedit[numparams[i]:numparams[i+1]]
                else:
                    ymodel *= myfuncs[i](pedit[numparams[i]:numparams[i+1]], funcx[i], fit[k].etc[m])
                m      += 1
            # Calculate chi^2
            if isrednoise == False:
                nextchisq  += mc.chisq(ymodel, y[k], sigma[k])
            else:
                nextchisq  += noisefunc(noisepars, ymodel-y[k], wavelet)
            # Apply prior, if one exists
            if len(fit[k].ipriors) > 0:
                pbar   = fit[k].priorvals[:,0]  #prior mean
                psigma = np.zeros(len(pbar))    #prior standard deviation
                # Determine psigma based on which side of asymmetric Gaussian nextp is on
                for i in range(len(fit[k].ipriors)):
                    if nextp[fit[j].ipriors[i]] < pbar[i]:
                        psigma[i] = fit[k].priorvals[i,1]
                    else:
                        psigma[i] = fit[k].priorvals[i,2]
                    nextchisq += ((nextp[fit[k].ipriors[i]] - pbar[i])/psigma[i])**2
                del(pbar,psigma)
            del(ymodel, m)
        
        #CALCULATE ACCEPTANCE PROBABILITY
        accept = np.exp(0.5 * (currchisq - nextchisq))
        if (accept >= 1) or (unif[j] <= accept):
            numaccept += 1
            params  = np.copy(nextp)
            currchisq  = nextchisq
            if (currchisq < bestchisq):
                bestp     = np.copy(params)
                bestchisq = currchisq
      	
        allparams[:, j] = params
        #PRINT INTERMEDIATE INFO
        if ((j+1) % intsteps == 0) and (j > 0): 
            print("\n" + time.ctime())
            print("Number of times parameter tries to step outside its prior:")
            print(outside)
            print("Current Best Parameters: ")
            print(bestp)
            #Flush allknots to file
            for k in range(numevents):
                for i in range(cummodels[k],cummodels[k+1]):
                    if functype[i] == 'ipmap':
                        fit[k].allknotpid.flush()
            
            #Apply Gelman-Rubin statistic
            if nchains > 0:
                psrf, meanpsrf = gr.convergetest(allparams[inotfixed,0:j+1],nchains)
                numconv = np.sum(psrf < 1.01)
                print("Gelman-Rubin statistic for free parameters:")
                print(psrf)
                if numconv == numnotfixed and j >= 1e4:
                    print("All parameters have converged to within 1% of unity. Halting MCMC.")
                    allparams = allparams[:,0:j+1]
                    break
            
        
        clock.check(j+1)
    
    return allparams, bestp, numaccept, j+1
