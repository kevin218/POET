#! usr/bin/env python
import orbit
import numpy as np
import models
import matplotlib.pyplot as plt
import mcmc

nparams = 7

def loadtr(obj):
	f = open(obj, "r")
	f.readline()
	raw = f.read().split('\n')
	data = np.zeros((len(raw)-1, 3))
	for i in range(0, len(raw)-1):
		data[i] = raw[i].split('\t')
	return np.transpose(data)

def loadparams(obj):
	f = open(obj, "r")
	f.readline()
	
	raw = f.read().split('\n')
	params = np.zeros((nparams, 4))
	for i in range(0, nparams):
		for j in range(0, 4):
			params[i][j] = raw[i].split('\t')[j+1]
	return np.transpose(params)

def execute(saveloc, numit, mag=False):
	rundir = "/home/esp01/events/dynamics/transits/"
	time, flux, error = loadtr(rundir+saveloc+"/data.txt")
	if mag == True:
		flux = 10**(-flux/2.5)
		error = 100**(-error/5.)-1
		print("Converting to relative flux...", error.shape)
	trdata = time, flux, error

	init, pmin, pmax, step = loadparams(rundir+saveloc+"/params.txt")		
	numparams = []
	numparams.append(0)
	numparams.append(len(init))
	allparams, best, numaccept = mcmc.mcmc(flux, time, init, pmin, pmax, step, numit, error, numparams, 1, [quadlimb], [time], 0)

	save(saveloc+"/output.npz", allparams, best, numaccept, trdata)
	print np.float(numaccept)/numit*100
	return allparams, best, numaccept, trdata

def save(output, allparams, best, numaccept, trdata):
	rundir = "/home/esp01/events/dynamics/transits/"
	np.savez(rundir+output, allparams=allparams, best=best, numaccept=numaccept, trdata=trdata)

def load(output):
	rundir = "/home/esp01/events/dynamics/transits/"
	foo = np.load(rundir+output)
	allparams = foo['allparams']
	best = foo['best']
	numaccept = foo['numaccept']
	trdata = foo['trdata']
	return allparams, best, numaccept, trdata

def quadlimb(params, x):
	midpt, width, rp_rs, b, flux, c1, c2 = params
	ingress = orbit.limbtime(b, width, 1, rp_rs)[0]
	trpars = np.array([midpt, width, rp_rs**2, ingress, ingress, flux])
	flux =  models.mandelecl(trpars, x)
	mu =  np.sqrt(1- b**2+((x-midpt)/width)**2)[np.where(np.abs(x-midpt) < width/2.)[0]]
	flux[np.where(np.abs(x-midpt) < width/2.)[0]] *= 1-c1*(1-mu)-c2*(1-mu)**2
	return flux

def analyze(output, saveloc=''):
	rundir = "/home/esp01/events/dynamics/transits/"
	allparams, best, numaccept = output[:3]
	time, flux, error = output[3]
	midpt, width, rp_rs, b, flux, c1, c2 = best
	names = [r'$T_0$', r'$D$', r'$R_p/R_s$', '$b$', r'$F$', r'$c_1$', r'$c_2$']
	print("\tBest\t\tMedian\t\tError")

	freeparams = list()
	numfp = 0
	for i in range(0, len(best)):
		print names[i], "	%f&	%f&	%f&" % (best[i], np.median(allparams[i]), np.std(allparams[i]))
		if np.median(allparams[i]) != best[i]:			
			freeparams.append([allparams[i], best[i], names[i]])
			numfp +=1

	model = quadlimb(best, output[3][0])

	chisq  = sum((model-flux)**2 / error**2)
	bic = chisq + numfp*np.log(len(model))
	print
	print "Free Parameters:", numfp
	print "MCMC Acceptance Rate:", 100*numaccept/float(allparams.shape[1]),"%"
	print "Chi-squared: ", chisq
	print "Reduced chi-squared: ", chisq/(len(model)-numfp)
	print "BIC: ", bic
	#print "RMS Residual:", rms

	#Plot data and model
	plt.figure(1)
	plt.errorbar(output[3][0], output[3][1], output[3][2], fmt='.')
	t = np.linspace(plt.xlim()[0], plt.xlim()[1], 1000)
	plt.plot(output[3][0], quadlimb(output[1], output[3][0]))
	plt.xlabel("Time (JD)")
	plt.ylabel("Relative Flux")
	if saveloc != '':
		plt.savefig(rundir+saveloc+'/data_model')

	#Histograms
	fp = 0
	dim = np.ceil(numfp/3.)
	j = 1
	skip = np.round(np.sqrt(len(allparams[0])))
	plt.figure(2, (8, 8))
	plt.clf()
	for i in range(0, numfp):
		plt.subplot(dim,3,i+1)
		plt.hist(freeparams[i][0][::skip], 20, label=str(freeparams[i][2]))
		plt.xticks(size=10,rotation=20)
		plt.yticks(size=7, rotation=0)
		#plt.axvline(best[i], label="Best-Fit Value")
		plt.xlabel(freeparams[i][2])
		#plt.legend()
	plt.show()
	#plt.subplots_adjust(left=0.08,right=0.92,bottom=0.07,top=0.95,hspace=0.5)
	if saveloc != '':
		plt.savefig(rundir+saveloc+'/histograms')
		plt.savefig(rundir+saveloc+'/histograms', format="pdf")
	plt.figure(3, (8, 8))
	plt.clf()
	j = 1
	for i in range(0, numfp):
		plt.subplot(numfp, 1, j)
		plt.plot(freeparams[i][0][::skip])
		plt.xlim((0, len(freeparams[1][0])/skip))
		plt.yticks(size=8)
		plt.xticks(size=8)
		plt.ylabel(freeparams[i][2])
		j+=1
		#plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95,hspace=0.5,wspace=0.2)
	plt.show()
	if saveloc != '':
		plt.savefig(rundir+saveloc+'/parameters')
		plt.savefig(rundir+saveloc+'/parameters', format="pdf")
	#plt.show()
	#Correlation plot
	k = 1
	plt.figure(4)
	plt.clf()
	for i in range(1, numfp):
		for j in range(0, numfp-1):
			if i > j:
				plt.subplot(numfp-1, numfp-1, k)
				plt.plot(freeparams[j][0][::skip], freeparams[i][0][::skip], ',')
				plt.xticks(size=0,rotation=90)
				plt.yticks(size=0)
				if j == 0:
					#plt.yticks(size=8)
					s = freeparams[i][2]
					plt.ylabel(s, size = 8)
				if i == numfp-1:
					#plt.xticks(size=8)
					s = freeparams[j][2]
					plt.xlabel(s, size = 8)
			k+=1
	#plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95,hspace=0.2,wspace=0.2)
	plt.show()
	if saveloc != '':
		plt.savefig(rundir+saveloc+'/correlations')
		plt.savefig(rundir+saveloc+'/correlations', format="pdf")
	
	
