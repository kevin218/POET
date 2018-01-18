
# Centroid analysis
# This function performs image differencing between in- and out-of-transit frames
# Result is a PSF at the true spatial location of the transit source
# Amplitude equals the transit depth times the image intensity for target
# Centroid of direct image - centroid of difference image should be ~0
# 3*formal error (3 sigma) on PSF fit of differenced image = arcsec limit on background eclipsing binaries.

import sys, os
r = os.getcwd().split("/")
maindir = "/".join(r[:r.index("run")])
sys.path.append(maindir + '/lib/')

import manageevent  as me
import centerdriver as cd
tempevent = me.loadevent('../../gj436cp22_ctr', load=['data','uncd','mask'])
event[0].data = tempevent.data
del tempevent

'''
hist2d, xedges, yedges = np.histogram2d(event.fp.x[0,ipretr[0]:ipretr[1]], event.fp.y[0,ipretr[0]:ipretr[1]],20)
plt.figure(2)
plt.clf()
plt.suptitle('Pre-Transit')
a=plt.subplot(111)
a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.2f'))
plt.imshow(hist2d.T,extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]), cmap=cm.gray_r, 
                            aspect='auto', origin='lower')
plt.colorbar()

hist2d, xedges, yedges = np.histogram2d(event.fp.x[0,iintr[0]:iintr[1]], event.fp.y[0,iintr[0]:iintr[1]],20)
plt.figure(3)
plt.clf()
plt.suptitle('Transit')
a=plt.subplot(111)
a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.2f'))
plt.imshow(hist2d.T,extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]), cmap=cm.gray_r, 
                            aspect='auto', origin='lower')
plt.colorbar()

hist2d, xedges, yedges = np.histogram2d(event.fp.x[0,iposttr[0]:iposttr[1]], event.fp.y[0,iposttr[0]:iposttr[1]],20)
plt.figure(4)
plt.clf()
plt.suptitle('Post-Transit')
a=plt.subplot(111)
a.yaxis.set_major_formatter(plt.matplotlib.ticker.FormatStrFormatter('%0.2f'))
plt.imshow(hist2d.T,extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]), cmap=cm.gray_r, 
                            aspect='auto', origin='lower')
plt.colorbar()

'''

#****NEED TO DEAL WITH POSITION-DEPENDENT FLUX****

def centroidAnal():

pos     = 0

pretr   = [2455772.775, 2455772.785]
intr    = [2455772.795, 2455772.815]
posttr  = [2455772.825, 2455772.835]

# Range of indices for pre-transit, transit, and post-transit
ipretr  = np.where(np.bitwise_and(event[0].bjdutc[pos] >= pretr[0],  event[0].bjdutc[pos] <= pretr[1]))[0]
ipretr  = [ipretr[0],ipretr[-1]]
iintr   = np.where(np.bitwise_and(event[0].bjdutc[pos] >= intr[0],   event[0].bjdutc[pos] <= intr[1]))[0]
iintr   = [iintr[0],iintr[-1]]
iposttr = np.where(np.bitwise_and(event[0].bjdutc[pos] >= posttr[0], event[0].bjdutc[pos] <= posttr[1]))[0]
iposttr = [iposttr[0],iposttr[-1]]

#ipretr  = [6000,12000]
#iintr   = [13000,19000]
#iposttr = [20000,26000]
binsize = [0.005,0.008] #[y,x]
#binsize = [0.005,0.010] #[y,x]
srcest  = [event[0].y[pos,ipretr[0]:iposttr[1]].mean(),event[0].x[pos,ipretr[0]:iposttr[1]].mean()]
#[14.9769, 14.7046]

# Apply BLISS map correction to frames
bestmip = np.ones(event[0].nimpos[pos])
#bestmip[np.where(event[0].good[pos])] = event[0].fit[0].bestmipuc
"""
predata  = (event[0].data[ipretr[0]:ipretr[1],:,:,pos].T/bestmip[ipretr[0]:ipretr[1]]).T
predata  = predata[np.where(event[0].good[pos,ipretr[0]:ipretr[1]])[0]]

postdata = (event[0].data[iposttr[0]:iposttr[1],:,:,pos].T/bestmip[iposttr[0]:iposttr[1]]).T
postdata = postdata[np.where(event[0].good[pos,iposttr[0]:iposttr[1]])[0]]

trdata = (event[0].data[iintr[0]:iintr[1],:,:,pos].T/bestmip[iintr[0]:iintr[1]]).T
trdata = trdata[np.where(event[0].good[pos,iintr[0]:iintr[1]])[0]]
"""

# Range of pixel positions to consider in above ranges
# Used to minimized position-dependent flux variations
yrange  = [srcest[0]-binsize[0]/2.,srcest[0]+binsize[0]/2.]
xrange  = [srcest[1]-binsize[1]/2.,srcest[1]+binsize[1]/2.]

iyrange = np.where(np.bitwise_and(event[0].y[pos,ipretr[0]:ipretr[1]] >= yrange[0],event[0].y[pos,ipretr[0]:ipretr[1]] <= yrange[1]))[0]
irange  = iyrange[np.where(np.bitwise_and(event[0].x[pos,ipretr[0]:ipretr[1]][iyrange] >= xrange[0],event[0].x[pos,ipretr[0]:ipretr[1]][iyrange] <= xrange[1]))]
predata = (event[0].data[ipretr[0]:ipretr[1]][irange,:,:,pos].T/bestmip[ipretr[0]:ipretr[1]][irange]).T
'''
print(np.round(srcest[0]-event[0].y[0,ipretr[0]:ipretr[1]][irange].mean(),4), event[0].y[0,ipretr[0]:ipretr[1]][irange].std())
print(np.round(srcest[1]-event[0].x[0,ipretr[0]:ipretr[1]][irange].mean(),4), event[0].x[0,ipretr[0]:ipretr[1]][irange].std())

'''
# Post-transit
iyrange  = np.where(np.bitwise_and(event[0].y[pos,iposttr[0]:iposttr[1]] >= yrange[0],event[0].y[pos,iposttr[0]:iposttr[1]] <= yrange[1]))[0]
irange   = iyrange[np.where(np.bitwise_and(event[0].x[pos,iposttr[0]:iposttr[1]][iyrange] >= xrange[0],event[0].x[pos,iposttr[0]:iposttr[1]][iyrange] <= xrange[1]))]
postdata = (event[0].data[iposttr[0]:iposttr[1]][irange,:,:,pos].T/bestmip[iposttr[0]:iposttr[1]][irange]).T
'''
print(np.round(srcest[0]-event[0].y[0,iposttr[0]:iposttr[1]][irange].mean(),4), event[0].y[0,iposttr[0]:iposttr[1]][irange].std())
print(np.round(srcest[1]-event[0].x[0,iposttr[0]:iposttr[1]][irange].mean(),4), event[0].x[0,iposttr[0]:iposttr[1]][irange].std())

'''
# Transit
iyrange = np.where(np.bitwise_and(event[0].y[pos,iintr[0]:iintr[1]] >= yrange[0],event[0].y[pos,iintr[0]:iintr[1]] <= yrange[1]))[0]
irange  = iyrange[np.where(np.bitwise_and(event[0].x[pos,iintr[0]:iintr[1]][iyrange] >= xrange[0],event[0].x[pos,iintr[0]:iintr[1]][iyrange] <= xrange[1]))]
trdata = (event[0].data[iintr[0]:iintr[1]][irange,:,:,pos].T/bestmip[iintr[0]:iintr[1]][irange]).T
'''
print(np.round(srcest[0]-event[0].y[0,iintr[0]:iintr[1]][irange].mean(),4), event[0].y[0,iintr[0]:iintr[1]][irange].std())
print(np.round(srcest[1]-event[0].x[0,iintr[0]:iintr[1]][irange].mean(),4), event[0].x[0,iintr[0]:iintr[1]][irange].std())

'''

"""
# Range of pixel positions to consider in above ranges
# Used to minimized position-dependent flux variations
yrange  = [srcest[0]-binsize[0]/2.,srcest[0]+binsize[0]/2.]
xrange  = [srcest[1]-binsize[1]/2.,srcest[1]+binsize[1]/2.]

# Calculate indices and data frames for above range constraints
# Pre-transit
iyrange = np.where(np.bitwise_and(event.fp.y[pos,ipretr[0]:ipretr[1]] >= yrange[0],event.fp.y[pos,ipretr[0]:ipretr[1]] <= yrange[1]))[0]
irange  = iyrange[np.where(np.bitwise_and(event.fp.x[pos,iyrange] >= xrange[0],event.fp.x[pos,[iyrange]] <= xrange[1]))[1]]
predata = event.data[irange,:,:,pos]
'''
print(event.fp.y[0,irange].mean(), event.fp.y[0,irange].std())
(14.989822580073428, 0.011886811778217814)
print(event.fp.x[0,irange].mean(), event.fp.x[0,irange].std())
(14.709697130397215, 0.0021216186389019735)
'''
# Post-transit
iyrange = np.where(np.bitwise_and(event.fp.y[pos,iposttr[0]:iposttr[1]] >= yrange[0],event.fp.y[pos,iposttr[0]:iposttr[1]] <= yrange[1]))[0]
irange  = iyrange[np.where(np.bitwise_and(event.fp.x[pos,iyrange] >= xrange[0],event.fp.x[pos,[iyrange]] <= xrange[1]))[1]]
postdata = event.data[irange,:,:,pos]
'''
print(event.fp.y[0,irange].mean(), event.fp.y[0,irange].std())
(14.984410812461936, 0.010755583543844819)
print(event.fp.x[0,irange].mean(), event.fp.x[0,irange].std())
(14.709008995749002, 0.0016089546347219279)
'''
# Transit
iyrange = np.where(np.bitwise_and(event.fp.y[pos,iintr[0]:iintr[1]] >= yrange[0],event.fp.y[pos,iintr[0]:iintr[1]] <= yrange[1]))[0]
irange  = iyrange[np.where(np.bitwise_and(event.fp.x[pos,iyrange] >= xrange[0],event.fp.x[pos,[iyrange]] <= xrange[1]))[1]]
trdata = event.data[irange,:,:,pos]
'''
print(event.fp.y[0,irange].mean(), event.fp.y[0,irange].std())
(14.987706659772504, 0.015604372475966339)
print(event.fp.x[0,irange].mean(), event.fp.x[0,irange].std())
(14.70951534100972, 0.0024232686852807416)
'''
"""

# Compute difference (frame subtraction)
trout   = np.mean(np.concatenate((predata,postdata)),axis=0)
trin    = np.mean(trdata,axis=0)
diff    = trout - trin

# Plot
plt.figure(pos+1, figsize=(8,3))
plt.clf()
plt.subplot(121)
plt.title('Direct Image')
plt.imshow(trout/1000., cmap=plt.cm.gray_r, origin='lower')
plt.colorbar()
plt.subplot(122)
plt.title('Difference Image - GJ 436c')
plt.imshow(diff/1000., cmap=plt.cm.gray_r, origin='lower')
plt.colorbar()

# Perform centroiding
posout  = cd.centerdriver('fgc', trout, srcest, 4, None, None,
                          mask=None, uncd=None, fitbg=1, maskstar=True,
                          expand=5.0, psf=None, psfctr=None)[0]
posin   = cd.centerdriver('fgc', trin, srcest, 4, None, None,
                          mask=None, uncd=None, fitbg=1, maskstar=True,
                          expand=5.0, psf=None, psfctr=None)[0]
# Use Least Asymmetry due to weak signal
posdiff = cd.centerdriver('lag', diff, srcest, 0, radius=7.0, size=6.0,
                          mask=None, uncd=None, fitbg=1, maskstar=True,
                          expand=5.0, psf=None, psfctr=None)[0]

'''
posdiff = cd.centerdriver('fgc', diff, srcest, 4, radius=7.0, size=4.0,
                          mask=None, uncd=None, fitbg=1, maskstar=True,
                          expand=5.0, psf=None, psfctr=None)[0]
'''
'''
print(posout-posdiff)
[-0.0201039  -0.28450996]
# Find magnitude
np.sqrt(np.sum((posout-posdiff)**2))
0.28521936178381485
# Assuming 1.2 "/pix for Spitzer, centroid difference = 0.34226 "/pix
#**** NEED UNCERTAINTY CALC***
# Amplitude equals the transit depth times the image intensity for target
'''
import apphot as ap

fluxout  = ap.apphot(trout,  posout, photap=5.25, skyin=10., skyout=30., expand=5)
fluxin   = ap.apphot( trin,   posin, photap=5.25, skyin=10., skyout=30., expand=5)
fluxdiff = ap.apphot( diff, posdiff, photap=5.25, skyin=10., skyout=30., expand=5)
#fluxdiff = ap.apphot( diff, srcest, photap=5.25, skyin=10., skyout=30., expand=5)

print(fluxdiff/fluxout*1e6)


