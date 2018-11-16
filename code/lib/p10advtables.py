#! /usr/bin/env python

# $Author: kevin $
# $Revision: 549 $
# $Date: 2011-08-24 12:40:35 -0400 (Wed, 24 Aug 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/pipeline/trunk/p10advtables.py $
# $Id: p10advtables.py 549 2011-08-24 16:40:35Z kevin $

'''
To run this program, please fill in all the fields commented
'ENTER MANUALLY'.  Save and close the file then type 'p10advtables'
to run the program.  This version is specially designed for BLISS 
mapping and exponential & polynomial ramp models.  If you do not 
have a fit that uses these, please contact ccampo/kevin to get a 
version that works for your model fits.
'''

import os, sys
import cPickle as pickle
import numpy             as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
sys.path.append('../')
import run
plt.ion()

# import param files
modeldir = "/home/esp01/ancil/modelparams/"
obj      = ['hd149bs11','hd149bs21','hd149bs31','hd149bs32','hd149bs33', \
            'hd149bs41','hd149bs42','hd149bs43','hd149bs44','hd149bs51','hd149bs52']
sys.path.append(modeldir)
for i in range(len(obj)):
    exec 'import ' + obj[i] + 'params'

#RESTORE ALL SAVE FILES
event = run.p7Restore(filedir='.',topdir='.')
#event = run.p6Restore(filedir='.',topdir='.')

'''
#MANUALLY RESTORE SAVE FILES
sys.path.append('lib/')
loadfile = []
for fname in os.listdir('.'):
    if (fname.endswith("6model.dat")):
        loadfile.append(fname)

loadfile.sort()
numevents = len(loadfile)
event    = []
for i in range(5,len(loadfile)):
    print("Loading " + loadfile[i])
    handle      = open(loadfile[i], 'r')
    event.append(pickle.load(handle))
    handle.close()

'''

ephoffset = 2450000     #MJD - BJD in table, ENTER MANUALLY

# individual fits - ENTER MANUALLY
hd149bs11, hd149bs21, hd149bs31, hd149bs32, hd149bs33, \
hd149bs41, hd149bs42, hd149bs43, hd149bs44, hd149bs51, hd149bs52 = event
#fit = [event[0].fit[0], event[1].fit[0], event[2].fit[0], event[3].fit[0]]
fit    = [hd149bs11.fit[0], 
          hd149bs21.fit[0], 
          hd149bs31.fit[0],
          hd149bs32.fit[0],
          hd149bs33.fit[0],
          hd149bs41.fit[0],
          hd149bs42.fit[0],
          hd149bs43.fit[0],
          hd149bs44.fit[0],
          hd149bs51.fit[0],
          hd149bs52.fit[0]]

# is this a joint fit? If so, specify which models were run jointly. ENTER MANUALLY
joint     = True
jointruns = [[0],
             [1],
             [2, 3, 4],
             [5, 6, 7, 8],
             [9, 10]]

# wavelengths/labels of each channel - ENTER MANUALLY
#wavelen = ["3.6 um", "4.5 um", "5.8 um", "8.0 um"]
wavelen = ['HD149bs11', 'HD149bs21', 'HD149bs31', 'HD149bs32', 'HD149bs33',
           'HD149bs41', 'HD149bs42', 'HD149bs43', 'HD149bs44', 'HD149bs51', 'HD149bs52']

# filename to save as - ENTER MANUALLY
# set to `None` if you do not wish to save a file
fname = "hd149026b_lightcurve_fit_table-2011-04-09.txt"

    

#### DO NOT EDIT THE SECTION BELOW UNLESS YOU KNOW WHAT YOU'RE DOING ####
#### CONTINUE TO THE END OF THE CODE TO SELECT TABLE SIZE.
'''
# generate position consistencies
def pconsist(fit):
    pcosx = []
    pcosy = []
    for i in range(len(fit)):
        xbar = fit[i].x[:-1] - fit[i].x[1:]
        ybar = fit[i].y[:-1] - fit[i].y[1:]
        px   = np.sqrt(np.sum(xbar**2)/xbar.size)
        py   = np.sqrt(np.sum(ybar**2)/ybar.size)
        pcosx.append(px)
        pcosy.append(py)
    return np.array(pcosx), np.array(pcosy)
'''
#CALCULATE ecltimeutc and ecltimetdb
def ephcalc(event, fit, midpt):
    if hasattr(event, 'timestd'):
        print('Verify that the ephemeris is reported in BJD_' + event.timestd + '!')
        offset = event.bjdtdb.flat[0]-event.bjdutc.flat[0]
        if   event.timestd == 'utc':
            ephtimeutc = event.ephtime
            ephtimetdb = event.ephtime + offset
        elif event.timestd == 'tdb':
            ephtimetdb = event.ephtime
            ephtimeutc = event.ephtime - offset
        else:
            print('Assuming that the ephemeris is reported in BJD_UTC!')
            ephtimeutc = event.ephtime
            ephtimetdb = event.ephtime + offset
    else:
        print('Assuming that the ephemeris is reported in BJD_UTC!')
        offset = event.bjdtdb.flat[0]-event.bjdutc.flat[0]
        ephtimeutc = event.ephtime
        ephtimetdb = event.ephtime + offset
    fit.ecltimeutc = (np.floor((event.bjdutc.flat[0] - ephtimeutc)/event.period) + \
                      fit.bestp[midpt]) * event.period + ephtimeutc
    fit.ecltimetdb = (np.floor((event.bjdtdb.flat[0] - ephtimetdb)/event.period) + \
                      fit.bestp[midpt]) * event.period + ephtimetdb
    return

# MAKE ARRAYS OF EACH PARAMETER FOR EACH CHANNEL
nfit    = len(fit)
xpos    = []
ypos    = []
pcosx   = []
pcosy   = []
#pcosx   = np.round(pconsist(fit)[0], decimals=3).tolist()
#pcosy   = np.round(pconsist(fit)[1], decimals=3).tolist()
photap  = []
skyin   = []
skyout  = []
flux    = []
depth   = []
temp    = []
midpt   = []
rprs    = []
cosi    = []
ars     = []
midtimeutc = []
midtimetdb = []
width   = []
t12     = []
t34     = []
ramp    = []
r0      = []
r1      = []
r2      = []
r3      = []
r4      = []
r5      = []
r6      = []
r7      = []
t0      = []
isbliss = []
minnumpts = []
totfrm  = []
goodfrm = []
rejfrm  = []
freepar = []
N       = []
aic     = []
bic     = []
sdnr    = []
chifact = []
photonSNR = []
# get each parameter with the proper ammount of decimal places and error
for j in range(nfit):
    xpos    .append( np.round(fit[j].xuc.mean(), decimals=2) )
    ypos    .append( np.round(fit[j].yuc.mean(), decimals=2) )
    pcosx   .append( np.round(event[j].xprecision, decimals=3) )
    pcosy   .append( np.round(event[j].yprecision, decimals=3) )
    photap  .append( event[j].photap )
    skyin   .append( event[j].skyin )
    skyout  .append( event[j].skyout )
    if hasattr(fit[j],'tbm'):
        temp    .append( str(np.round(fit[j].tbm)) + " & " + str(np.round(fit[j].tbsd)) )
    if hasattr(fit[j].i,'midpt'):
        midpt   .append( str(np.round(fit[j].bestp[fit[j].i.midpt], decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.midpt], decimals=4)) )
        depth   .append( str(np.round(fit[j].bestp[fit[j].i.depth]*100, decimals=3)) + " & " + str(np.round(fit[j].typerr[fit[j].i.depth]*100, decimals=3)) )
        width   .append( str(np.round(event[j].period*24*fit[j].bestp[fit[j].i.width], decimals=2)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.width], decimals=2)) )
        t12     .append( str(np.round(24*event[j].period*fit[j].bestp[fit[j].i.t12], decimals=3)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.t12], decimals=3)) )
        t34     .append( str(np.round(24*event[j].period*fit[j].bestp[fit[j].i.t34], decimals=3)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.t34], decimals=3)) )
        flux    .append( str(np.round(fit[j].bestp[fit[j].i.flux])) + " & " + str(np.round(fit[j].typerr[fit[j].i.flux])) )
    elif hasattr(fit[j].i,'midpt2'):
        midpt   .append( str(np.round(fit[j].bestp[fit[j].i.midpt2], decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.midpt2], decimals=4)) )
        depth   .append( str(np.round(fit[j].bestp[fit[j].i.depth2]*100, decimals=3)) + " & " + str(np.round(fit[j].typerr[fit[j].i.depth2]*100, decimals=3)) )
        width   .append( str(np.round(event[j].period*24*fit[j].bestp[fit[j].i.width2], decimals=2)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.width2], decimals=2)) )
        t12     .append( str(np.round(24*event[j].period*fit[j].bestp[fit[j].i.t122], decimals=3)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.t122], decimals=3)) )
        t34     .append( str(np.round(24*event[j].period*fit[j].bestp[fit[j].i.t342], decimals=3)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.t342], decimals=3)) )
        flux    .append( str(np.round(fit[j].bestp[fit[j].i.flux2])) + " & " + str(np.round(fit[j].typerr[fit[j].i.flux2])) )
    elif hasattr(fit[j].i,'midpt3'):
        midpt   .append( str(np.round(fit[j].bestp[fit[j].i.midpt3], decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.midpt3], decimals=4)) )
        depth   .append( str(np.round(fit[j].bestp[fit[j].i.depth3]*100, decimals=3)) + " & " + str(np.round(fit[j].typerr[fit[j].i.depth3]*100, decimals=3)) )
        width   .append( str(np.round(event[j].period*24*fit[j].bestp[fit[j].i.width3], decimals=2)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.width3], decimals=2)) )
        t12     .append( str(np.round(24*event[j].period*fit[j].bestp[fit[j].i.t123], decimals=3)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.t123], decimals=3)) )
        t34     .append( str(np.round(24*event[j].period*fit[j].bestp[fit[j].i.t343], decimals=3)) + " & "
                         + str(np.round(24*event[j].period*fit[j].typerr[fit[j].i.t343], decimals=3)) )
        flux    .append( str(np.round(fit[j].bestp[fit[j].i.flux3])) + " & " + str(np.round(fit[j].typerr[fit[j].i.flux3])) )
    elif hasattr(fit[j].i,'trspmid'):
        mjdoffset = event[j].params.tuoffset - ephoffset
        midtimeutc.append( str(np.round(fit[j].bestp[fit[j].i.trspmid]+mjdoffset, decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.trspmid], decimals=4)) )
        midtimetdb.append( str(np.round(fit[j].bestp[fit[j].i.trspmid]+mjdoffset+65.183/86400, decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.trspmid], decimals=4)) )
        rprs    .append( str(np.round(fit[j].bestp[fit[j].i.rprs], decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.rprs], decimals=4)) )
        cosi    .append( str(np.round(fit[j].bestp[fit[j].i.cosi], decimals=6)) + " & "
                         + str(np.round(fit[j].typerr[fit[j].i.cosi], decimals=6)) )
        ars     .append( str(np.round(fit[j].bestp[fit[j].i.ars], decimals=3)) + " & "
                         + str(np.round(fit[j].typerr[fit[j].i.ars], decimals=3)) )
        flux    .append( str(np.round(fit[j].bestp[fit[j].i.trspf])) + " & " + str(np.round(fit[j].typerr[fit[j].i.trspf])) )
    elif hasattr(fit[j].i,'trspmid2'):
        mjdoffset = event[j].params.tuoffset - ephoffset
        midtimeutc.append( str(np.round(fit[j].bestp[fit[j].i.trspmid2]+mjdoffset, decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.trspmid2], decimals=4)) )
        midtimetdb.append( str(np.round(fit[j].bestp[fit[j].i.trspmid2]+mjdoffset+65.183/86400, decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.trspmid2], decimals=4)) )
        rprs    .append( str(np.round(fit[j].bestp[fit[j].i.rprs2], decimals=4)) + " & " + str(np.round(fit[j].typerr[fit[j].i.rprs2], decimals=4)) )
        cosi    .append( str(np.round(fit[j].bestp[fit[j].i.cosi2], decimals=6)) + " & "
                         + str(np.round(fit[j].typerr[fit[j].i.cosi2], decimals=6)) )
        ars     .append( str(np.round(fit[j].bestp[fit[j].i.ars2], decimals=3)) + " & "
                         + str(np.round(fit[j].typerr[fit[j].i.ars2], decimals=3)) )
        flux    .append( str(np.round(fit[j].bestp[fit[j].i.trspf2])) + " & " + str(np.round(fit[j].typerr[fit[j].i.trspf2])) )
    if hasattr(fit[j],'ecltime'):
        ephcalc(event[j], fit[j], fit[j].i.midpt)
        del fit[j].ecltime
    elif hasattr(fit[j],'ecltime2'):
        ephcalc(event[j], fit[j], fit[j].i.midpt2)
        del fit[j].ecltime2
    elif hasattr(fit[j],'ecltime3'):
        ephcalc(event[j], fit[j], fit[j].i.midpt3)
        del fit[j].ecltime3
    if hasattr(fit[j],'ecltimeutc'):
        midtimeutc.append( str(np.round(fit[j].ecltimeutc, decimals=4)-ephoffset) + " & " + str(np.round(fit[j].ecltimeerr, decimals=4)) )
    elif hasattr(fit[j],'ecltimeutc2'):
        midtimeutc.append( str(np.round(fit[j].ecltimeutc2, decimals=4)-ephoffset) + " & " + str(np.round(fit[j].ecltimeerr2, decimals=4)) )
    elif hasattr(fit[j],'ecltimeutc3'):
        midtimeutc.append( str(np.round(fit[j].ecltimeutc3, decimals=4)-ephoffset) + " & " + str(np.round(fit[j].ecltimeerr3, decimals=4)) )
    if hasattr(fit[j],'ecltimetdb'):
        midtimetdb.append( str(np.round(fit[j].ecltimetdb, decimals=4)-ephoffset) + " & " + str(np.round(fit[j].ecltimeerr, decimals=4)) )
    elif hasattr(fit[j],'ecltimetdb2'):
        midtimetdb.append( str(np.round(fit[j].ecltimetdb2, decimals=4)-ephoffset) + " & " + str(np.round(fit[j].ecltimeerr2, decimals=4)) )
    elif hasattr(fit[j],'ecltimetdb3'):
        midtimetdb.append( str(np.round(fit[j].ecltimetdb3, decimals=4)-ephoffset) + " & " + str(np.round(fit[j].ecltimeerr3, decimals=4)) )
    ramp    .append( str(fit[j].model[2]) )
    #r0
    if hasattr(fit[j].i,'rem'):
        r0  .append( str(np.round(fit[j].bestp[fit[j].i.rem],2)) + " & " + str(np.round(fit[j].typerr[fit[j].i.rem],2)) )
    elif hasattr(fit[j].i,'fem'):
        r0  .append( str(np.round(fit[j].bestp[fit[j].i.fem],2)) + " & " + str(np.round(fit[j].typerr[fit[j].i.fem],2)) )
    elif hasattr(fit[j].i,'re2m1'):
        r0  .append( str(np.round(fit[j].bestp[fit[j].i.re2m1],2)) + " & " + str(np.round(fit[j].typerr[fit[j].i.re2m1],2)) )
    elif hasattr(fit[j].i,'ser0'):
        r0  .append( str(np.round(fit[j].bestp[fit[j].i.ser0],2)) + " & " + str(np.round(fit[j].typerr[fit[j].i.ser0],2)) )
    elif hasattr(fit[j].i,'selr0'):
        r0  .append( str(np.round(fit[j].bestp[fit[j].i.selr0],2)) + " & " + str(np.round(fit[j].typerr[fit[j].i.selr0],2)) )
    elif hasattr(fit[j].i,'seqr0'):
        r0  .append( str(np.round(fit[j].bestp[fit[j].i.seqr0],2)) + " & " + str(np.round(fit[j].typerr[fit[j].i.seqr0],2)) )
    elif hasattr(fit[j].i,'se2r0'):
        r0  .append( str(np.round(fit[j].bestp[fit[j].i.se2r0],2)) + " & " + str(np.round(fit[j].typerr[fit[j].i.se2r0],2)) )
    else:
        r0  .append( int(0) )
    #r1
    if hasattr(fit[j].i,'ser1'):
        r1  .append( str(np.round(fit[j].bestp[fit[j].i.ser1],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.ser1],6)) )
    elif hasattr(fit[j].i,'selr1'):
        r1  .append( str(np.round(fit[j].bestp[fit[j].i.selr1],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.selr1],6)) )
    elif hasattr(fit[j].i,'seqr1'):
        r1  .append( str(np.round(fit[j].bestp[fit[j].i.seqr1],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.seqr1],6)) )
    elif hasattr(fit[j].i,'se2r1'):
        r1  .append( str(np.round(fit[j].bestp[fit[j].i.se2r1],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.se2r1],6)) )
    elif hasattr(fit[j].i,'ret0'):
        r1 .append( str(np.round(fit[j].bestp[fit[j].i.ret0], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.ret0], decimals=6)) )
    elif hasattr(fit[j].i,'fet0'):
        r1 .append( str(np.round(fit[j].bestp[fit[j].i.fet0], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.fet0], decimals=6)) )
    elif hasattr(fit[j].i,'re2t1'):
        r1 .append( str(np.round(fit[j].bestp[fit[j].i.re2t1], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.re2t1], decimals=6)) )
    else:
        r1 .append( int(0) )
    #r2
    if hasattr(fit[j].i,'linr2'):
        r2 .append( str(np.round(fit[j].bestp[fit[j].i.linr2], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.linr2], decimals=6)) )
    elif hasattr(fit[j].i,'seqr2'):
        r2  .append( str(np.round(fit[j].bestp[fit[j].i.seqr2],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.seqr2],6)) )
    elif hasattr(fit[j].i,'selr2'):
        r2  .append( str(np.round(fit[j].bestp[fit[j].i.selr2],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.selr2],6)) )
    elif hasattr(fit[j].i,'qrr2'):
        r2 .append( str(np.round(fit[j].bestp[fit[j].i.qrr2], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.qrr2], decimals=6)) )
    else:
	    r2 .append( int(0) )
    #r3
    if hasattr(fit[j].i,'qrr3'):
        r3 .append( str(np.round(fit[j].bestp[fit[j].i.qrr3], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.qrr3], decimals=6)) )
    elif hasattr(fit[j].i,'seqr3'):
        r3  .append( str(np.round(fit[j].bestp[fit[j].i.seqr3],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.seqr3],6)) )
    else:
	    r3 .append( int(0) )
    #r4
    if hasattr(fit[j].i,'re2m2'):
        r4  .append( str(np.round(fit[j].bestp[fit[j].i.re2m2],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.re2m2],6)) )
    elif hasattr(fit[j].i,'se2r4'):
        r4  .append( str(np.round(fit[j].bestp[fit[j].i.se2r4],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.se2r4],6)) )
    else:
        r4  .append( int(0) )
    #r5
    if hasattr(fit[j].i,'re2t2'):
        r5 .append( str(np.round(fit[j].bestp[fit[j].i.re2t2], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.re2t2], decimals=6)) )
    elif hasattr(fit[j].i,'se2r5'):
        r5  .append( str(np.round(fit[j].bestp[fit[j].i.se2r5],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.se2r5],6)) )
    else:
        r5  .append( int(0) )
    #r6
    if hasattr(fit[j].i,'llr6'):
        r6 .append( str(np.round(fit[j].bestp[fit[j].i.llr6], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.llr6], decimals=6)) )
    elif hasattr(fit[j].i,'lqr6'):
        r6  .append( str(np.round(fit[j].bestp[fit[j].i.lqr6],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.lqr6],6)) )
    elif hasattr(fit[j].i,'logr6'):
        r6  .append( str(np.round(fit[j].bestp[fit[j].i.logr6],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.logr6],6)) )
    elif hasattr(fit[j].i,'l4qr6'):
        r6  .append( str(np.round(fit[j].bestp[fit[j].i.l4qr6],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.l4qr6],6)) )
    else:
        r6  .append( int(0) )
    #r7
    if hasattr(fit[j].i,'lqr7'):
        r7  .append( str(np.round(fit[j].bestp[fit[j].i.lqr7],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.lqr7],6)) )
    elif hasattr(fit[j].i,'logr7'):
        r7  .append( str(np.round(fit[j].bestp[fit[j].i.logr7],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.logr7],6)) )
    elif hasattr(fit[j].i,'l4qr7'):
        r7  .append( str(np.round(fit[j].bestp[fit[j].i.l4qr7],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.l4qr7],6)) )
    else:
        r7  .append( int(0) )
    #t0
    if hasattr(fit[j].i,'llt0'):
        t0 .append( str(np.round(fit[j].bestp[fit[j].i.llt0], decimals=6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.llt0], decimals=6)) )
    elif hasattr(fit[j].i,'lqt0'):
        t0  .append( str(np.round(fit[j].bestp[fit[j].i.lqt0],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.lqt0],6)) )
    elif hasattr(fit[j].i,'logt0'):
        t0  .append( str(np.round(fit[j].bestp[fit[j].i.logt0],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.logt0],6)) )
    elif hasattr(fit[j].i,'l4qt0'):
        t0  .append( str(np.round(fit[j].bestp[fit[j].i.l4qt0],6)) + " & " + str(np.round(fit[j].typerr[fit[j].i.l4qt0],6)) )
    else:
        t0  .append( int(0) )
    #BLISS map
    if hasattr(fit[j].i,'blip'):
        isbliss .append("Yes")
        minnumpts.append(str(event[j].params.minnumpts[0]))
    else:
        isbliss .append("No")
        minnumpts.append("-")
    totfrm  .append( event[j].phase.size )
    goodfrm .append( fit[j].nobj )
    rejfrm  .append( 100.*(totfrm[j] - sum(event[j].good.flatten()))/totfrm[j] )
    freepar .append( fit[j].numfreepars )
    N       .append( fit[j].nobj )
    aic     .append( fit[j].aic )
    bic     .append( fit[j].bic )
    sdnr    .append( fit[j].normresiduals.std() )
    chifact .append( np.mean(fit[j].newsigma/fit[j].sigma) )
    photonSNR.append(str(np.round(100*fit[j].sn/fit[j].photsn, decimals=1)))

# CHECK FOR JOINT FITS AND UPDATE BIC AND NUMBER OF FREE PARAMETERS
if joint:
    for i in range(len(jointruns)):
        a    = 0
        b    = 0
        nobj = 0
        numfreepars = 0
        for j in jointruns[i]:
            a  += sum((fit[j].residuals/event[j].fit[0].newsigma)**2)
            b  += sum((fit[j].residuals/event[j].fit[0].newsigma)**2)
            nobj += fit[j].nobj
            numfreepars += fit[j].numfreepars
        for j in jointruns[i]:
            aic[j] = a + 2*numfreepars
            bic[j] = b + numfreepars*np.log(nobj)
        
# add zeros for missing channels
allparams = [
    wavelen ,
    xpos    ,
    ypos    ,
    pcosx   ,
    pcosy   ,
    photap  ,
    skyin   ,
    skyout  ,
    flux    ,
    depth   ,
    temp    ,
    midpt   ,
    rprs    ,
    cosi    ,
    ars     ,
    midtimeutc ,
    midtimetdb ,
    width   ,
    t12     ,
    t34     ,
    ramp    ,
    r0      ,
    r1      ,
    r2      ,
    r3      ,
    r4      ,
    r5      ,
    r6      ,
    r7      ,
    t0      ,
    isbliss ,
    minnumpts,
    totfrm  ,
    goodfrm ,
    rejfrm  ,
    freepar ,
    aic     ,
    bic     ,
    sdnr    ,
    chifact ,
    photonSNR
    ]

'''
for i in range(len(allparams)):
    nmissing = 4 - nfit
    for j in range(nmissing):
        allparams[i].append(0)
'''

# remove all errors from fixed parameters
# also make all 0.0 just 0
for param in allparams:
    for j in range(nfit):
        try:
            if param[j][-6:] == " & 0.0":
                param[j] = param[j][:-6]
        except:
            continue
        if param[j] == '0.0':
            param[j] = '0'

# generate table with 3 columns
def gentable3(offset):
    partable ="""
\\colhead{{Parameter}}                                                  &   \mctc{{{ch1wavelen:^25}}}   &   \mctc{{{ch2wavelen:^25}}}   &   \mctc{{{ch3wavelen:^25}}}      \\\\

Array Position (\math{{\\bar{{x}}}}, pix)                               &   \mctc{{{ch1xpos:^18.4}}}   &   \mctc{{{ch2xpos:^18.4}}}   &   \mctc{{{ch3xpos:^18.4}}}  \\\\
Array Position (\math{{\\bar{{y}}}}, pix)                               &   \mctc{{{ch1ypos:^18.4}}}   &   \mctc{{{ch2ypos:^18.4}}}   &   \mctc{{{ch3ypos:^18.4}}}  \\\\
Position Consistency\\tablenotemark{{a}} (\math{{\delta\sb{{x}}}}, pix) &   \mctc{{{ch1pcosx:^18}}}   &   \mctc{{{ch2pcosx:^18}}}   &   \mctc{{{ch3pcosx:^18}}}   \\\\
Position Consistency\\tablenotemark{{a}} (\math{{\delta\sb{{y}}}}, pix) &   \mctc{{{ch1pcosy:^18}}}   &   \mctc{{{ch2pcosy:^18}}}   &   \mctc{{{ch3pcosy:^18}}}   \\\\
Aperture Size (pix)                                                     &   \mctc{{{ch1photap:^25}}}   &   \mctc{{{ch2photap:^25}}}   &   \mctc{{{ch3photap:^25}}}     \\\\
Sky Annulus Inner Radius (pix)                                          &   \mctc{{{ch1skyin:^25}}}   &   \mctc{{{ch2skyin:^25}}}   &   \mctc{{{ch3skyin:^25}}}      \\\\
Sky Annulus Outer Radius (pix)                                          &   \mctc{{{ch1skyout:^25}}}   &   \mctc{{{ch2skyout:^25}}}   &   \mctc{{{ch3skyout:^25}}}      \\\\
System Flux \\math{{F\sb{{s}}}} (\\micro Jy)                            &   {ch1flux:^25}   &   {ch2flux:^25}   &   {ch3flux:^25}      \\\\
Transit Midpoint (BJD\\sb{{UTC}} - 2,450,000)                           &   {ch1midtimeutc:^25}   &   {ch2midtimeutc:^25}   &   {ch3midtimeutc:^25}     \\\\
Transit Midpoint (BJD\sb{{TDB}} - 2,450,000)                            &   {ch1midtimetdb:^25}   &   {ch2midtimetdb:^25}   &   {ch3midtimetdb:^25}    \\\\
\math{{R\\sb{{p}}/R\\sb{{\star}}}}                                      &   {ch1rprs:^25}   &   {ch2rprs:^25}   &   {ch3rprs:^25}      \\\\
cos i                                                                   &   \mctc{{{ch1cosi:^25}}}   &   \mctc{{{ch2cosi:^25}}}   &   \mctc{{{ch3cosi:^25}}}      \\\\
\math{{a/R\\sb{{\\star}}}}                                              &   \mctc{{{ch1ars:^25}}}   &   \mctc{{{ch2ars:^25}}}   &   \mctc{{{ch3ars:^25}}}      \\\\
Ramp Name                                                               &   \mctc{{{ch1ramp:^25}}}   &   \mctc{{{ch2ramp:^25}}}   &   \mctc{{{ch3ramp:^25}}}      \\\\
Ramp, \math{{r\sb{{0}}}}                                                &   {ch1r0:^25}   &   {ch2r0:^25}   &   {ch3r0:^25}   \\\\
Ramp, \math{{r\sb{{1}}}}                                                &   {ch1r1:^25}   &   {ch2r1:^25}   &   {ch3r1:^25}   \\\\
Ramp, \math{{r\sb{{2}}}}                                                &   {ch1r2:^25}   &   {ch2r2:^25}   &   {ch3r2:^25}   \\\\
Ramp, \math{{r\sb{{3}}}}                                                &   {ch1r3:^25}   &   {ch2r3:^25}   &   {ch3r3:^25}   \\\\
Ramp, \math{{r\sb{{4}}}}                                                &   {ch1r4:^25}   &   {ch2r4:^25}   &   {ch3r4:^25}   \\\\
Ramp, \math{{r\sb{{5}}}}                                                &   {ch1r5:^25}   &   {ch2r5:^25}   &   {ch3r5:^25}   \\\\
Ramp, \math{{r\sb{{6}}}}                                                &   {ch1r6:^25}   &   {ch2r6:^25}   &   {ch3r6:^25}   \\\\
Ramp, \math{{r\sb{{7}}}}                                                &   {ch1r7:^25}   &   {ch2r7:^25}   &   {ch3r7:^25}   \\\\
Ramp, \math{{t\sb{{0}}}}                                                &   {ch1t0:^25}   &   {ch2t0:^25}   &   {ch3t0:^25}   \\\\
BLISS Map (\math{{M(x,y)}})                                             &   \mctc{{{ch1isbliss:^18}}}   &   \mctc{{{ch2isbliss:^18}}}   &   \mctc{{{ch3isbliss:^18}}}    \\\\
Minimum Number of Points Per Bin                                        &   \mctc{{{ch1minnumpts:^18}}}   &   \mctc{{{ch2minnumpts:^18}}}   &   \mctc{{{ch3minnumpts:^18}}}   \\\\
Total Frames                                                            &   \mctc{{{ch1totfrm:^18}}}   &   \mctc{{{ch2totfrm:^18}}}   &   \mctc{{{ch3totfrm:^18}}}     \\\\
Frames Used                                                             &   \mctc{{{ch1goodfrm:^18}}}   &   \mctc{{{ch2goodfrm:^18}}}   &   \mctc{{{ch3goodfrm:^18}}}    \\\\
Rejected Frames (\\%)                                                   &   \mctc{{{ch1rejfrm:^18}}}   &   \mctc{{{ch2rejfrm:^18}}}   &   \mctc{{{ch3rejfrm:^18}}}   \\\\
Free Parameters                                                         &   \mctc{{{ch1freepar:^18}}}   &   \mctc{{{ch2freepar:^18}}}   &   \mctc{{{ch3freepar:^18}}}    \\\\
Number of Data Points in Fit                                            &   \mctc{{{ch1N:^18}}}   &   \mctc{{{ch2N:^18}}}   &   \mctc{{{ch3N:^18}}}    \\\\
AIC Value                                                               &   \mctc{{{ch1aic:^18}}}   &   \mctc{{{ch2aic:^18}}}   &   \mctc{{{ch3aic:^18}}}    \\\\
BIC Value                                                               &   \mctc{{{ch1bic:^18}}}   &   \mctc{{{ch2bic:^18}}}   &   \mctc{{{ch3bic:^18}}}   \\\\
SDNR                                                                    &   \mctc{{{ch1sdnr:^18}}}   &   \mctc{{{ch2sdnr:^18}}}   &   \mctc{{{ch3sdnr:^18}}}   \\\\
Uncertainty Scaling Factor                                              &   \mctc{{{ch1chifact:^18}}}   &   \mctc{{{ch2chifact:^18}}}   &   \mctc{{{ch3chifact:^18}}}  \\\\
Percentage of Photon-Limited S/N (\\%)                                  &   \mctc{{{ch1photonSNR:^18}}}   &   \mctc{{{ch2photonSNR:^18}}}   &   \mctc{{{ch3photonSNR:^18}}}  \\\\
""".format(
    ch1wavelen = wavelen[0+offset], # wavelength (um)
    ch2wavelen = wavelen[1+offset],
    ch3wavelen = wavelen[2+offset],
    ch1xpos    = xpos[0+offset],    # x position (pix)
    ch2xpos    = xpos[1+offset],
    ch3xpos    = xpos[2+offset],
    ch1ypos    = ypos[0+offset],    # y position (pix)     
    ch2ypos    = ypos[1+offset],
    ch3ypos    = ypos[2+offset],
    ch1pcosx   = pcosx[0+offset],   # x consistency (pix)
    ch2pcosx   = pcosx[1+offset],
    ch3pcosx   = pcosx[2+offset],
    ch1pcosy   = pcosy[0+offset],   # y consistency (pix)
    ch2pcosy   = pcosy[1+offset],
    ch3pcosy   = pcosy[2+offset],
    ch1photap  = photap[0+offset],  # aperture radius (pix)
    ch2photap  = photap[1+offset],
    ch3photap  = photap[2+offset],
    ch1skyin   = skyin[0+offset],   # sky annulus inner radius (pix)
    ch2skyin   = skyin[1+offset],
    ch3skyin   = skyin[2+offset],
    ch1skyout  = skyout[0+offset],  # sky annulus outer radius (pix)
    ch2skyout  = skyout[1+offset],
    ch3skyout  = skyout[2+offset],
    ch1flux    = flux[0+offset],    # system flux (uJy)
    ch2flux    = flux[1+offset],
    ch3flux    = flux[2+offset],
    ch1midtimeutc = midtimeutc[0+offset],    # midpoint (BJD_utc)
    ch2midtimeutc = midtimeutc[1+offset],
    ch3midtimeutc = midtimeutc[2+offset],
    ch1midtimetdb = midtimetdb[0+offset],    # midpoint (BJD_tdb)
    ch2midtimetdb = midtimetdb[1+offset],
    ch3midtimetdb = midtimetdb[2+offset],
    ch1rprs   = rprs[0+offset],      # duration (hrs)
    ch2rprs   = rprs[1+offset],  
    ch3rprs   = rprs[2+offset],  
    ch1cosi     = cosi[0+offset],        # ingress time (hrs)
    ch2cosi     = cosi[1+offset],    
    ch3cosi     = cosi[2+offset],    
    ch1ars     = ars[0+offset],        # egress time (hrs)
    ch2ars     = ars[1+offset],        
    ch3ars     = ars[2+offset],    
    ch1ramp    = ramp[0+offset],       # ramp name
    ch2ramp    = ramp[1+offset],   
    ch3ramp    = ramp[2+offset], 
    ch1r0      = r0[0+offset],      # ramp r0
    ch2r0      = r0[1+offset],  
    ch3r0      = r0[2+offset],   
    ch1r1      = r1[0+offset],      # ramp r1
    ch2r1      = r1[1+offset],  
    ch3r1      = r1[2+offset],   
    ch1r2      = r2[0+offset],      # ramp r2
    ch2r2      = r2[1+offset],  
    ch3r2      = r2[2+offset],   
    ch1r3      = r3[0+offset],      # ramp r3
    ch2r3      = r3[1+offset],  
    ch3r3      = r3[2+offset], 
    ch1r4      = r4[0+offset],      # ramp r4
    ch2r4      = r4[1+offset],  
    ch3r4      = r4[2+offset], 
    ch1r5      = r5[0+offset],      # ramp r5
    ch2r5      = r5[1+offset],  
    ch3r5      = r5[2+offset], 
    ch1r6      = r6[0+offset],      # ramp r6
    ch2r6      = r6[1+offset],  
    ch3r6      = r6[2+offset], 
    ch1r7      = r7[0+offset],      # ramp r7
    ch2r7      = r7[1+offset],  
    ch3r7      = r7[2+offset], 
    ch1t0      = t0[0+offset],      # ramp t0
    ch2t0      = t0[1+offset],  
    ch3t0      = t0[2+offset], 
    ch1isbliss  = isbliss[0+offset],            # use BLISS map
    ch2isbliss  = isbliss[1+offset], 
    ch3isbliss  = isbliss[2+offset], 
    ch1minnumpts  = minnumpts[0+offset],            # Min # of points/bin
    ch2minnumpts  = minnumpts[1+offset], 
    ch3minnumpts  = minnumpts[2+offset], 
    ch1totfrm  = totfrm[0+offset],            # total frames
    ch2totfrm  = totfrm[1+offset], 
    ch3totfrm  = totfrm[2+offset], 
    ch1goodfrm = goodfrm[0+offset],           # good frames
    ch2goodfrm = goodfrm[1+offset],
    ch3goodfrm = goodfrm[2+offset],
    ch1rejfrm  = rejfrm[0+offset],            # rejected frames
    ch2rejfrm  = rejfrm[1+offset], 
    ch3rejfrm  = rejfrm[2+offset], 
    ch1freepar = freepar[0+offset],           # free parameters
    ch2freepar = freepar[1+offset],
    ch3freepar = freepar[2+offset],
    ch1N       = N[0+offset],                 # number of data points
    ch2N       = N[1+offset],
    ch3N       = N[2+offset],
    ch1aic     = aic[0+offset],               # aic
    ch2aic     = aic[1+offset],    
    ch3aic     = aic[2+offset],     
    ch1bic     = bic[0+offset],               # bic
    ch2bic     = bic[1+offset],    
    ch3bic     = bic[2+offset],    
    ch1sdnr    = sdnr[0+offset],              # sdnr
    ch2sdnr    = sdnr[1+offset],   
    ch3sdnr    = sdnr[2+offset],    
    ch1chifact = chifact[0+offset],           # chisq rescale factor
    ch2chifact = chifact[1+offset],
    ch3chifact = chifact[2+offset],
    ch1photonSNR  = photonSNR[0+offset],            # photon-limited SNR
    ch2photonSNR  = photonSNR[1+offset], 
    ch3photonSNR  = photonSNR[2+offset])
    return partable

# generate table with 4 columns
def gentable4(offset):
    partable ="""
\\colhead{{Parameter}}                           &   {ch1wavelen:^25}   &   {ch2wavelen:^25}   &   {ch3wavelen:^25}   &   {ch4wavelen:^25}   \\\\

Array Position (x, pix)                       &   {ch1xpos:^25}   &   {ch2xpos:^25}   &   {ch3xpos:^25}   &   {ch4xpos:^25}   \\\\
Array Position (y, pix)                       &   {ch1ypos:^25}   &   {ch2ypos:^25}   &   {ch3ypos:^25}   &   {ch4ypos:^25}   \\\\
Position Consistency (x, pix)                 &   {ch1pcosx:^25}   &   {ch2pcosx:^25}   &   {ch3pcosx:^25}   &   {ch4pcosx:^25}   \\\\
Position Consistency (y, pix)                 &   {ch1pcosy:^25}   &   {ch2pcosy:^25}   &   {ch3pcosy:^25}   &   {ch4pcosy:^25}   \\\\
Aperture Size (pix)                           &   {ch1photap:^25}   &   {ch2photap:^25}   &   {ch3photap:^25}   &   {ch4photap:^25}   \\\\
Sky Annulus Inner Radius (pix)                &   {ch1skyin:^25}   &   {ch2skyin:^25}   &   {ch3skyin:^25}   &   {ch4skyin:^25}   \\\\
Sky Annulus Outer Radius (pix)                &   {ch1skyout:^25}   &   {ch2skyout:^25}   &   {ch3skyout:^25}   &   {ch4skyout:^25}   \\\\
System Flux (uJy)                             &   {ch1flux:^25}   &   {ch2flux:^25}   &   {ch3flux:^25}   &   {ch4flux:^25}   \\\\
Eclipse Depth (\\%)                            &   {ch1depth:^25}   &   {ch2depth:^25}   &   {ch3depth:^25}   &   {ch4depth:^25}   \\\\
Brightness Temperature (K)                    &   {ch1temp:^25}   &   {ch2temp:^25}   &   {ch3temp:^25}   &   {ch4temp:^25}   \\\\
Eclipse Midpoint (orbits)                     &   {ch1midpt:^25}   &   {ch2midpt:^25}   &   {ch3midpt:^25}   &   {ch4midpt:^25}   \\\\
Eclipse Midpoint (BJD\\sb{UTC} - 2,450,000)    &   {ch1midtimeutc:^25}   &   {ch2midtimeutc:^25}   &   {ch3midtimeutc:^25}   &   {ch4midtimeutc:^25}   \\\\
Eclipse Midpoint (BJD\\sb{TDB} - 2,450,000)    &   {ch1midtimetdb:^25}   &   {ch2midtimetdb:^25}   &   {ch3midtimetdb:^25}   &   {ch4midtimetdb:^25}   \\\\
Eclipse Duration (hrs)                        &   {ch1width:^25}   &   {ch2width:^25}   &   {ch3width:^25}   &   {ch4width:^25}   \\\\
Ingress Time (hrs)                            &   {ch1t12:^25}   &   {ch2t12:^25}   &   {ch3t12:^25}   &   {ch4t12:^25}   \\\\
Egress Time (hrs)                             &   {ch1t34:^25}   &   {ch2t34:^25}   &   {ch3t34:^25}   &   {ch4t34:^25}   \\\\
Ramp Name                                     &   {ch1ramp:^25}   &   {ch2ramp:^25}   &   {ch3ramp:^25}   &   {ch4ramp:^25}   \\\\
Ramp, \math{{r\sb{{0}}}}                      &   {ch1r0:^25}   &   {ch2r0:^25}   &   {ch3r0:^25}   &   {ch4r0:^25}   \\\\
Ramp, \math{{r\sb{{1}}}}                      &   {ch1r1:^25}   &   {ch2r1:^25}   &   {ch3r1:^25}   &   {ch4r1:^25}   \\\\
Ramp, \math{{r\sb{{2}}}}                      &   {ch1r2:^25}   &   {ch2r2:^25}   &   {ch3r2:^25}   &   {ch4r2:^25}   \\\\
Ramp, \math{{r\sb{{3}}}}                      &   {ch1r3:^25}   &   {ch2r3:^25}   &   {ch3r3:^25}   &   {ch4r3:^25}   \\\\
Ramp, \math{{r\sb{{4}}}}                      &   {ch1r4:^25}   &   {ch2r4:^25}   &   {ch3r4:^25}   &   {ch4r4:^25}   \\\\
Ramp, \math{{r\sb{{5}}}}                      &   {ch1r5:^25}   &   {ch2r5:^25}   &   {ch3r5:^25}   &   {ch4r5:^25}   \\\\
Ramp, \math{{r\sb{{6}}}}                      &   {ch1r6:^25}   &   {ch2r6:^25}   &   {ch3r6:^25}   &   {ch4r6:^25}   \\\\
Ramp, \math{{r\sb{{7}}}}                      &   {ch1r7:^25}   &   {ch2r7:^25}   &   {ch3r7:^25}   &   {ch4r7:^25}   \\\\
Ramp, \math{{t\sb{{0}}}}                      &   {ch1t0:^25}   &   {ch2t0:^25}   &   {ch3t0:^25}   &   {ch4t0:^25}   \\\\
BLISS Map (\math{{M(x,y)}})                     &   {ch1isbliss:^25}   &   {ch2isbliss:^25}   &   {ch3isbliss:^25}   &   {ch4isbliss:^25}   \\\\
Minimum Number of Points Per Bin              &   {ch1minnumpts:^25}   &   {ch2minnumpts:^25}   &   {ch3minnumpts:^25}   &   {ch4minnumpts:^25}   \\\\
Total Frames                                  &   {ch1totfrm:^25}   &   {ch2totfrm:^25}   &   {ch3totfrm:^25}   &   {ch4totfrm:^25}   \\\\
Frames Used                                   &   {ch1goodfrm:^25}   &   {ch2goodfrm:^25}   &   {ch3goodfrm:^25}   &   {ch4goodfrm:^25}   \\\\
Rejected Frames (\\%)                          &   {ch1rejfrm:^25}   &   {ch2rejfrm:^25}   &   {ch3rejfrm:^25}   &   {ch4rejfrm:^25}   \\\\
Free Parameters                               &   {ch1freepar:^25}   &   {ch2freepar:^25}   &   {ch3freepar:^25}   &   {ch4freepar:^25}   \\\\
Number of Data Points in Fit                  &   {ch1N:^25}   &   {ch2N:^25}   &   {ch3N:^25}   &   {ch4N:^25}   \\\\
AIC Value                                     &   {ch1aic:^25}   &   {ch2aic:^25}   &   {ch3aic:^25}   &   {ch4aic:^25}   \\\\
BIC Value                                     &   {ch1bic:^25}   &   {ch2bic:^25}   &   {ch3bic:^25}   &   {ch4bic:^25}   \\\\
Standard Deviation of Normalized Residuals    &   {ch1sdnr:^25}   &   {ch2sdnr:^25}   &   {ch3sdnr:^25}   &   {ch4sdnr:^25}   \\\\
Chi-squared Normalization Factor              &   {ch1chifact:^25}   &   {ch2chifact:^25}   &   {ch3chifact:^25}   &   {ch4chifact:^25}   \\\\
Percentage of Photon-Limited S/N (\\%)         &   {ch1photonSNR:^25}   &   {ch2photonSNR:^25}   &   {ch3photonSNR:^25}   &   {ch4photonSNR:^25}   \\\\
""".format(
    ch1wavelen = wavelen[0+offset], # wavelength (um)
    ch2wavelen = wavelen[1+offset],
    ch3wavelen = wavelen[2+offset],
    ch4wavelen = wavelen[3+offset],
    ch1xpos    = xpos[0+offset],    # x position (pix)
    ch2xpos    = xpos[1+offset],
    ch3xpos    = xpos[2+offset],
    ch4xpos    = xpos[3+offset],
    ch1ypos    = ypos[0+offset],    # y position (pix)     
    ch2ypos    = ypos[1+offset],
    ch3ypos    = ypos[2+offset],
    ch4ypos    = ypos[3+offset],
    ch1pcosx   = pcosx[0+offset],   # x consistency (pix)
    ch2pcosx   = pcosx[1+offset],
    ch3pcosx   = pcosx[2+offset],
    ch4pcosx   = pcosx[3+offset],
    ch1pcosy   = pcosy[0+offset],   # y consistency (pix)
    ch2pcosy   = pcosy[1+offset],
    ch3pcosy   = pcosy[2+offset],
    ch4pcosy   = pcosy[3+offset],
    ch1photap  = photap[0+offset],  # aperture radius (pix)
    ch2photap  = photap[1+offset],
    ch3photap  = photap[2+offset],
    ch4photap  = photap[3+offset],
    ch1skyin   = skyin[0+offset],   # sky annulus inner radius (pix)
    ch2skyin   = skyin[1+offset],
    ch3skyin   = skyin[2+offset],
    ch4skyin   = skyin[3+offset],
    ch1skyout  = skyout[0+offset],  # sky annulus outer radius (pix)
    ch2skyout  = skyout[1+offset],
    ch3skyout  = skyout[2+offset],
    ch4skyout  = skyout[3+offset],
    ch1flux    = flux[0+offset],    # system flux (uJy)
    ch2flux    = flux[1+offset],
    ch3flux    = flux[2+offset],
    ch4flux    = flux[3+offset],
    ch1depth   = depth[0+offset],   # eclipse depth (%)
    ch2depth   = depth[1+offset],
    ch3depth   = depth[2+offset],
    ch4depth   = depth[3+offset],
    ch1temp    = temp[0+offset],       # brightness temperature (K)
    ch2temp    = temp[1+offset],   
    ch3temp    = temp[2+offset],   
    ch4temp    = temp[3+offset],   
    ch1midpt   = midpt[0+offset],      # midpoint (phase)     
    ch2midpt   = midpt[1+offset], 
    ch3midpt   = midpt[2+offset],  
    ch4midpt   = midpt[3+offset],  
    ch1midtimeutc = midtimeutc[0+offset],    # midpoint (BJD_utc)
    ch2midtimeutc = midtimeutc[1+offset],
    ch3midtimeutc = midtimeutc[2+offset],
    ch4midtimeutc = midtimeutc[3+offset],
    ch1midtimetdb = midtimetdb[0+offset],    # midpoint (BJD_tdb)
    ch2midtimetdb = midtimetdb[1+offset],
    ch3midtimetdb = midtimetdb[2+offset],
    ch4midtimetdb = midtimetdb[3+offset],
    ch1width   = width[0+offset],      # duration (hrs)
    ch2width   = width[1+offset],  
    ch3width   = width[2+offset],  
    ch4width   = width[3+offset],  
    ch1t12     = t12[0+offset],        # ingress time (hrs)
    ch2t12     = t12[1+offset],    
    ch3t12     = t12[2+offset],    
    ch4t12     = t12[3+offset],    
    ch1t34     = t34[0+offset],        # egress time (hrs)
    ch2t34     = t34[1+offset],        
    ch3t34     = t34[2+offset],    
    ch4t34     = t34[3+offset],    
    ch1ramp    = ramp[0+offset],       # ramp name
    ch2ramp    = ramp[1+offset],   
    ch3ramp    = ramp[2+offset],   
    ch4ramp    = ramp[3+offset],   
    ch1r0      = r0[0+offset],      # ramp r0
    ch2r0      = r0[1+offset],  
    ch3r0      = r0[2+offset],  
    ch4r0      = r0[3+offset],   
    ch1r1      = r1[0+offset],      # ramp r1
    ch2r1      = r1[1+offset],  
    ch3r1      = r1[2+offset],    
    ch4r1      = r1[3+offset], 
    ch1r2      = r2[0+offset],      # ramp r2
    ch2r2      = r2[1+offset],  
    ch3r2      = r2[2+offset],    
    ch4r2      = r2[3+offset], 
    ch1r3      = r3[0+offset],      # ramp r3
    ch2r3      = r3[1+offset],  
    ch3r3      = r3[2+offset],  
    ch4r3      = r3[3+offset], 
    ch1r4      = r4[0+offset],      # ramp r4
    ch2r4      = r4[1+offset],  
    ch3r4      = r4[2+offset],  
    ch4r4      = r4[3+offset], 
    ch1r5      = r5[0+offset],      # ramp r5
    ch2r5      = r5[1+offset],  
    ch3r5      = r5[2+offset],  
    ch4r5      = r5[3+offset], 
    ch1r6      = r6[0+offset],      # ramp r6
    ch2r6      = r6[1+offset],  
    ch3r6      = r6[2+offset],  
    ch4r6      = r6[3+offset], 
    ch1r7      = r7[0+offset],      # ramp r7
    ch2r7      = r7[1+offset],  
    ch3r7      = r7[2+offset],  
    ch4r7      = r7[3+offset], 
    ch1t0      = t0[0+offset],      # ramp t0
    ch2t0      = t0[1+offset],  
    ch3t0      = t0[2+offset],  
    ch4t0      = t0[3+offset],
    ch1isbliss  = isbliss[0+offset],            # use BLISS map
    ch2isbliss  = isbliss[1+offset], 
    ch3isbliss  = isbliss[2+offset], 
    ch4isbliss  = isbliss[3+offset], 
    ch1minnumpts  = minnumpts[0+offset],            # Min # of points/bin
    ch2minnumpts  = minnumpts[1+offset], 
    ch3minnumpts  = minnumpts[2+offset], 
    ch4minnumpts  = minnumpts[3+offset], 
    ch1totfrm  = totfrm[0+offset],            # total frames
    ch2totfrm  = totfrm[1+offset], 
    ch3totfrm  = totfrm[2+offset], 
    ch4totfrm  = totfrm[3+offset], 
    ch1goodfrm = goodfrm[0+offset],           # good frames
    ch2goodfrm = goodfrm[1+offset],
    ch3goodfrm = goodfrm[2+offset],
    ch4goodfrm = goodfrm[3+offset],
    ch1rejfrm  = rejfrm[0+offset],            # rejected frames
    ch2rejfrm  = rejfrm[1+offset], 
    ch3rejfrm  = rejfrm[2+offset], 
    ch4rejfrm  = rejfrm[3+offset], 
    ch1freepar = freepar[0+offset],           # free parameters
    ch2freepar = freepar[1+offset],
    ch3freepar = freepar[2+offset],
    ch4freepar = freepar[3+offset],
    ch1N       = N[0+offset],                 # number of data points
    ch2N       = N[1+offset],
    ch3N       = N[2+offset],
    ch4N       = N[3+offset],
    ch1aic     = aic[0+offset],               # aic
    ch2aic     = aic[1+offset],    
    ch3aic     = aic[2+offset],    
    ch4aic     = aic[3+offset],    
    ch1bic     = bic[0+offset],               # bic
    ch2bic     = bic[1+offset],    
    ch3bic     = bic[2+offset],    
    ch4bic     = bic[3+offset],    
    ch1sdnr    = sdnr[0+offset],              # sdnr
    ch2sdnr    = sdnr[1+offset],   
    ch3sdnr    = sdnr[2+offset],   
    ch4sdnr    = sdnr[3+offset],   
    ch1chifact = chifact[0+offset],           # chisq rescale factor
    ch2chifact = chifact[1+offset],
    ch3chifact = chifact[2+offset],
    ch4chifact = chifact[3+offset],
    ch1photonSNR  = photonSNR[0+offset],            # photon-limited SNR
    ch2photonSNR  = photonSNR[1+offset], 
    ch3photonSNR  = photonSNR[2+offset], 
    ch4photonSNR  = photonSNR[3+offset])
    return partable
    
# generate table with 5 columns
def gentable5(offset):
    partable ="""
\\colhead{{Parameter}}                                              &   \mctc{{{ch1wavelen:^18}}}   &   \mctc{{{ch2wavelen:^18}}}   &   \mctc{{{ch3wavelen:^18}}}   &   \mctc{{{ch4wavelen:^18}}}   &   \mctc{{{ch5wavelen:^18}}}   \\\\

Array Position (\math{{\\bar{{x}}}}, pix)                             &   \mctc{{{ch1xpos:^18.4}}}   &   \mctc{{{ch2xpos:^18.4}}}   &   \mctc{{{ch3xpos:^18.4}}}   &   \mctc{{{ch4xpos:^18.4}}}   &   \mctc{{{ch5xpos:^18.4}}}   \\\\
Array Position (\math{{\\bar{{y}}}}, pix)                             &   \mctc{{{ch1ypos:^18.4}}}   &   \mctc{{{ch2ypos:^18.4}}}   &   \mctc{{{ch3ypos:^18.4}}}   &   \mctc{{{ch4ypos:^18.4}}}   &   \mctc{{{ch5ypos:^18.4}}}   \\\\
Position Consistency\\tablenotemark{{a}} (\math{{\delta\sb{{x}}}}, pix) &   \mctc{{{ch1pcosx:^18}}}   &   \mctc{{{ch2pcosx:^18}}}   &   \mctc{{{ch3pcosx:^18}}}   &   \mctc{{{ch4pcosx:^18}}}   &   \mctc{{{ch5pcosx:^18}}}   \\\\
Position Consistency\\tablenotemark{{a}} (\math{{\delta\sb{{y}}}}, pix) &   \mctc{{{ch1pcosy:^18}}}   &   \mctc{{{ch2pcosy:^18}}}   &   \mctc{{{ch3pcosy:^18}}}   &   \mctc{{{ch4pcosy:^18}}}   &   \mctc{{{ch5pcosy:^18}}}   \\\\
Aperture Size (pix)                                              &   \mctc{{{ch1photap:^18}}}   &   \mctc{{{ch2photap:^18}}}   &   \mctc{{{ch3photap:^18}}}   &   \mctc{{{ch4photap:^18}}}   &   \mctc{{{ch5photap:^18}}}   \\\\
Sky Annulus Inner Radius (pix)                                   &   \mctc{{{ch1skyin:^18}}}   &   \mctc{{{ch2skyin:^18}}}   &   \mctc{{{ch3skyin:^18}}}   &   \mctc{{{ch4skyin:^18}}}   &   \mctc{{{ch5skyin:^18}}}   \\\\
Sky Annulus Outer Radius (pix)                                   &   \mctc{{{ch1skyout:^18}}}   &   \mctc{{{ch2skyout:^18}}}   &   \mctc{{{ch3skyout:^18}}}   &   \mctc{{{ch4skyout:^18}}}   &   \mctc{{{ch5skyout:^18}}}   \\\\
System Flux \math{{F\sb{{s}}}} (\micro Jy)                           &   {ch1flux:^25}   &   {ch2flux:^25}   &   {ch3flux:^25}   &   {ch4flux:^25}   &   {ch5flux:^25}   \\\\
Eclipse Depth (\\%)                                               &   {ch1depth:^25}   &   {ch2depth:^25}   &   {ch3depth:^25}   &   {ch4depth:^25}   &   {ch5depth:^25}   \\\\
Brightness Temperature (K)                                       &   {ch1temp:^25}   &   {ch2temp:^25}   &   {ch3temp:^25}   &   {ch4temp:^25}   &   {ch5temp:^25}   \\\\
Eclipse Midpoint (orbits)                                        &   {ch1midpt:^25}   &   {ch2midpt:^25}   &   {ch3midpt:^25}   &   {ch4midpt:^25}   &   {ch5midpt:^25}   \\\\
Eclipse Midpoint (BJD\\sb{{UTC}} - 2,450,000)                       &   {ch1midtimeutc:^25}   &   {ch2midtimeutc:^25}   &   {ch3midtimeutc:^25}   &   {ch4midtimeutc:^25}   &   {ch5midtimeutc:^25}   \\\\
Eclipse Midpoint (BJD\\sb{{TDB}} - 2,450,000)                       &   {ch1midtimetdb:^25}   &   {ch2midtimetdb:^25}   &   {ch3midtimetdb:^25}   &   {ch4midtimetdb:^25}   &   {ch5midtimetdb:^25}   \\\\
Eclipse Duration (\math{{t\sb{{\\rm 4-1}}}}, hrs)                     &   {ch1width:^25}   &   {ch2width:^25}   &   {ch3width:^25}   &   {ch4width:^25}   &   {ch5width:^25}   \\\\
Ingress Time (\math{{t\sb{{\\rm 2-1}}}}, hrs)                         &   {ch1t12:^25}   &   {ch2t12:^25}   &   {ch3t12:^25}   &   {ch4t12:^25}   &   {ch5t12:^25}   \\\\
Egress Time (\math{{t\sb{{\\rm 4-3}}}}, hrs)                          &   {ch1t34:^25}   &   {ch2t34:^25}   &   {ch3t34:^25}   &   {ch4t34:^25}   &   {ch5t34:^25}   \\\\
Ramp Name                                                        &   \mctc{{{ch1ramp:^18}}}   &   \mctc{{{ch2ramp:^18}}}   &   \mctc{{{ch3ramp:^18}}}   &   \mctc{{{ch4ramp:^18}}}   &   \mctc{{{ch5ramp:^18}}}   \\\\
Ramp, \math{{r\sb{{0}}}}                                                &   {ch1r0:^25}   &   {ch2r0:^25}   &   {ch3r0:^25}   &   {ch4r0:^25}   &   {ch5r0:^25}   \\\\
Ramp, \math{{r\sb{{1}}}}                                                &   {ch1r1:^25}   &   {ch2r1:^25}   &   {ch3r1:^25}   &   {ch4r1:^25}   &   {ch5r1:^25}   \\\\
Ramp, \math{{r\sb{{2}}}}                                                &   {ch1r2:^25}   &   {ch2r2:^25}   &   {ch3r2:^25}   &   {ch4r2:^25}   &   {ch5r2:^25}   \\\\
Ramp, \math{{r\sb{{3}}}}                                                &   {ch1r3:^25}   &   {ch2r3:^25}   &   {ch3r3:^25}   &   {ch4r3:^25}   &   {ch5r3:^25}   \\\\
Ramp, \math{{r\sb{{4}}}}                                                &   {ch1r4:^25}   &   {ch2r4:^25}   &   {ch3r4:^25}   &   {ch4r4:^25}   &   {ch5r4:^25}   \\\\
Ramp, \math{{r\sb{{5}}}}                                                &   {ch1r5:^25}   &   {ch2r5:^25}   &   {ch3r5:^25}   &   {ch4r5:^25}   &   {ch5r5:^25}   \\\\
Ramp, \math{{r\sb{{6}}}}                                                &   {ch1r6:^25}   &   {ch2r6:^25}   &   {ch3r6:^25}   &   {ch4r6:^25}   &   {ch5r6:^25}   \\\\
Ramp, \math{{r\sb{{7}}}}                                                &   {ch1r7:^25}   &   {ch2r7:^25}   &   {ch3r7:^25}   &   {ch4r7:^25}   &   {ch5r7:^25}   \\\\
Ramp, \math{{t\sb{{0}}}}                                                &   {ch1t0:^25}   &   {ch2t0:^25}   &   {ch3t0:^25}   &   {ch4t0:^25}   &   {ch5t0:^25}   \\\\
BLISS Map (\math{{M(x,y)}})                                        &   \mctc{{{ch1isbliss:^18}}}   &   \mctc{{{ch2isbliss:^18}}}   &   \mctc{{{ch3isbliss:^18}}}   &   \mctc{{{ch4isbliss:^18}}}   &   \mctc{{{ch5isbliss:^18}}}  \\\\
Minimum Number of Points Per Bin                                 &   \mctc{{{ch1minnumpts:^18}}}   &   \mctc{{{ch2minnumpts:^18}}}   &   \mctc{{{ch3minnumpts:^18}}}   &   \mctc{{{ch4minnumpts:^18}}}   &   \mctc{{{ch5minnumpts:^18}}}   \\\\
Total Frames                                                     &   \mctc{{{ch1totfrm:^18}}}   &   \mctc{{{ch2totfrm:^18}}}   &   \mctc{{{ch3totfrm:^18}}}   &   \mctc{{{ch4totfrm:^18}}}   &   \mctc{{{ch5totfrm:^18}}}   \\\\
Frames Used                                                      &   \mctc{{{ch1goodfrm:^18}}}   &   \mctc{{{ch2goodfrm:^18}}}   &   \mctc{{{ch3goodfrm:^18}}}   &   \mctc{{{ch4goodfrm:^18}}}   &   \mctc{{{ch5goodfrm:^18}}}   \\\\
Rejected Frames (\\%)                                             &   \mctc{{{ch1rejfrm:^18}}}   &   \mctc{{{ch2rejfrm:^18}}}   &   \mctc{{{ch3rejfrm:^18}}}   &   \mctc{{{ch4rejfrm:^18}}}   &   \mctc{{{ch5rejfrm:^18}}}   \\\\
Free Parameters                                                  &   \mctc{{{ch1freepar:^18}}}   &   \mctc{{{ch2freepar:^18}}}   &   \mctc{{{ch3freepar:^18}}}   &   \mctc{{{ch4freepar:^18}}}   &   \mctc{{{ch5freepar:^18}}}   \\\\
Number of Data Points in Fit                                     &   \mctc{{{ch1N:^18}}}   &   \mctc{{{ch2N:^18}}}   &   \mctc{{{ch3N:^18}}}   &   \mctc{{{ch4N:^18}}}   &   \mctc{{{ch5N:^18}}}   \\\\
AIC Value                                                        &   \mctc{{{ch1aic:^18}}}   &   \mctc{{{ch2aic:^18}}}   &   \mctc{{{ch3aic:^18}}}   &   \mctc{{{ch4aic:^18}}}   &   \mctc{{{ch5aic:^18}}}   \\\\
BIC Value                                                        &   \mctc{{{ch1bic:^18}}}   &   \mctc{{{ch2bic:^18}}}   &   \mctc{{{ch3bic:^18}}}   &   \mctc{{{ch4bic:^18}}}   &   \mctc{{{ch5bic:^18}}}   \\\\
SDNR                                                             &   \mctc{{{ch1sdnr:^18}}}   &   \mctc{{{ch2sdnr:^18}}}   &   \mctc{{{ch3sdnr:^18}}}   &   \mctc{{{ch4sdnr:^18}}}   &   \mctc{{{ch5sdnr:^18}}}   \\\\
Uncertainty Scaling Factor                                       &   \mctc{{{ch1chifact:^18}}}   &   \mctc{{{ch2chifact:^18}}}   &   \mctc{{{ch3chifact:^18}}}   &   \mctc{{{ch4chifact:^18}}}   &   \mctc{{{ch5chifact:^18}}}   \\\\
Percentage of Photon-Limited S/N (\\%)                            &   \mctc{{{ch1photonSNR:^18}}}   &   \mctc{{{ch2photonSNR:^18}}}   &   \mctc{{{ch3photonSNR:^18}}}   &   \mctc{{{ch4photonSNR:^18}}}   &   \mctc{{{ch5photonSNR:^18}}}   \\\\
""".format(
    ch1wavelen = wavelen[0+offset], # wavelength (um)
    ch2wavelen = wavelen[1+offset],
    ch3wavelen = wavelen[2+offset],
    ch4wavelen = wavelen[3+offset],
    ch5wavelen = wavelen[4+offset],
    ch1xpos    = xpos[0+offset],    # x position (pix)
    ch2xpos    = xpos[1+offset],
    ch3xpos    = xpos[2+offset],
    ch4xpos    = xpos[3+offset],
    ch5xpos    = xpos[4+offset],
    ch1ypos    = ypos[0+offset],    # y position (pix)     
    ch2ypos    = ypos[1+offset],
    ch3ypos    = ypos[2+offset],
    ch4ypos    = ypos[3+offset],
    ch5ypos    = ypos[4+offset],
    ch1pcosx   = pcosx[0+offset],   # x consistency (pix)
    ch2pcosx   = pcosx[1+offset],
    ch3pcosx   = pcosx[2+offset],
    ch4pcosx   = pcosx[3+offset],
    ch5pcosx   = pcosx[4+offset],
    ch1pcosy   = pcosy[0+offset],   # y consistency (pix)
    ch2pcosy   = pcosy[1+offset],
    ch3pcosy   = pcosy[2+offset],
    ch4pcosy   = pcosy[3+offset],
    ch5pcosy   = pcosy[4+offset],
    ch1photap  = photap[0+offset],  # aperture radius (pix)
    ch2photap  = photap[1+offset],
    ch3photap  = photap[2+offset],
    ch4photap  = photap[3+offset],
    ch5photap  = photap[4+offset],
    ch1skyin   = skyin[0+offset],   # sky annulus inner radius (pix)
    ch2skyin   = skyin[1+offset],
    ch3skyin   = skyin[2+offset],
    ch4skyin   = skyin[3+offset],
    ch5skyin   = skyin[4+offset],
    ch1skyout  = skyout[0+offset],  # sky annulus outer radius (pix)
    ch2skyout  = skyout[1+offset],
    ch3skyout  = skyout[2+offset],
    ch4skyout  = skyout[3+offset],
    ch5skyout  = skyout[4+offset],
    ch1flux    = flux[0+offset],    # system flux (uJy)
    ch2flux    = flux[1+offset],
    ch3flux    = flux[2+offset],
    ch4flux    = flux[3+offset],
    ch5flux    = flux[4+offset],
    ch1depth   = depth[0+offset],   # eclipse depth (%)
    ch2depth   = depth[1+offset],
    ch3depth   = depth[2+offset],
    ch4depth   = depth[3+offset],
    ch5depth   = depth[4+offset],
    ch1temp    = temp[0+offset],       # brightness temperature (K)
    ch2temp    = temp[1+offset],   
    ch3temp    = temp[2+offset],   
    ch4temp    = temp[3+offset],   
    ch5temp    = temp[4+offset],  
    ch1midpt   = midpt[0+offset],      # midpoint (phase)     
    ch2midpt   = midpt[1+offset], 
    ch3midpt   = midpt[2+offset],  
    ch4midpt   = midpt[3+offset],  
    ch5midpt   = midpt[4+offset],  
    ch1midtimeutc = midtimeutc[0+offset],    # midpoint (BJD_utc)
    ch2midtimeutc = midtimeutc[1+offset],
    ch3midtimeutc = midtimeutc[2+offset],
    ch4midtimeutc = midtimeutc[3+offset],
    ch5midtimeutc = midtimeutc[4+offset],
    ch1midtimetdb = midtimetdb[0+offset],    # midpoint (BJD_tdb)
    ch2midtimetdb = midtimetdb[1+offset],
    ch3midtimetdb = midtimetdb[2+offset],
    ch4midtimetdb = midtimetdb[3+offset],
    ch5midtimetdb = midtimetdb[4+offset],
    ch1width   = width[0+offset],      # duration (hrs)
    ch2width   = width[1+offset],  
    ch3width   = width[2+offset],  
    ch4width   = width[3+offset],  
    ch5width   = width[4+offset], 
    ch1t12     = t12[0+offset],        # ingress time (hrs)
    ch2t12     = t12[1+offset],    
    ch3t12     = t12[2+offset],    
    ch4t12     = t12[3+offset],     
    ch5t12     = t12[4+offset],   
    ch1t34     = t34[0+offset],        # egress time (hrs)
    ch2t34     = t34[1+offset],        
    ch3t34     = t34[2+offset],    
    ch4t34     = t34[3+offset],      
    ch5t34     = t34[4+offset],  
    ch1ramp    = ramp[0+offset],       # ramp name
    ch2ramp    = ramp[1+offset],   
    ch3ramp    = ramp[2+offset],   
    ch4ramp    = ramp[3+offset],    
    ch5ramp    = ramp[4+offset],   
    ch1r0      = r0[0+offset],      # ramp r0
    ch2r0      = r0[1+offset],  
    ch3r0      = r0[2+offset],  
    ch4r0      = r0[3+offset],   
    ch5r0      = r0[4+offset], 
    ch1r1      = r1[0+offset],      # ramp r1
    ch2r1      = r1[1+offset],  
    ch3r1      = r1[2+offset],    
    ch4r1      = r1[3+offset], 
    ch5r1      = r1[4+offset], 
    ch1r2      = r2[0+offset],      # ramp r2
    ch2r2      = r2[1+offset],  
    ch3r2      = r2[2+offset],    
    ch4r2      = r2[3+offset], 
    ch5r2      = r2[4+offset], 
    ch1r3      = r3[0+offset],      # ramp r3
    ch2r3      = r3[1+offset],  
    ch3r3      = r3[2+offset],  
    ch4r3      = r3[3+offset], 
    ch5r3      = r3[4+offset], 
    ch1r4      = r4[0+offset],      # ramp r4
    ch2r4      = r4[1+offset],  
    ch3r4      = r4[2+offset],  
    ch4r4      = r4[3+offset], 
    ch5r4      = r4[4+offset], 
    ch1r5      = r5[0+offset],      # ramp r5
    ch2r5      = r5[1+offset],  
    ch3r5      = r5[2+offset],  
    ch4r5      = r5[3+offset], 
    ch5r5      = r5[4+offset], 
    ch1r6      = r6[0+offset],      # ramp r6
    ch2r6      = r6[1+offset],  
    ch3r6      = r6[2+offset],  
    ch4r6      = r6[3+offset], 
    ch5r6      = r6[4+offset], 
    ch1r7      = r7[0+offset],      # ramp r7
    ch2r7      = r7[1+offset],  
    ch3r7      = r7[2+offset],  
    ch4r7      = r7[3+offset], 
    ch5r7      = r7[4+offset], 
    ch1t0      = t0[0+offset],      # ramp t0
    ch2t0      = t0[1+offset],  
    ch3t0      = t0[2+offset],  
    ch4t0      = t0[3+offset],
    ch5t0      = t0[4+offset],
    ch1isbliss  = isbliss[0+offset],            # use BLISS map
    ch2isbliss  = isbliss[1+offset], 
    ch3isbliss  = isbliss[2+offset], 
    ch4isbliss  = isbliss[3+offset], 
    ch5isbliss  = isbliss[4+offset], 
    ch1minnumpts  = minnumpts[0+offset],            # Min # of points/bin
    ch2minnumpts  = minnumpts[1+offset], 
    ch3minnumpts  = minnumpts[2+offset], 
    ch4minnumpts  = minnumpts[3+offset], 
    ch5minnumpts  = minnumpts[4+offset], 
    ch1totfrm  = totfrm[0+offset],            # total frames
    ch2totfrm  = totfrm[1+offset], 
    ch3totfrm  = totfrm[2+offset], 
    ch4totfrm  = totfrm[3+offset], 
    ch5totfrm  = totfrm[4+offset], 
    ch1goodfrm = goodfrm[0+offset],           # good frames
    ch2goodfrm = goodfrm[1+offset],
    ch3goodfrm = goodfrm[2+offset],
    ch4goodfrm = goodfrm[3+offset],
    ch5goodfrm = goodfrm[4+offset],
    ch1rejfrm  = rejfrm[0+offset],            # rejected frames
    ch2rejfrm  = rejfrm[1+offset], 
    ch3rejfrm  = rejfrm[2+offset], 
    ch4rejfrm  = rejfrm[3+offset], 
    ch5rejfrm  = rejfrm[4+offset], 
    ch1freepar = freepar[0+offset],           # free parameters
    ch2freepar = freepar[1+offset],
    ch3freepar = freepar[2+offset],
    ch4freepar = freepar[3+offset],
    ch5freepar = freepar[4+offset],
    ch1N       = N[0+offset],                 # number of data points
    ch2N       = N[1+offset],
    ch3N       = N[2+offset],
    ch4N       = N[3+offset],
    ch5N       = N[4+offset],
    ch1aic     = aic[0+offset],               # aic
    ch2aic     = aic[1+offset],    
    ch3aic     = aic[2+offset],    
    ch4aic     = aic[3+offset],     
    ch5aic     = aic[4+offset], 
    ch1bic     = bic[0+offset],               # bic
    ch2bic     = bic[1+offset],    
    ch3bic     = bic[2+offset],    
    ch4bic     = bic[3+offset],     
    ch5bic     = bic[4+offset],   
    ch1sdnr    = sdnr[0+offset],              # sdnr
    ch2sdnr    = sdnr[1+offset],   
    ch3sdnr    = sdnr[2+offset],   
    ch4sdnr    = sdnr[3+offset],   
    ch5sdnr    = sdnr[4+offset],  
    ch1chifact = chifact[0+offset],           # chisq rescale factor
    ch2chifact = chifact[1+offset],
    ch3chifact = chifact[2+offset],
    ch4chifact = chifact[3+offset],
    ch5chifact = chifact[4+offset],
    ch1photonSNR  = photonSNR[0+offset],            # photon-limited SNR
    ch2photonSNR  = photonSNR[1+offset], 
    ch3photonSNR  = photonSNR[2+offset], 
    ch4photonSNR  = photonSNR[3+offset],
    ch5photonSNR  = photonSNR[4+offset])
    return partable
    
# generate table with 6 columns
def gentable6(offset):
    partable ="""
\\colhead{{Parameter}}                                              &   \mctc{{{ch1wavelen:^18}}}   &   \mctc{{{ch2wavelen:^18}}}   &   \mctc{{{ch3wavelen:^18}}}   &   \mctc{{{ch4wavelen:^18}}}   &   \mctc{{{ch5wavelen:^18}}}   &   \mctc{{{ch6wavelen:^18}}}   \\\\

Array Position (\math{{\\bar{{x}}}}, pix)                             &   \mctc{{{ch1xpos:^18.4}}}   &   \mctc{{{ch2xpos:^18.4}}}   &   \mctc{{{ch3xpos:^18.4}}}   &   \mctc{{{ch4xpos:^18.4}}}   &   \mctc{{{ch5xpos:^18.4}}}   &   \mctc{{{ch6xpos:^18.4}}}   \\\\
Array Position (\math{{\\bar{{y}}}}, pix)                             &   \mctc{{{ch1ypos:^18.4}}}   &   \mctc{{{ch2ypos:^18.4}}}   &   \mctc{{{ch3ypos:^18.4}}}   &   \mctc{{{ch4ypos:^18.4}}}   &   \mctc{{{ch5ypos:^18.4}}}    &   \mctc{{{ch6ypos:^18.4}}}  \\\\
Position Consistency\\tablenotemark{{a}} (\math{{\delta\sb{{x}}}}, pix) &   \mctc{{{ch1pcosx:^18}}}   &   \mctc{{{ch2pcosx:^18}}}   &   \mctc{{{ch3pcosx:^18}}}   &   \mctc{{{ch4pcosx:^18}}}   &   \mctc{{{ch5pcosx:^18}}}   &   \mctc{{{ch6pcosx:^18}}}   \\\\
Position Consistency\\tablenotemark{{a}} (\math{{\delta\sb{{y}}}}, pix) &   \mctc{{{ch1pcosy:^18}}}   &   \mctc{{{ch2pcosy:^18}}}   &   \mctc{{{ch3pcosy:^18}}}   &   \mctc{{{ch4pcosy:^18}}}   &   \mctc{{{ch5pcosy:^18}}}   &   \mctc{{{ch6pcosy:^18}}}   \\\\
Aperture Size (pix)                                              &   \mctc{{{ch1photap:^18}}}   &   \mctc{{{ch2photap:^18}}}   &   \mctc{{{ch3photap:^18}}}   &   \mctc{{{ch4photap:^18}}}   &   \mctc{{{ch5photap:^18}}}   &   \mctc{{{ch6photap:^18}}}   \\\\
Sky Annulus Inner Radius (pix)                                   &   \mctc{{{ch1skyin:^18}}}   &   \mctc{{{ch2skyin:^18}}}   &   \mctc{{{ch3skyin:^18}}}   &   \mctc{{{ch4skyin:^18}}}   &   \mctc{{{ch5skyin:^18}}}   &   \mctc{{{ch6skyin:^18}}}   \\\\
Sky Annulus Outer Radius (pix)                                   &   \mctc{{{ch1skyout:^18}}}   &   \mctc{{{ch2skyout:^18}}}   &   \mctc{{{ch3skyout:^18}}}   &   \mctc{{{ch4skyout:^18}}}   &   \mctc{{{ch5skyout:^18}}}   &   \mctc{{{ch6skyout:^18}}}   \\\\
System Flux \math{{F\sb{{s}}}} (\micro Jy)                           &   {ch1flux:^25}   &   {ch2flux:^25}   &   {ch3flux:^25}   &   {ch4flux:^25}   &   {ch5flux:^25}   &   {ch6flux:^25}   \\\\
Eclipse Depth (\\%)                                               &   {ch1depth:^25}   &   {ch2depth:^25}   &   {ch3depth:^25}   &   {ch4depth:^25}   &   {ch5depth:^25}   &   {ch6depth:^25}   \\\\
Brightness Temperature (K)                                       &   {ch1temp:^25}   &   {ch2temp:^25}   &   {ch3temp:^25}   &   {ch4temp:^25}   &   {ch5temp:^25}   &   {ch6temp:^25}   \\\\
Eclipse Midpoint (orbits)                                        &   {ch1midpt:^25}   &   {ch2midpt:^25}   &   {ch3midpt:^25}   &   {ch4midpt:^25}   &   {ch5midpt:^25}   &   {ch6midpt:^25}   \\\\
Eclipse Midpoint (BJD\\sb{{UTC}} - 2,450,000)                       &   {ch1midtimeutc:^25}   &   {ch2midtimeutc:^25}   &   {ch3midtimeutc:^25}   &   {ch4midtimeutc:^25}   &   {ch5midtimeutc:^25}   &   {ch6midtimeutc:^25}   \\\\
Eclipse Midpoint (BJD\\sb{{TDB}} - 2,450,000)                       &   {ch1midtimetdb:^25}   &   {ch2midtimetdb:^25}   &   {ch3midtimetdb:^25}   &   {ch4midtimetdb:^25}   &   {ch5midtimetdb:^25}   &   {ch6midtimetdb:^25}   \\\\
Eclipse Duration (\math{{t\sb{{\\rm 4-1}}}}, hrs)                     &   {ch1width:^25}   &   {ch2width:^25}   &   {ch3width:^25}   &   {ch4width:^25}   &   {ch5width:^25}   &   {ch6width:^25}   \\\\
Ingress Time (\math{{t\sb{{\\rm 2-1}}}}, hrs)                         &   {ch1t12:^25}   &   {ch2t12:^25}   &   {ch3t12:^25}   &   {ch4t12:^25}   &   {ch5t12:^25}   &   {ch6t12:^25}   \\\\
Egress Time (\math{{t\sb{{\\rm 4-3}}}}, hrs)                          &   {ch1t34:^25}   &   {ch2t34:^25}   &   {ch3t34:^25}   &   {ch4t34:^25}   &   {ch5t34:^25}   &   {ch6t34:^25}   \\\\
Ramp Name                                                        &   \mctc{{{ch1ramp:^18}}}   &   \mctc{{{ch2ramp:^18}}}   &   \mctc{{{ch3ramp:^18}}}   &   \mctc{{{ch4ramp:^18}}}   &   \mctc{{{ch5ramp:^18}}}   &   \mctc{{{ch6ramp:^18}}}   \\\\
Ramp, \math{{r\sb{{0}}}}                                                &   {ch1r0:^25}   &   {ch2r0:^25}   &   {ch3r0:^25}   &   {ch4r0:^25}   &   {ch5r0:^25}   &   {ch6r0:^25}   \\\\
Ramp, \math{{r\sb{{1}}}}                                                &   {ch1r1:^25}   &   {ch2r1:^25}   &   {ch3r1:^25}   &   {ch4r1:^25}   &   {ch5r1:^25}   &   {ch6r1:^25}   \\\\
Ramp, \math{{r\sb{{2}}}}                                                &   {ch1r2:^25}   &   {ch2r2:^25}   &   {ch3r2:^25}   &   {ch4r2:^25}   &   {ch5r2:^25}   &   {ch6r2:^25}   \\\\
Ramp, \math{{r\sb{{3}}}}                                                &   {ch1r3:^25}   &   {ch2r3:^25}   &   {ch3r3:^25}   &   {ch4r3:^25}   &   {ch5r3:^25}   &   {ch6r3:^25}   \\\\
Ramp, \math{{r\sb{{4}}}}                                                &   {ch1r4:^25}   &   {ch2r4:^25}   &   {ch3r4:^25}   &   {ch4r4:^25}   &   {ch5r4:^25}   &   {ch6r4:^25}   \\\\
Ramp, \math{{r\sb{{5}}}}                                                &   {ch1r5:^25}   &   {ch2r5:^25}   &   {ch3r5:^25}   &   {ch4r5:^25}   &   {ch5r5:^25}   &   {ch6r5:^25}   \\\\
Ramp, \math{{r\sb{{6}}}}                                                &   {ch1r6:^25}   &   {ch2r6:^25}   &   {ch3r6:^25}   &   {ch4r6:^25}   &   {ch5r6:^25}   &   {ch6r6:^25}   \\\\
Ramp, \math{{r\sb{{7}}}}                                                &   {ch1r7:^25}   &   {ch2r7:^25}   &   {ch3r7:^25}   &   {ch4r7:^25}   &   {ch5r7:^25}   &   {ch6r7:^25}   \\\\
Ramp, \math{{t\sb{{0}}}}                                                &   {ch1t0:^25}   &   {ch2t0:^25}   &   {ch3t0:^25}   &   {ch4t0:^25}   &   {ch5t0:^25}   &   {ch6t0:^25}   \\\\
BLISS Map (\math{{M(x,y)}})                                        &   \mctc{{{ch1isbliss:^18}}}   &   \mctc{{{ch2isbliss:^18}}}   &   \mctc{{{ch3isbliss:^18}}}   &   \mctc{{{ch4isbliss:^18}}}   &   \mctc{{{ch5isbliss:^18}}}   &   \mctc{{{ch6isbliss:^18}}}   \\\\
Minimum Number of Points Per Bin                                 &   \mctc{{{ch1minnumpts:^18}}}   &   \mctc{{{ch2minnumpts:^18}}}   &   \mctc{{{ch3minnumpts:^18}}}   &   \mctc{{{ch4minnumpts:^18}}}   &   \mctc{{{ch5minnumpts:^18}}}   &   \mctc{{{ch6minnumpts:^18}}}   \\\\
Total Frames                                                     &   \mctc{{{ch1totfrm:^18}}}   &   \mctc{{{ch2totfrm:^18}}}   &   \mctc{{{ch3totfrm:^18}}}   &   \mctc{{{ch4totfrm:^18}}}   &   \mctc{{{ch5totfrm:^18}}}   &   \mctc{{{ch6totfrm:^18}}}   \\\\
Frames Used                                                      &   \mctc{{{ch1goodfrm:^18}}}   &   \mctc{{{ch2goodfrm:^18}}}   &   \mctc{{{ch3goodfrm:^18}}}   &   \mctc{{{ch4goodfrm:^18}}}   &   \mctc{{{ch5goodfrm:^18}}}   &   \mctc{{{ch6goodfrm:^18}}}   \\\\
Rejected Frames (\\%)                                             &   \mctc{{{ch1rejfrm:^18}}}   &   \mctc{{{ch2rejfrm:^18}}}   &   \mctc{{{ch3rejfrm:^18}}}   &   \mctc{{{ch4rejfrm:^18}}}   &   \mctc{{{ch5rejfrm:^18}}}   &   \mctc{{{ch6rejfrm:^18}}}   \\\\
Free Parameters                                                  &   \mctc{{{ch1freepar:^18}}}   &   \mctc{{{ch2freepar:^18}}}   &   \mctc{{{ch3freepar:^18}}}   &   \mctc{{{ch4freepar:^18}}}   &   \mctc{{{ch5freepar:^18}}}   &   \mctc{{{ch6freepar:^18}}}   \\\\
Number of Data Points in Fit                                     &   \mctc{{{ch1N:^18}}}   &   \mctc{{{ch2N:^18}}}   &   \mctc{{{ch3N:^18}}}   &   \mctc{{{ch4N:^18}}}   &   \mctc{{{ch5N:^18}}}   &   \mctc{{{ch6N:^18}}}   \\\\
AIC Value                                                        &   \mctc{{{ch1aic:^18}}}   &   \mctc{{{ch2aic:^18}}}   &   \mctc{{{ch3aic:^18}}}   &   \mctc{{{ch4aic:^18}}}   &   \mctc{{{ch5aic:^18}}}   &   \mctc{{{ch6aic:^18}}}   \\\\
BIC Value                                                        &   \mctc{{{ch1bic:^18}}}   &   \mctc{{{ch2bic:^18}}}   &   \mctc{{{ch3bic:^18}}}   &   \mctc{{{ch4bic:^18}}}   &   \mctc{{{ch5bic:^18}}}   &   \mctc{{{ch6bic:^18}}}   \\\\
SDNR                                                             &   \mctc{{{ch1sdnr:^18}}}   &   \mctc{{{ch2sdnr:^18}}}   &   \mctc{{{ch3sdnr:^18}}}   &   \mctc{{{ch4sdnr:^18}}}   &   \mctc{{{ch5sdnr:^18}}}   &   \mctc{{{ch6sdnr:^18}}}   \\\\
Uncertainty Scaling Factor                                       &   \mctc{{{ch1chifact:^18}}}   &   \mctc{{{ch2chifact:^18}}}   &   \mctc{{{ch3chifact:^18}}}   &   \mctc{{{ch4chifact:^18}}}   &   \mctc{{{ch5chifact:^18}}}   &   \mctc{{{ch6chifact:^18}}}  \\\\
Percentage of Photon-Limited S/N (\\%)                            &   \mctc{{{ch1photonSNR:^18}}}   &   \mctc{{{ch2photonSNR:^18}}}   &   \mctc{{{ch3photonSNR:^18}}}   &   \mctc{{{ch4photonSNR:^18}}}   &   \mctc{{{ch5photonSNR:^18}}}   &   \mctc{{{ch6photonSNR:^18}}}   \\\\
""".format(
    ch1wavelen = wavelen[0+offset], # wavelength (um)
    ch2wavelen = wavelen[1+offset],
    ch3wavelen = wavelen[2+offset],
    ch4wavelen = wavelen[3+offset],
    ch5wavelen = wavelen[4+offset],
    ch6wavelen = wavelen[5+offset],
    ch1xpos    = xpos[0+offset],    # x position (pix)
    ch2xpos    = xpos[1+offset],
    ch3xpos    = xpos[2+offset],
    ch4xpos    = xpos[3+offset],
    ch5xpos    = xpos[4+offset],
    ch6xpos    = xpos[5+offset],
    ch1ypos    = ypos[0+offset],    # y position (pix)     
    ch2ypos    = ypos[1+offset],
    ch3ypos    = ypos[2+offset],
    ch4ypos    = ypos[3+offset],
    ch5ypos    = ypos[4+offset],
    ch6ypos    = ypos[5+offset],
    ch1pcosx   = pcosx[0+offset],   # x consistency (pix)
    ch2pcosx   = pcosx[1+offset],
    ch3pcosx   = pcosx[2+offset],
    ch4pcosx   = pcosx[3+offset],
    ch5pcosx   = pcosx[4+offset],
    ch6pcosx   = pcosx[5+offset],
    ch1pcosy   = pcosy[0+offset],   # y consistency (pix)
    ch2pcosy   = pcosy[1+offset],
    ch3pcosy   = pcosy[2+offset],
    ch4pcosy   = pcosy[3+offset],
    ch5pcosy   = pcosy[4+offset],
    ch6pcosy   = pcosy[5+offset],
    ch1photap  = photap[0+offset],  # aperture radius (pix)
    ch2photap  = photap[1+offset],
    ch3photap  = photap[2+offset],
    ch4photap  = photap[3+offset],
    ch5photap  = photap[4+offset],
    ch6photap  = photap[5+offset],
    ch1skyin   = skyin[0+offset],   # sky annulus inner radius (pix)
    ch2skyin   = skyin[1+offset],
    ch3skyin   = skyin[2+offset],
    ch4skyin   = skyin[3+offset],
    ch5skyin   = skyin[4+offset],
    ch6skyin   = skyin[5+offset],
    ch1skyout  = skyout[0+offset],  # sky annulus outer radius (pix)
    ch2skyout  = skyout[1+offset],
    ch3skyout  = skyout[2+offset],
    ch4skyout  = skyout[3+offset],
    ch5skyout  = skyout[4+offset],
    ch6skyout  = skyout[5+offset],
    ch1flux    = flux[0+offset],    # system flux (uJy)
    ch2flux    = flux[1+offset],
    ch3flux    = flux[2+offset],
    ch4flux    = flux[3+offset],
    ch5flux    = flux[4+offset],
    ch6flux    = flux[5+offset],
    ch1depth   = depth[0+offset],   # eclipse depth (%)
    ch2depth   = depth[1+offset],
    ch3depth   = depth[2+offset],
    ch4depth   = depth[3+offset],
    ch5depth   = depth[4+offset],
    ch6depth   = depth[5+offset],
    ch1temp    = temp[0+offset],       # brightness temperature (K)
    ch2temp    = temp[1+offset],   
    ch3temp    = temp[2+offset],   
    ch4temp    = temp[3+offset],   
    ch5temp    = temp[4+offset],  
    ch6temp    = temp[5+offset], 
    ch1midpt   = midpt[0+offset],      # midpoint (phase)     
    ch2midpt   = midpt[1+offset], 
    ch3midpt   = midpt[2+offset],  
    ch4midpt   = midpt[3+offset],  
    ch5midpt   = midpt[4+offset],  
    ch6midpt   = midpt[5+offset],  
    ch1midtimeutc = midtimeutc[0+offset],    # midpoint (BJD_utc)
    ch2midtimeutc = midtimeutc[1+offset],
    ch3midtimeutc = midtimeutc[2+offset],
    ch4midtimeutc = midtimeutc[3+offset],
    ch5midtimeutc = midtimeutc[4+offset],
    ch6midtimeutc = midtimeutc[5+offset],
    ch1midtimetdb = midtimetdb[0+offset],    # midpoint (BJD_tdb)
    ch2midtimetdb = midtimetdb[1+offset],
    ch3midtimetdb = midtimetdb[2+offset],
    ch4midtimetdb = midtimetdb[3+offset],
    ch5midtimetdb = midtimetdb[4+offset],
    ch6midtimetdb = midtimetdb[5+offset],
    ch1width   = width[0+offset],      # duration (hrs)
    ch2width   = width[1+offset],  
    ch3width   = width[2+offset],  
    ch4width   = width[3+offset],  
    ch5width   = width[4+offset], 
    ch6width   = width[5+offset], 
    ch1t12     = t12[0+offset],        # ingress time (hrs)
    ch2t12     = t12[1+offset],    
    ch3t12     = t12[2+offset],    
    ch4t12     = t12[3+offset],     
    ch5t12     = t12[4+offset],   
    ch6t12     = t12[5+offset],   
    ch1t34     = t34[0+offset],        # egress time (hrs)
    ch2t34     = t34[1+offset],        
    ch3t34     = t34[2+offset],    
    ch4t34     = t34[3+offset],      
    ch5t34     = t34[4+offset],     
    ch6t34     = t34[5+offset],  
    ch1ramp    = ramp[0+offset],       # ramp name
    ch2ramp    = ramp[1+offset],   
    ch3ramp    = ramp[2+offset],   
    ch4ramp    = ramp[3+offset],    
    ch5ramp    = ramp[4+offset],   
    ch6ramp    = ramp[5+offset],  
    ch1r0      = r0[0+offset],      # ramp r0
    ch2r0      = r0[1+offset],  
    ch3r0      = r0[2+offset],  
    ch4r0      = r0[3+offset],   
    ch5r0      = r0[4+offset], 
    ch6r0      = r0[5+offset], 
    ch1r1      = r1[0+offset],      # ramp r1
    ch2r1      = r1[1+offset],  
    ch3r1      = r1[2+offset],    
    ch4r1      = r1[3+offset], 
    ch5r1      = r1[4+offset], 
    ch6r1      = r1[5+offset], 
    ch1r2      = r2[0+offset],      # ramp r2
    ch2r2      = r2[1+offset],  
    ch3r2      = r2[2+offset],    
    ch4r2      = r2[3+offset], 
    ch5r2      = r2[4+offset], 
    ch6r2      = r2[5+offset], 
    ch1r3      = r3[0+offset],      # ramp r3
    ch2r3      = r3[1+offset],  
    ch3r3      = r3[2+offset],  
    ch4r3      = r3[3+offset], 
    ch5r3      = r3[4+offset], 
    ch6r3      = r3[5+offset], 
    ch1r4      = r4[0+offset],      # ramp r4
    ch2r4      = r4[1+offset],  
    ch3r4      = r4[2+offset],  
    ch4r4      = r4[3+offset], 
    ch5r4      = r4[4+offset], 
    ch6r4      = r4[5+offset], 
    ch1r5      = r5[0+offset],      # ramp r5
    ch2r5      = r5[1+offset],  
    ch3r5      = r5[2+offset],  
    ch4r5      = r5[3+offset], 
    ch5r5      = r5[4+offset], 
    ch6r5      = r5[5+offset], 
    ch1r6      = r6[0+offset],      # ramp r6
    ch2r6      = r6[1+offset],  
    ch3r6      = r6[2+offset],  
    ch4r6      = r6[3+offset], 
    ch5r6      = r6[4+offset], 
    ch6r6      = r6[5+offset], 
    ch1r7      = r7[0+offset],      # ramp r7
    ch2r7      = r7[1+offset],  
    ch3r7      = r7[2+offset],  
    ch4r7      = r7[3+offset], 
    ch5r7      = r7[4+offset], 
    ch6r7      = r7[5+offset], 
    ch1t0      = t0[0+offset],      # ramp t0
    ch2t0      = t0[1+offset],  
    ch3t0      = t0[2+offset],  
    ch4t0      = t0[3+offset],
    ch5t0      = t0[4+offset],
    ch6t0      = t0[5+offset], 
    ch1isbliss  = isbliss[0+offset],            # use BLISS map
    ch2isbliss  = isbliss[1+offset], 
    ch3isbliss  = isbliss[2+offset], 
    ch4isbliss  = isbliss[3+offset], 
    ch5isbliss  = isbliss[4+offset], 
    ch6isbliss  = isbliss[5+offset],
    ch1minnumpts  = minnumpts[0+offset],            # Min # of points/bin
    ch2minnumpts  = minnumpts[1+offset], 
    ch3minnumpts  = minnumpts[2+offset], 
    ch4minnumpts  = minnumpts[3+offset], 
    ch5minnumpts  = minnumpts[4+offset], 
    ch6minnumpts  = minnumpts[5+offset], 
    ch1totfrm  = totfrm[0+offset],            # total frames
    ch2totfrm  = totfrm[1+offset], 
    ch3totfrm  = totfrm[2+offset], 
    ch4totfrm  = totfrm[3+offset], 
    ch5totfrm  = totfrm[4+offset], 
    ch6totfrm  = totfrm[5+offset], 
    ch1goodfrm = goodfrm[0+offset],           # good frames
    ch2goodfrm = goodfrm[1+offset],
    ch3goodfrm = goodfrm[2+offset],
    ch4goodfrm = goodfrm[3+offset],
    ch5goodfrm = goodfrm[4+offset],
    ch6goodfrm = goodfrm[5+offset],
    ch1rejfrm  = rejfrm[0+offset],            # rejected frames
    ch2rejfrm  = rejfrm[1+offset], 
    ch3rejfrm  = rejfrm[2+offset], 
    ch4rejfrm  = rejfrm[3+offset], 
    ch5rejfrm  = rejfrm[4+offset], 
    ch6rejfrm  = rejfrm[5+offset], 
    ch1freepar = freepar[0+offset],           # free parameters
    ch2freepar = freepar[1+offset],
    ch3freepar = freepar[2+offset],
    ch4freepar = freepar[3+offset],
    ch5freepar = freepar[4+offset],
    ch6freepar = freepar[5+offset],
    ch1N       = N[0+offset],                 # number of data points
    ch2N       = N[1+offset],
    ch3N       = N[2+offset],
    ch4N       = N[3+offset],
    ch5N       = N[4+offset],
    ch6N       = N[5+offset],
    ch1aic     = aic[0+offset],               # aic
    ch2aic     = aic[1+offset],    
    ch3aic     = aic[2+offset],    
    ch4aic     = aic[3+offset],     
    ch5aic     = aic[4+offset],    
    ch6aic     = aic[5+offset],
    ch1bic     = bic[0+offset],               # bic
    ch2bic     = bic[1+offset],    
    ch3bic     = bic[2+offset],    
    ch4bic     = bic[3+offset],     
    ch5bic     = bic[4+offset],     
    ch6bic     = bic[5+offset],   
    ch1sdnr    = sdnr[0+offset],              # sdnr
    ch2sdnr    = sdnr[1+offset],   
    ch3sdnr    = sdnr[2+offset],   
    ch4sdnr    = sdnr[3+offset],   
    ch5sdnr    = sdnr[4+offset],  
    ch6sdnr    = sdnr[5+offset], 
    ch1chifact = chifact[0+offset],           # chisq rescale factor
    ch2chifact = chifact[1+offset],
    ch3chifact = chifact[2+offset],
    ch4chifact = chifact[3+offset],
    ch5chifact = chifact[4+offset],
    ch6chifact = chifact[5+offset],
    ch1photonSNR  = photonSNR[0+offset],            # photon-limited SNR
    ch2photonSNR  = photonSNR[1+offset], 
    ch3photonSNR  = photonSNR[2+offset], 
    ch4photonSNR  = photonSNR[3+offset],
    ch5photonSNR  = photonSNR[4+offset],
    ch6photonSNR  = photonSNR[5+offset])
    return partable

#ENTER MANUALLY
#Generate table with 4 columns
#Start at fit[offset]
#partable = gentable4(offset=0)
#print(partable1)

#Generate table with 5 columns
#Start at fit[offset]
partable1 = gentable5(offset=0)
#print(partable1)

#Generate table with 6 columns
#Start at fit[offset]
partable2 = gentable6(offset=5)
#print(partable2)

if fname:
    outfile = file(fname, 'w')
    outfile.write(partable1)
    outfile.write(partable2)
    outfile.close()


