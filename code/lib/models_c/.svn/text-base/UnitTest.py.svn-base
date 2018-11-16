#! /usr/bin/env python


from __future__ import print_function
import numpy as np
import models_c as mc
import models as mp
import os
os.environ['OMP_NUM_THREADS']='1'
import matplotlib.pyplot as plt

#Ramp parameters
t      = np.arange(0,1,0.0001)
t0     = 0.5
t1     = 0.5
t2     = 0.5
goal   = 1.0
m      = 20.
m1     = 10.
m2     = 100.
a      = 1.0
b      = 1.0
c      = 1.0
d      = 1.0
e      = 1.0
f      = 1.0
g      = 1.0
knots  = [1., 1., 1.]
#Eclipse parameters
midpt = 0.25
width = 0.1
depth = 0.01
t12   = 0.01
t34   = 0.01
flux  = 1.
#Transit parameters
rprs  = 0.3
cosi  = 0.01
ars   = 5.
p     = 0.5
c1    = 0.0
c2    = 0.0
c3    = 0.0
c4    = 0.25
#Intrapixel parameters
nobj     = 100
ysize    = 10
xsize    = 10
ymin     = 0
ymax     = ysize
xmin     = 0
xmax     = xsize
ygrid, ystep  = np.linspace(ymin, ymax, ysize, retstep=True)
xgrid, xstep  = np.linspace(xmin, xmax, xsize, retstep=True)
xygrid   = np.meshgrid(xgrid, ygrid)

ipparams         = np.array([0.0])
y                = np.random.normal(4.5, 1, nobj)
x                = np.random.normal(4.5, 1, nobj)
q                = np.random.randint(0, 4, nobj)
wbfipmask        = []
binloc           = np.zeros((2, nobj), dtype=int) - 1
numpts           = np.zeros((ysize, xsize), dtype=int)
binfluxmask      = np.zeros(nobj, dtype=int)
for i in range(ysize):
    wbftemp   = np.where(np.abs(y-ygrid[i]) < (ystep/2.))[0]
    for j in range(xsize):
        wbf = wbftemp[np.where((np.abs(x[wbftemp]-xgrid[j]) < (xstep/2.)))]
        wbfipmask.append(wbf)
        numpts[i,j] = len(wbf)
        if numpts[i,j] >= 1:
            binloc[0, wbf] = i*xsize + j
            binfluxmask[i*xsize + j] = 1

griddist       = np.ones((4, nobj))
for i in range(ysize-1):
    wherey = np.where(np.bitwise_and(y > ygrid[i  ], 
                                     y < ygrid[i+1]))[0]
    for j in range(xsize-1):
        wherexy = wherey[np.where(np.bitwise_and(x[wherey] > xgrid[j  ],
                                                 x[wherey] < xgrid[j+1]))[0]]
        if len(wherexy) > 0:
            binloc[1, wherexy] = gridpt = i*xsize + j
            #IF THERE ARE NO POINTS IN ONE OR MORE BINS...
            if (len(wbfipmask[gridpt        ]) < 1e-10) or \
               (len(wbfipmask[gridpt      +1]) < 1e-10) or \
               (len(wbfipmask[gridpt+xsize  ]) < 1e-10) or \
               (len(wbfipmask[gridpt+xsize+1]) < 1e-10):
                #SET griddist = NEAREST BIN (USE NEAREST NEIGHBOR INTERPOLATION)
                for loc in wherexy:
                    if   loc in wbfipmask[gridpt        ]:
                        griddist[0, loc] = 0
                        griddist[2, loc] = 0
                    elif loc in wbfipmask[gridpt      +1]:
                        griddist[0, loc] = 0
                        griddist[3, loc] = 0
                    elif loc in wbfipmask[gridpt+xsize  ]:
                        griddist[1, loc] = 0
                        griddist[2, loc] = 0
                    elif loc in wbfipmask[gridpt+xsize+1]:
                        griddist[1, loc] = 0
                        griddist[3, loc] = 0
            else:
                #CALCULATE griddist NORMALLY FOR BILINEAR INTERPOLATION
                griddist[0, wherexy] = np.array((y[wherexy]-ygrid[i])  /ystep)
                griddist[1, wherexy] = np.array((ygrid[i+1]-y[wherexy])/ystep)
                griddist[2, wherexy] = np.array((x[wherexy]-xgrid[j])  /xstep)
                griddist[3, wherexy] = np.array((xgrid[j+1]-x[wherexy])/xstep)

aplev            = np.random.normal(1,0.1,nobj)
etc              = np.ones(nobj)
kernel           = np.zeros((ysize, xsize))
tup1             = [0, 0, 0.0, 0.0]
issmoothing      = False

### bilinit()

posflux = [y, x, aplev, wbfipmask, binfluxmask, kernel, \
          tup1, binloc, griddist,  \
          xygrid[0].shape, issmoothing]

yp  =  mp.bilinint(ipparams, posflux, etc)
yc = mc.bilinint(ipparams, posflux, etc)

if sum(abs(yp-yc)) < 1e-10:
    print("Bilinint:   PASS")
else:
    print("Bilinint:   FAIL", sum(abs(yp-yc)))

### expramp()
rampparams = np.array((goal, m, t0))

yp  = mp.expramp(rampparams, t, [])
yc = mc.expramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Expramp:    PASS")
else:
    print("Expramp:    FAIL", sum(abs(yp-yc)))


### fallingexp()
# Variables 'goal', 'm', 't0' in expramp()

rampparams = np.array((goal, m, t0))

yp  = mp.fallingexp(rampparams, t, [])
yc = mc.fallingexp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Fallingexp: PASS")
else:
    print("Fallingexp: FAIL", sum(abs(yp-yc)))


### felramp()
# Variables 'goal', 'm', 't0' in fallingexp(), 't' in expramp()

rampparams = np.array((goal, m, t0, a, t1))

yp  = mp.felramp(rampparams, t, [])
yc = mc.felramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Felramp:    PASS")
else:
    print("Felramp:    FAIL", sum(abs(yp-yc)))


### linramp()
# Variables 't0' in expramp(), 't' in fallingexp()

rampparams = np.array((a, b, t0))

yp  = mp.linramp(rampparams, t, [])
yc = mc.linramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Linramp:    PASS")
else:
    print("Linramp:    FAIL", sum(abs(yp-yc)))


### llramp()
# Variables 't0' in expramp(), 't' in fallingexp()

rampparams = np.array((t0, a, b, c, t1))

yp  = mp.llramp(rampparams, t, [])
yc = mc.llramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Llramp:     PASS")
else:
    print("Llramp:     FAIL", sum(abs(yp-yc)))


### log4qramp()
# Variables 't0' in expramp(), 't' in fallingexp(), 't1' in llramp()

rampparams = np.array((t0, a, b, c, d, e, f, g, t1))

yp  = mp.log4qramp(rampparams, t, [])
yc = mc.log4qramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Log4qramp:  PASS")
else:
    print("Log4qramp:  FAIL", sum(abs(yp-yc)))


### logramp()
# Variables 't0' in expramp(), 't' in fallingexp()

rampparams = np.array((t0, a, b, c, d, e))

yp  = mp.logramp(rampparams, t, [])
yc = mc.logramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Logramp:    PASS")
else:
    print("Logramp:    FAIL", sum(abs(yp-yc)))


### lqramp()
# Variables 't0' in expramp(), 't' in fallingexp(), 't1' in llramp()

rampparams = np.array((t0, a, b, c, d, t1))

yp  = mp.lqramp(rampparams, t, [])
yc = mc.lqramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Lqramp:     PASS")
else:
    print("Lqramp:     FAIL", sum(abs(yp-yc)))


### mandelecl()

eclparams = np.array((midpt, width, depth, t12, t34, flux))

yp  = mp.mandelecl(eclparams, t, [])
yc = mc.mandelecl(eclparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Mandelecl:  PASS")
else:
    print("Mandelecl:  FAIL", sum(abs(yp-yc)))


### mandeltr()
# Variables 'midpt', 'flux', 't' in mandelecl()

params = np.array((midpt, rprs, cosi, ars, flux, p))

yp  = mp.mandeltr(params, t, [])
yc = mc.mandeltr(params, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Mandeltr:   PASS")
else:
    print ("Mandeltr:  FAIL", sum(abs(yp-yc)))


### nnint()
# Variables in 'posflux', 'ipparams' in bilinint()

yp  = mp.nnint(ipparams, posflux, etc)
yc = mc.nnint(ipparams, posflux, etc)
#print(yp-yc)

if sum(abs(yp-yc)) < 1e-10:
    print("Nnint:      PASS")
else:
    print("Nnint:      FAIL", sum(abs(yp-yc)))

### ortho()
'''
params   = np.arange(4.0)
invtrans = np.matrix(np.identity(4))
origin   = np.zeros(4) + 10

yp  = mp.orthoInvtrans(params, invtrans, origin)
yc = mc.orthoInvTrans(params, invtrans, origin)

if sum(abs(yp-yc)) < 1e-10:
    print("ortho:      PASS")
else:
    print("ortho:      FAIL", sum(abs(yp-yc)))
'''
### quadip()

ipparams = np.array((a, b, c, d, e, f))
position = np.array((y, x, q))

yp  = mp.quadip(ipparams, position, [])
yc = mc.quadip(ipparams, position, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Quadip:     PASS")
else:
    print("Quadip:     FAIL", sum(abs(yp-yc)))


### quadramp()
# Variables 'midpt', 'width', 'depth' in mandelecl(), 't0' in expramp(), 't' in fallingexp()

rampparams = np.array((a, b, c, t0))

yp  = mp.quadramp(rampparams, t, [])
yc = mc.quadramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Quadramp:   PASS")
else:
    print("Quadramp:   FAIL", sum(abs(yp-yc)))


### re2ramp()
# Variables 'goal', 't' in expramp()

rampparams = np.array((goal, a, m1, t1, b, m2, t2))

yp  = mp.re2ramp(rampparams, t, [])
yc = mc.re2ramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Re2ramp:    PASS")
else:
    print("Re2ramp:    FAIL", sum(abs(yp-yc)))


### relramp()
# Variables 'goal', 'm', 't0' in expramp(), 't' in fallingexp(), 't1' in llramp()

rampparams = np.array((goal, m, t0, a, b, t1))

yp  = mp.relramp(rampparams, t, [])
yc = mc.relramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Relramp:    PASS")
else:
    print("Relramp:    FAIL", sum(abs(yp-yc)))


### reqramp()
# Variables 'goal', 'm', 't0' in expramp(), 't' in fallingexp(), 't1' in llramp()

rampparams = np.array((goal, m, t0, a, b, c, t1))

yp  = mp.relramp(rampparams, t, [])
yc = mc.relramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Reqramp:    PASS")
else:
    print("Reqramp:    FAIL", sum(abs(yp-yc)))


### risingexp()
# Variables 'goal', 'm', 't0' in expramp(), 't' in fallingexp()

rampparams = np.array((goal, m, t0))

yp  = mp.risingexp(rampparams, t, [])
yc = mc.risingexp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Risingexp:  PASS")
else:
    print("Risingexp:  FAIL", sum(abs(yp-yc)))


### seramp()
rampparams = np.array((goal, m, m*t0, -1))

yp  = mp.seramp(rampparams, t, [])
yc = mc.seramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Seramp:     PASS")
else:
    print("Seramp:     FAIL", sum(abs(yp-yc)))

### selramp()
rampparams = np.array((goal, m, m*t0, a, t1, -1))

yp  = mp.selramp(rampparams, t, [])
yc = mc.selramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Selramp:    PASS")
else:
    print("Selramp:    FAIL", sum(abs(yp-yc)))

### seqramp()
rampparams = np.array((goal, m, m*t0, a, b, t1, -1))

yp  = mp.seqramp(rampparams, t, [])
yc = mc.seqramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Reqramp:    PASS")
else:
    print("Reqramp:    FAIL", sum(abs(yp-yc)))

### se2ramp()
rampparams = np.array((goal, m, m*t0, -1, m2, m2*t2, -1))

yp  = mp.se2ramp(rampparams, t, [])
yc = mc.se2ramp(rampparams, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Se2ramp:    PASS")
else:
    print("Se2ramp:    FAIL", sum(abs(yp-yc)))

### trnlldsp()
# Variables 'midpt', 'flux' in mandelecl, 'rprs', 'cosi', 'ars', 'p' in mandeltr, 't' in expramp()

params = np.array((midpt, rprs, cosi, ars, flux, p, c1, c2, c3, c4))

yp  = mp.trnlldsp(params, t, [])
yc = mc.trnlldsp(params, t, [])

if sum(abs(yp-yc)) < 1e-10:
    print("Trnlldsp:   PASS")
else:
    print("Trnlldsp:   FAIL", sum(abs(yp-yc)))


### vsll()

# p#  = # ??? <-- Not used in program?

visparams = np.array((t0, a, b, c, t1))

yp  = mp.vsll(visparams, [t, knots], [])
yc = mc.vsll(visparams, [t, knots], [])

if sum(abs(yp-yc)) < 1e-10:
    print("Vsll:       PASS")
else:
    print("Vsll:       FAIL", sum(abs(yp-yc)))


