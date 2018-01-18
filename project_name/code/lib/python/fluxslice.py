
condition = ((fit.x < -.39) & (fit.x > -.40))
fluxucslice = fit.fluxuc.compress(condition)
yslice = fit.y.compress(condition)
ipparams = pars[1][0]
a    = ipparams[0]
b    = ipparams[1]
c    = ipparams[2]
d    = ipparams[3]
e    = ipparams[4]
f    = ipparams[5]
g    = ipparams[6]
h    = ipparams[7]
i    = ipparams[8]
j    = ipparams[9]
y    = fit.y
x    = -0.4
ipc  = a*y**3 + b*x**3 + c*y**2*x + d*y*x**2 + e*y**2 + f*x**2 + g*y*x + h*y + i*x + j

plt.figure(2)
plt.plot(yslice,fluxucslice,'.',ms=1)
plt.plot(y, ipc, 'k-')

condition = ((fit.y < -.04) & (fit.y > -.05))
fluxucslice = fit.fluxuc.compress(condition)
xslice = fit.x.compress(condition)
ipparams = pars[1][0]
a    = ipparams[0]
b    = ipparams[1]
c    = ipparams[2]
d    = ipparams[3]
e    = ipparams[4]
f    = ipparams[5]
g    = ipparams[6]
h    = ipparams[7]
i    = ipparams[8]
j    = ipparams[9]
y    = -0.04
x    = fit.x
ipc  = a*y**3 + b*x**3 + c*y**2*x + d*y*x**2 + e*y**2 + f*x**2 + g*y*x + h*y + i*x + j

plt.figure(4)
plt.plot(xslice,fluxucslice,'.',ms=1)
plt.plot(x, ipc, 'k.')

condition = ((y < -.04) & (y > -.05))
fluxucslice = fluxuc.compress(condition)
xslice = x.compress(condition)
ipparams = pars[1][0]
a    = ipparams[0]
b    = ipparams[1]
c    = ipparams[2]
d    = ipparams[3]
e    = ipparams[4]
y2    = 0.384
ipq = a*y2**2 + b*x**2 + c*y2 + d*x + e
plt.figure(2)
plt.plot(xslice,fluxucslice,'.',ms=1)
plt.plot(x, ipq, 'r.')
