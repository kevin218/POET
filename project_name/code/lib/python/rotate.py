
import numpy as np

#Rotate frame of reference about the specified axis by some angle theta
def xaxis(y,z,theta):
    yy = y*np.cos(theta) - z*np.sin(theta)
    zz = y*np.sin(theta) + z*np.cos(theta)
    return yy,zz

def yaxis(x,z,theta):
    xx =  x*np.cos(theta) + z*np.sin(theta)
    zz = -x*np.sin(theta) + z*np.cos(theta)
    return xx,zz

def zaxis(x,y,theta):
    xx = x*np.cos(theta) - y*np.sin(theta)
    yy = x*np.sin(theta) + y*np.cos(theta)
    return xx,yy
