"""
Limb Darkening Laws:
From Transiting Exoplanets, C. Haswell

Example:
--------
import limbdark as ld
import numpy as np
import matplotlib.pyplot as plt

mu = np.linspace(0.02, 1.0, 50)
lin  = ld.linear(mu, 0.215)
log  = ld.log   (mu, 0.14,  -0.12)
quad = ld.quad  (mu, 0.29,  -0.13)

plt.figure(0)
plt.clf()
plt.plot(mu, lin,  'r', label='linear')
plt.plot(mu, log,  'g', label='log')
plt.plot(mu, quad, 'b', label='quad')
plt.legend(loc='best')

mu0 = np.cos(80*np.pi/180)
print(ld.linear(mu0, 0.215),
      ld.log(mu0, 0.14,  -0.12),
      ld.quad(mu0, 0.29,  -0.13) )

"""

import numpy as np

def linear(mu, u):
  """
  This routine calculates I(mu)/I(mu=1) with a linear LD law.
 
  Parameters:
  -----------
  mu: 1D ndarray
      mu=cos(gamma), where gamma is the traveling angle of the photon
      in the interior of the sun with respect to the line of sight.
  
  """
  return 1 - u*(1-mu)


def log(mu, (u, v)):
  """
  This routine calculates I(mu)/I(mu=1) with a logarithmic LD law.
  """
  return 1 - u*(1-mu) - v*mu*np.log(mu)


def quad(mu, (u, v)):
  """
  This routine calculates I(mu)/I(mu=1) with a quadratic LD law.
  """
  return 1 - u*(1-mu) - v*(1-mu)**2.0


def cubic(mu, (u, v)):
  """
  This routine calculates I(mu)/I(mu=1) with a cubic LD law.
  """
  return 1 - u*(1-mu) - v*(1-mu)**3.0


def nl(mu, (c0, c1, c2, c3)):
  """
  This routine calculates I(mu)/I(mu=1) with a non-linear LD law.
  """
  return 1 - c0*(1-mu**0.5) - c1*(1-mu) - c2*(1-mu**1.5) - c3*(1-mu**2)


