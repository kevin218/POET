from numarray import *

#
# Copyright (c) 1994-2005, Research Systems, Inc.  All rights reserved.
#       Unauthorized reproduction prohibited.
# Modifications by Joseph Harrington are in the public domain.
#+
# NAME:
#       FX_ROOT_JH
#
# PURPOSE:
#       This function computes real and complex roots (zeros) of
#       a univariate nonlinear function.  This version improves on
#       that in the IDL release by offering _EXTRA, STATUS, and a
#       sanity check on TOL.
#
# CATEGORY:
#       Nonlinear Equations/Root Finding
#
# CALLING SEQUENCE:
#       Result = FX_ROOT(X, Func)
#
# INPUTS:
#       X :      A 3-element initial guess vector of type real or complex.
#                Real initial guesses may result in real or complex roots.
#                Complex initial guesses will result in complex roots.
#
#       Func:    A scalar string specifying the name of a user-supplied IDL
#                function that defines the univariate nonlinear function.
#                This function must accept the vector argument X.
#
# KEYWORD PARAMETERS:
#       DOUBLE:  If set to a non-zero value, computations are done in
#                double precision arithmetic.
#
#       ITMAX:   Set this keyword to specify the maximum number of iterations
#                The default is 100.
#
#       STOP:    Set this keyword to specify the stopping criterion used to
#                judge the accuracy of a computed root, r(k).
#                STOP = 0 implements an absolute error criterion between two
#                successively-computed roots, |r(k) - r(k+1)|.
#                STOP = 1 implements a functional error criterion at the
#                current root, |Func(r(k))|. The default is 0.
#
#       TOL:     Set this keyword to specify the stopping error tolerance.
#                If the STOP keyword is set to 0, the algorithm stops when
#                |x(k) - x(k+1)| < TOL.
#                If the STOP keyword is set to 1, the algorithm stops when
#                |Func(x(k))| < TOL. The default is 1.0e-4.
#		 Tol is limited to machine precision.  If set below
#		 precision, it will be reset to precision IN THE
#		 CALLER.
#
#	STATUS:	 (returned) Set to 0 if the algorithm did not
#		 converge, 1 if it did IN THE CALLER.
#
#	_EXTRA:	 Structure containing parameters to pass to FUNC.
#
# EXAMPLE:
#       Define an IDL function named FUNC.
#         function FUNC, x
#           return, exp(sin(x)^2 + cos(x)^2 - 1) - 1
#         end
#
#       Define a real 3-element initial guess vector.
#         x = [0.0, -!pi/2, !pi]
#
#       Compute a root of the function using double-precision arithmetic.
#         root = FX_ROOT(x, 'FUNC', /double)
#
#       Check the accuracy of the computed root.
#         print, exp(sin(root)^2 + cos(root)^2 - 1) - 1
#
#       Define a complex 3-element initial guess vector.
#         x = [complex(-!pi/3, 0), complex(0, !pi), complex(0, -!pi/6)]
#
#       Compute a root of the function.
#         root = FX_ROOT(x, 'FUNC')
#
#       Check the accuracy of the computed root.
#         print, exp(sin(root)^2 + cos(root)^2 - 1) - 1
#
# PROCEDURE:
#       FX_ROOT implements an optimal Muller's method using complex
#       arithmetic only when necessary.
#
# SIDE EFFECTS:
#	Sets STATUS and may set TOL IN THE CALLER.
#
# REFERENCE:
#       Numerical Recipes, The Art of Scientific Computing (Second Edition)
#       Cambridge University Press
#       ISBN 0-521-43108-5
#
# MODIFICATION HISTORY:
#       Written by:  GGS, RSI, March 1994
#       Modified:    GGS, RSI, September 1994
#                    Added support for double-precision complex inputs.
#	2005-02-07 jh	Added _extra.
#	2005-02-12 jh	Added status, tol protection.  Fixed indentation.
#-

def fx_root_jh(xi, func, double=None, itmax=None, stop=None, tol=None, status=None, extra=None):

#on_error, 2 ;Return to caller if error occurs.

   e = extra
   
   status = 0
   
   x = xi + 0.0 #Create an internal floating-point variable, x.
   sx = size(x)
   if sx[1] != 3:   
      message('x must be a 3-element initial guess vector.')
   
   #Initialize keyword parameters.
   """
   if (double is not None) != 0:   
      if bitwise_or(sx[2] == 4, sx[2] == 5):   
         x = x + 0.0e0
      else:   
         x = dcomplex(x)
   """
   tn = size(x, tnam=True)
   if (itmax is not None) == 0:   
      itmax = 100
   if (stop is not None) == 0:   
      stop = 0
   if (tol is not None) == 0:   
      tol = 1.0e-4
   # protect against division by zero from too small a tol
   if bitwise_or(tn == 'DOUBLE', tn == 'DCOMPLEX'):   
      tol = maximum(tol, (machar(d=True)).eps)
   else:   
      tol = maximum(tol, (machar()).eps)
   
   #Initialize stopping criterion and iteration count.
   cond = 0
   it = 0
   
   #Begin to iteratively compute a root of the nonlinear function.
   while (it < itmax and cond != 1):
      q = (x[2] - x[1]) / (x[1] - x[0])
      pls = (1 + q)
      f = call_function(func, x, extra=e)
      a = q * f[2] - q * pls * f[1] + q ** 2 * f[0]
      b = (2 * q + 1) * f[2] - pls ** 2 * f[1] + q ** 2 * f[0]
      c = pls * f[2]
      disc = b ** 2 - 4 * a * c
      roc = size(disc)  #Real or complex discriminant?
      if bitwise_and(roc[1] != 6, roc[1] != 9):    #Proceed toward real root.
         if disc < 0:    #Switch to complex root.
            #Single-precision complex.
            if bitwise_and((double is not None) == 0, sx[2] != 9):   
               r0 = b + complex(0, sqrt(abs(disc)))
               r1 = b - complex(0, sqrt(abs(disc)))
            else:       #Double-precision complex.
               r0 = b + dcomplex(0, sqrt(abs(disc)))
               r1 = b - dcomplex(0, sqrt(abs(disc)))
            if abs(r0) > abs(r1):   
               div = r0
            else:   
               div = r1
         else:           # real root
            rr0 = b + sqrt(disc)
            rr1 = b - sqrt(disc)
            div = ((abs(rr0) >= abs(rr1)) and [rr0] or [rr1])[0]
      else:               #Proceed toward complex root.
         c0 = b + sqrt(disc)
         c1 = b - sqrt(disc)
         if abs(c0) > abs(c1):   
            div = c0
         else:   
            div = c1
      root = x[2] - (x[2] - x[1]) * (2 * c / div)
      #Absolute error tolerance.
      if bitwise_and(stop == 0, abs(root - x[2]) <= tol):   
         cond = 1
      else:   
         evalfunc = call_function(func, root, extra=e)
         #Functional error tolerance.
         if bitwise_and(stop != 0, abs(evalfunc) <= tol):   
            cond = 1
         else:   
            if evalfunc == 0:   
               cond = 1
      x = concatenate([x[1], x[2], root])
      it = it + 1
   if bitwise_and(it >= itmax, cond == 0):   
      print('Algorithm failed to converge within given parameters.')
   else:   
      status = 1
   
   return root

