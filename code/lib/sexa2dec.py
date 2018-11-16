# $Author: carthik $
# $Revision: 267 $
# $Date: 2010-06-08 22:33:22 -0400 (Tue, 08 Jun 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/sexa2dec.py $
# $Id: sexa2dec.py 267 2010-06-09 02:33:22Z carthik $

import numpy as np

def sexa2dec(sexa):
  """
    Convert sexagesimal (time format) to decimal

    Parameters
    ----------
    sexa : String or list of strings of the form 'Snnn:mm:ss.sssss...' 
           where S is the sign. S can be '+' or '-' or can be ommited.

    Returns
    -------
    Float or numpy array of floats, with the value(s) of the input in
    decimal scale.

    Examples
    --------
    >>> import sexa2dec as s2d
    >>> print(s2d.sexa2dec('- 2:20:33.334'))
    -2.34259277778

    >>> print(33.0/60.0 + 32.334/3600.0)
    0.558981666667

    >>> print(s2d.sexa2dec('-00:33:32.334'))
    -0.558981666667

    >>> print(s2d.sexa2dec('- 0:33:32.334'))
    -0.558981666667

    >>> print(s2d.sexa2dec('  0:33:32.334'))
    0.558981666667

    >>> print(s2d.sexa2dec('+ 0:33:32.334'))
    0.558981666667

    >>> print(s2d.sexa2dec('+00:33:32.334'))
    0.558981666667

    >>> print(s2d.sexa2dec(['- 0:33:32.334', '+4']))
    [-0.55898167  4.        ]

    Revisions
    ---------
    2010-07-14  Patricio  Initial version.         pcubillos@fulbrightmail.org
    2010-10-27  patricio  Added docstring.
  """

  # Make sure it is an array.
  dim = np.ndim(sexa)
  sexaarr = np.array([sexa]) if dim == 0   else np.array(sexa)

  ns = np.size(sexaarr)
  result = np.zeros(ns, np.double)

  for i in np.arange(ns):
    # Split the string
    sexait = sexaarr[i].strip().split(':')

    # Determine sign
    neg  = sexait[0].find('-')
    pos  = sexait[0].find('+')
    sign =   -1 if neg != -1   else 1

    if   neg != -1 or pos != -1:
      sexait[0] = sexait[0][1:]

    # Calculate the number 
    hr  = np.double(sexait[0])
    min = np.double(sexait[1]) if len(sexait) >= 2 else 0.0
    sec = np.double(sexait[2]) if len(sexait) == 3 else 0.0

    result[i] = sign * ( hr + ( min + sec / 60.0 ) / 60.0 )

  # Return a scalar if we received one
  return   result[0] if dim == 0   else result
