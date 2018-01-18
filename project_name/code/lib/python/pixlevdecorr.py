import numpy as np

def pld_box(data, center, pldhw, skylev):

  """
    Perform the first step in pixel-level decorrelation, a way of handling intra-pixel sensitivity

    Parameters
    ----------
    data    : 4D numpy array
              (time, y, x, npos) from event.data
    center  : 1x2 list
              (y, x) from event.targpos; target position on data frame
    pldhw   : integer
              half-width of box around center
    skylev  : 2D numpy array
              (npos, background subtraction for each frame)

    Returns
    -------
    This function returns a 3D numpy array (time, y, x) in which each
    pixel is normalized by the total flux in its frame

    Revisions
    ---------
    Written by: Hannah Diamond-Lowe, University of Chicago. 02-08-2014
               hdiamondlowe@uchicago.edu

  """


  time = len(data)
  size = 2*pldhw + 1

  ctr_y = int(np.round(center[0]))
  ctr_x = int(np.round(center[1]))
  pldhw = int(pldhw)

  # Create a box of width and length pldhw around the target positionon the frame
  pldbox = data[:, ctr_y-pldhw:ctr_y+pldhw+1, ctr_x-pldhw:ctr_x+pldhw+1, 0]

  # Subtract the background
  pldbox -= skylev[0,...,np.newaxis,np.newaxis]

  # Create the array of frames, with each pixel normalized by its frame's total flux
  P_hat = pldbox/np.sum(pldbox,axis=(1,2))[:,np.newaxis,np.newaxis]
  # if using an older version of numpy: P_hat = pldbox/np.sum(np.sum(pldbox,axis=1), axis=1)[:,np.newaxis,np.newaxis

  return P_hat
