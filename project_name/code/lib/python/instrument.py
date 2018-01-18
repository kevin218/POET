# $Author: patricio $
# $Revision: 304 $
# $Date: 2010-07-13 11:36:20 -0400 (Tue, 13 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/instrument.py $
# $Id: instrument.py 304 2010-07-13 15:36:20Z patricio $

import numpy as np

# set Spitzer-specific stuff in info structure
class Instrument:

  def __init__(self, chan):

    wavel = np.array([3.6, 4.5, 5.8, 8.0, 16.0, 24.0]) * 1e-6

    # Spitzer central wavelengths, in meters
    self.spitzwavl = wavel[chan-1]

    # Instrument name and Data channel directory
    if chan < 5:
      self.name    = 'irac'
      self.prefix  = 'I'
      self.channel = '/ch%d'%chan
    elif chan == 5:
      self.name    = 'irs'
      self.prefix  = 'S'
      self.channel = '/ch0'
    else:
      self.name    = 'mips'
      self.prefix  = 'M'
      self.channel = '/ch1'

    # Spitzer standard filename components
    self.bcddir  = self.channel + '/bcd/'
    self.rawdir  = self.channel + '/raw/'
    self.caldir  = self.channel + '/cal/'

    self.rawsuf   = '_dce.fits'       # raw data from spacecraft + header

    if   chan <= 4:
      self.bcdsuf   = '_bcd.fits'     # bcd image (basic calibrated data)
      self.buncsuf  = '_bunc.fits'    # bcd uncertainties
      self.bdmsksuf  = '_bimsk.fits'  # bcd outlier mask 
      self.bdmsksuf2 = '_bdmsk.fits'  # bcd outlier mask 
      self.brmsksuf = '_brmsk.fits'   # bcd outlier mask 
    elif chan == 5:
      self.bcdsuf   = '_bcdb.fits'    # bcd image (basic calibrated data)
      self.buncsuf  = '_uncb.fits'    # bcd uncertainties
      self.bdmsksuf = '_b_msk.fits'   # bcd outlier mask
      self.bdmsksuf2 = '_xxxx.fits'   # inelegant solution
      self.brmsksuf = '_mskb.fits'    # bcd outlier mask 
    else: # chan == 6:
      self.bcdsuf   = '_bcd.fits'     # bcd image (basic calibrated data)
      self.buncsuf  = '_bunc.fits'    # bcd uncertainties
      self.bdmsksuf  = '_bbmsk.fits'  # bcd outlier mask 
      self.bdmsksuf2 = '_xxxx.fits'   # inelegant solution
      self.brmsksuf = '_brmsk.fits'   # bcd outlier mask 

    self.msaicsuf = '_msaic.fits'     # pbcd mosaic image
    self.msuncsuf = '_msunc.fits'     # pbcd mosaic uncertainties
    self.mscovsuf = '_mscov.fits'     # pbcd mosaic coverage (number of images)
    self.irsasuf  = '_irsa.tbl'       # list of 2MASS sources in the field
    self.pmasksuf = '_pmask.fits'     # pointing-refinement-corrected keywords

    # Critical mask flags
    self.pcrit   = np.long(65535)  # in pmask (permanent bad-pixel
                                   # mask, IRACDH2.0, T4.1)
    if chan < 6:
      # in dmask (per-frame bad-pixel mask, IRACDH2.0, T4.2) added bit
      # 4 (decimal 16) since uncerts are high and flux is low in top
      # row, which has this flag
      self.dcrit   = np.int_(32560)  
    else:
      self.dcrit   = np.long(65024)

    # in IDL:
    #  pcrit   = 'FFFF'x
    #  dcrit   = '7F20'x
    #  dcrit  += '0010'x  irac/irs
    #  dcrit   = 'FE00'x  mips
