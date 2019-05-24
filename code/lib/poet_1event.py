# $Author: patricio $
# $Revision: 302 $
# $Date: 2010-07-10 01:35:30 -0400 (Sat, 10 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_1event.py $
# $Id: poet_1event.py 302 2010-07-10 05:35:30Z patricio $

import os, sys, re, time

#FINDME Megan added this
#if '/Users/megan/auto_aor' in sys.path:
#  sys.path.remove('/Users/megan/auto_aor')
  
import numpy  as np
from astropy.io import fits as pf
import matplotlib.pyplot as plt
import sexa2dec      as s2d
import poet_dataread as pdr
import reader3       as rd
import tepclass      as tc
import instrument    as inst
import logedit       as le
import timer         as t
from   manageevent import *
from univ import Univ



class Event(Univ):

  def __init__(self, eventpcf):
    tini = time.time()
    # Open new log based on the file name
    logname = eventpcf[:-4] + "ini.log"
    global log
    log = le.Logedit(logname)
    self.logname = logname

    # initialize Univ
    Univ.__init__(self)
    pcf = rd.read_pcf(eventpcf)
    self.initpars(pcf)
    self.calc(pcf)
    self.read()
    self.check()
    self.save()

    # Print time elapsed and close log:
    cwd = os.getcwd() + "/"
    log.writelog("\nOutput files:")
    log.writelog("Data:")
    log.writelog(" " + cwd + self.eventname + "_ini.dat")
    log.writelog(" " + cwd + self.eventname + "_ini.h5")
    log.writelog("Log:")
    log.writelog(" " + cwd + logname)
    log.writelog("Figures:")
    log.writelog(" " + cwd + self.eventname + "-fig101.png")

    dt = t.hms_time(time.time()-tini)
    log.writeclose('\nEnd init and read. Time (h:m:s):  %s'%dt)


  def initpars(self, pcf):
    """
      add docstring.
    """

    log.writelog( 'MARK: ' + time.ctime() + ' : New Event Started ')

    # Read Planet Parameters From TepFile
    # ccampo 5/27/2011:
    # origtep is a tepfile with it's original units (not converted)
    self.origtep      = tc.tepfile(pcf.tepfile.value[0], conv=False)
    tep               = tc.tepfile(pcf.tepfile.value[0])
    self.tep          = tep
    self.ra           = tep.ra.val
    self.dec          = tep.dec.val
    self.rstar        = tep.rs.val
    self.rstarerr     = tep.rs.uncert
    self.metalstar    = tep.feh.val
    self.metalstarerr = tep.feh.uncert
    self.tstar        = tep.ts.val
    self.tstarerr     = tep.ts.uncert
    self.logg         = tep.loggstar.val
    self.loggerr      = tep.loggstar.uncert
    self.rplan        = tep.rp.val
    self.rplanerr     = tep.rp.uncert
    self.semimaj      = tep.a.val
    self.semimajerr   = tep.a.uncert
    self.incl         = tep.i.val
    self.inclerr      = tep.i.uncert
    self.ephtime      = tep.ttrans.val
    self.ephtimeerr   = tep.ttrans.uncert
    self.period       = tep.period.val      / 86400.0  # conv to days
    self.perioderr    = tep.period.uncert   / 86400.0  # ditto
    self.transdur     = tep.transdur.val
    self.transdurerr  = tep.transdur.uncert



    # original code - tep conversions
    """
    self.tep          = tc.tepfile(pcf.tepfile.value[0])
    self.ra           = s2d.sexa2dec(tep.ra.val)  /  12.0 * np.pi
    self.dec          = s2d.sexa2dec(tep.dec.val) / 180.0 * np.pi
    self.rstar        = tep.rs.val    * self.rsun
    self.rstarerr     = tep.rs.uncert * self.rsun
    self.metalstar    = tep.feh.val
    self.metalstarerr = tep.feh.uncert
    self.tstar        = tep.ts.val
    self.tstarerr     = tep.ts.uncert
    self.logg         = tep.loggstar.val
    self.loggerr      = tep.loggstar.uncert
    self.rplan        = tep.rp.val    * self.rjup
    self.rplanerr     = tep.rp.uncert * self.rjup
    self.semimaj      = tep.a.val    * self.au
    self.semimajerr   = tep.a.uncert * self.au
    self.incl         = tep.i.val    * np.pi / 180.0
    self.inclerr      = tep.i.uncert * np.pi / 180.0
    self.ephtime      = tep.ttrans.val
    self.ephtimeerr   = tep.ttrans.uncert
    self.period       = tep.period.val
    self.perioderr    = tep.period.uncert
    self.transdur     = tep.transdur.val    * 86400.0
    self.transdurerr  = tep.transdur.uncert * 86400.0
    """

    self.arat    = (self.rplan / self.rstar)**2
    self.araterr = 2*(self.arat)*np.sqrt( (self.rplanerr/self.rplan)**2 +
                                          (self.rstarerr/self.rstar)**2   )

    # position corrections:
    if pcf.ra.get(0) != None:
      self.ra  = s2d.sexa2dec(pcf.ra.value[0] ) /  12.0 * np.pi

    if pcf.dec.get(0) != None:
      self.dec = s2d.sexa2dec(pcf.dec.value[0]) / 180.0 * np.pi

    # Convert units to uJy/pix
    self.fluxunits    = pcf.fluxunits.get(0)
    #FINDME: are we using converted??
    #self.converted    = pcf.getvalue(False)  # are units uJy/pix ?

    # Initialize control file parameters
    self.planetname = pcf.planetname.getarr()
    if np.size(self.planetname) > 1:
      self.planetname = '-'.join(self.planetname)
    else:
      self.planetname = self.planetname[0]
    self.planet       = pcf.planet.get(0)
    self.ecltype      = pcf.ecltype.get(0)
    self.photchan     = pcf.photchan.get(0)
    if self.photchan < 5:
      self.instrument = 'irac'
    elif self.photchan == 5:
      self.instrument = 'irs'
    else:
      self.instrument = 'mips'

    # Instrument contains the instrument parameters
    # FINDME: inherit Instrument instead of adding as a parameter ?
    self.inst = inst.Instrument(self.photchan)

    self.visit     = pcf.visit.get(0)
    self.sscver    = pcf.sscver.get(0)

    # FINDME: fix obsdate
    # FINDME2: hey, we are not using lcfile!
    self.lcfile    = ( self.planetname + '_' + self.ecltype          +
                       '_OBSDATE_Spitzer_'   + self.inst.name + '_'  +
                       '%.1f_microns.fits'%(self.inst.spitzwavl*1e6) )

    # Directories
    self.topdir    = pcf.topdir.get(0)
    self.datadir   = pcf.datadir.get(0)
    self.dpref     = ( self.topdir + '/'  + self.datadir + '/' +
                       self.sscver + '/r' )

    # aors
    self.aorname   = pcf.aorname.value  # get aorname as string
    self.aortype   = pcf.aortype.getarr()

    # Number of aors per event
    self.naor      = np.size(self.aorname[np.where(self.aortype == 0)])

    # Number of position and noddings
    self.npos      = pcf.npos.get(0)
    self.nnod      = pcf.nnod.get(0)

    # Variables added for calibration AORs
    # self.calnmcyc  = pcf.getvalue('calnmcyc')

    # Ancil files
    self.hordir    = pcf.hordir.get(0)
    self.kuruczdir = pcf.kuruczdir.get(0)
    self.filtdir   = pcf.filtdir.get(0)
    self.psfdir    = pcf.psfdir.get(0)
    self.leapdir   = pcf.leapdir.get(0)

    self.pmaskfile  = np.zeros(self.naor, '|S150')
    for i in np.arange(self.naor):
      self.pmaskfile[i] = (self.dpref + self.aorname[i] + self.inst.caldir +
                           pcf.pmaskfile.get(0) )

    self.horvecfile = self.topdir + self.hordir    + pcf.horfile.get(0)
    self.kuruczfile = self.topdir + self.kuruczdir + pcf.kuruczfile.get(0)

    filt = pcf.filtfile.getarr()
    if self.photchan < 5:
      filt = re.sub('CHAN', str(self.photchan), filt[0])
    elif self.photchan == 5:
      filt = filt[1]
    else: # self.photchan == 5:
      filt = filt[2]
    self.filtfile   = self.topdir + self.filtdir + filt

    # Default PSF file:
    if self.photchan < 5  and pcf.psffile.get(0) == "default":
      self.psffile = (self.topdir + self.psfdir + 'IRAC PSF/' +
                      'IRAC.%i.PRF.5X.070312.fits'%self.photchan )
    # User specified PSF file:
    else:
      self.psffile =  pcf.psffile.get(0)

    # Bad pixels:
    # Chunk size
    self.szchunk      = pcf.szchunk.get(0)
    # Sigma rejection threshold
    self.sigma        = pcf.sigma.getarr()
    # User rejected pixels
    self.userrej      = pcf.userrej.getarr()
    if self.userrej[0] != None:
      self.userrej = self.userrej.reshape(np.size(self.userrej)/2, 2)
    else:
      self.userrej = None

    # set event directory
    self.eventdir = os.getcwd()


  def calc(self, pcf):
    """
      Add docstring.
    """
    # Instrument Channel
    self.spitzchan = ( self.photchan if self.photchan <= 4 else
                                  (0 if self.photchan == 5 else 1) )

    # Name of the event
    self.eventname = ( self.planet           + self.ecltype +
                       np.str(self.photchan) + np.str(self.visit) )

    # Added to check whether ancillary data files exist
    self.ispmask  = np.zeros(self.naor, bool)
    for i in np.arange(self.naor):
      self.ispmask[i] = os.path.isfile(self.pmaskfile[i])
    self.ishorvec     = os.path.isfile(self.horvecfile)
    self.iskurucz     = os.path.isfile(self.kuruczfile)
    self.isfilt       = os.path.isfile(self.filtfile  )
    self.ispsf        = os.path.isfile(self.psffile   )
    self.isleapdir    = os.path.isdir (self.topdir + self.leapdir)

    # Calibration aors
    self.havecalaor     = 0 if np.all(self.aortype==0) else 1
    self.calnaor        = np.size(self.aorname[np.where(self.aortype==1)])
    self.precalaorname  = self.aorname[np.where(self.aortype==1)]
    self.postcalaorname = self.aorname[np.where(self.aortype==2)]

    # Array containing the number of expositions per AOR
    self.nexpid  = np.zeros(self.naor, np.long)


    # compile patterns: lines ending with each suffix
    bcdpattern    = re.compile("(.+" + self.inst.bcdsuf    + ")\n")
    bdmskpattern  = re.compile("(.+" + self.inst.bdmsksuf  + ")\n")
    bdmsk2pattern = re.compile("(.+" + self.inst.bdmsksuf2 + ")\n")

    # make list of files in each AOR
    self.bcdfiles = []
    for aornum in np.arange(self.naor):
      dir = self.dpref + self.aorname[aornum] + self.inst.bcddir
      frameslist = os.listdir(dir)
      framesstring = '\n'.join(frameslist) + '\n'

      # find the data files
      bcdfiles = bcdpattern.findall(framesstring)
      # and sort them
      self.bcdfiles.append(sorted(bcdfiles))

      # find bdmask suffix:
      if bdmskpattern.findall(framesstring) != []:
        self.masksuf = self.inst.bdmsksuf
      elif bdmsk2pattern.findall(framesstring) != []:
        self.masksuf = self.inst.bdmsksuf2
      else:
        log.writelog("No mask files found")
        self.masksuf = None
        #return

      # get first index of exposition ID, number of expID, and ndcenum
      #                    expid      dcenum     pipev
      first = re.search("_([0-9]{4})_([0-9]{4})_([0-9])",self.bcdfiles[-1][0])
      last  = re.search("_([0-9]{4})_([0-9]{4})_([0-9])",self.bcdfiles[-1][-1])

      self.expadj         = int(first.group(1))
      self.nexpid[aornum] = int(last.group(1)) + 1 - self.expadj
      self.ndcenum        = int(last.group(2)) + 1
      self.pipev          = int(last.group(3))

    # FINDME: delete this
    # self.calnexpid      = self.nexpid[np.where(self.aortype==1)]

    # pick a random image, not the first
    file = bcdfiles[-2]
    data, head = pf.getdata(dir + file, header=True)

    # data size
    shape = data.shape
    if data.ndim >= 3:
      self.nz = shape[0]
      self.ny = shape[1]
      self.nx = shape[2]
    else:
      self.nz = 1
      self.ny = shape[0]
      self.nx = shape[1]

    # number of lines in the header
    self.nh = len(head)

    # number of small, medium and big cycles
    if self.instrument == 'irac':
      self.nbcyc = 1
      self.nmcyc = np.sum(self.nexpid)
      self.nscyc = self.ndcenum if self.nz == 1 else self.nz
    elif self.instrument == 'irs':
      self.nbcyc = np.sum(self.nexpid)/self.nnod
      self.nmcyc = self.ndcenum
      self.nscyc = 1
    else: # self.instrument == 'mips'
      self.nbcyc = np.sum(self.nexpid)/self.nnod
      self.nmcyc = (self.ndcenum - 1)/7
      self.nscyc = 7

    # max. number of images per position
    if self.instrument == 'mips':
      self.maxnimpos  = self.nbcyc * (self.nmcyc + 1)
    else:
      self.maxnimpos  = ( np.sum(self.nexpid[np.where(self.aortype==0)]) *
                          self.ndcenum * self.nz / self.nnod )
      #check ineger number
    # calmaxnimpos
    self.calmaxnimpos = ( np.sum(self.nexpid[np.where(self.aortype==1)]) *
                          self.ndcenum )

    try:
      self.framtime = head['FRAMTIME']   # interval between exposure starts
    except:
      self.framtime = 0.0
    try:
      self.exptime  = head['EXPTIME']    # effective exposure time
    except:
      self.exptime  = None
    try:
      self.gain     = head['GAIN']       # e/DN conversion
    except:
      self.gain     = None

    self.bunit      = head['BUNIT']      # Units of image data
    self.fluxconv   = head['FLUXCONV']   # Flux Conv factor (MJy/Str per DN/sec)
    self.posscl       = np.zeros((2,self.npos))
    self.posscl[0, :] = np.abs(head['PXSCAL2']) # ["/pix] axis 2 @ CRPIX1,CRPIX2
    self.posscl[1, :] = np.abs(head['PXSCAL1']) # ["/pix] axis 1 @ CRPIX1,CRPIX2

    if self.spitzchan != head['CHNLNUM']:  # Spitzer photometry channel
      log.writelog( 'poet_calc: photometry channel unexpected')

    # Frequency calculated from wavelength
    self.freq = self.c / self.inst.spitzwavl


  def read(self):
    """
      add docstring.
      Read Data
    """
    pdr.poet_dataread(self, log=log)

  def check(self):
    """
      add docstring.
    """

    # Source estimated position
    self.srcest = np.zeros((2,self.npos))
    for p in np.arange(self.npos):
      self.srcest[0,p] = np.round(np.average(self.fp.heady[p,0:self.nimpos[p]]))
      self.srcest[1,p] = np.round(np.average(self.fp.headx[p,0:self.nimpos[p]]))

    # Plot a reference image
    image = np.zeros((self.ny, self.nx))
    for pos in np.arange(self.npos):
      image += self.data[0,:,:,pos]

    image[np.where(np.isfinite(image) != 1)] = 0
    plt.figure(101, (10,9))
    plt.clf()
    plt.imshow(image, interpolation='nearest', origin='ll', cmap=plt.cm.gray)
    plt.plot(self.srcest[1,:], self.srcest[0,:],'r+')
    plt.xlim(0,self.nx-0.5)
    plt.ylim(0,self.ny-0.5)
    plt.title(self.eventname + ' reference image')
    plt.savefig(self.eventname + "-fig101.png")

    # Throw a warning if the source estimate position lies outside of
    # the image.
    wrning = False
    if (np.any(self.srcest[1,:] < 0) or np.any(self.srcest[1,:] > self.nx)
     or np.any(self.srcest[0,:] < 0) or np.any(self.srcest[0,:] > self.ny)):
      wrning = True

    # Write to log
    log.writelog('\n%s: event %s'%(self.planetname, self.eventname))
    log.writelog('nexpid  = ' + np.str(self.nexpid))
    log.writelog('ndcenum = %d'%self.ndcenum)
    #log.writelog('you should see %d positions in the output image '%self.npos +
    #            '(red crosses)')
    log.writelog("Target guess position:\n" + str(self.srcest[0:2,:]))
    if wrning:
       log.writelog('Source position estimate out of bounds!')
    log.writelog('nimpos  = ' + np.str(self.nimpos))
    log.writelog('Read %d frames\n'%np.sum(self.nimpos))

    # Report files not found:
    print("Ancil Files:")
    if not self.ispmask[0]:
      log.writelog('Pmask:       File not found!')
    else:
      log.writelog("Pmask:       " + str(self.pmaskfile[0].decode('utf-8')))
    if not self.ishorvec:
      log.writelog('Horizon:     File not found!')
    else:
      log.writelog("Horizon:     " + self.horvecfile)
    if not self.iskurucz:
      log.writelog('Kurucz:      File not found!')
    else:
      log.writelog("Kurucz:      " + self.kuruczfile)
    if not self.isfilt:
      log.writelog('Filter:      File not found!')
    else:
      log.writelog("Filter:      " + self.filtfile)
    if not self.ispsf:
      log.writelog('PSF:         Not supplied.')
    else:
      log.writelog("PSF:         " + self.psffile)
    if not self.isleapdir:
      log.writelog('Leap Second: Directory not found!')
    else:
      log.writelog("Leap Second: " + self.topdir + self.leapdir)

    if self.exptime  == None:
      log.writelog("Exposure time undefined.")
    if self.gain     == None:
      log.writelog("Gain undefined.")

  def save(self):
    # what to load in p2_badpix
    self.loadnext = ['data', 'uncd', 'bdmskd']

    if self.instrument == 'mips':
      self.loadnext.append('brmskd')
      saveevent(self, self.eventname + "_ini",
                save=['data', 'uncd', 'head', 'bdmskd', 'brmskd'])
    else:
      saveevent(self, self.eventname + "_ini", delete=['brmskd'],
                save=['data', 'uncd', 'head', 'bdmskd'])
