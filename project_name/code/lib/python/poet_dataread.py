# $Author: patricio $
# $Revision: 302 $
# $Date: 2010-07-10 01:35:30 -0400 (Sat, 10 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_dataread.py $
# $Id: poet_dataread.py 302 2010-07-10 05:35:30Z patricio $


import numpy   as np
#import pyfits  as pf
from astropy.io import fits as pf
# import pywcs   as pw
from astropy import wcs as pw
import time, re
import univ

class FrameParameters:
  """
  class holder of the frame parameters.
  """
  def __init__(self):
    pass

def poet_dataread(event, type=0, log=None):
  """
    This function reads a set of IRAC AORS, (or IRAC Subarray AORS),
       sorting by dither position, if any.

    Parameters:
    ----------
    event : An event object.
    type  : integer
            Specifies the type of data to read.
            0 = data, 1 = precalibration data, 2 = postcalibration data.
    log   : A logedit object that keeps the log.

    Outputs:
    -------
    data  : [maxnimpos, nx, ny, npos] float array containing the data
            frames, sorted into dither positions.
    head  : header data
    uncd  : [maxnimpos, nx, ny, npos] float array, uncertainties
    bdmskd: [maxnimpos, nx, ny, npos] int array, per-pixel data flag
    nimpos: array like
            array containing the number of frames in each position.
    fp:     FrameParameters object containing [npos, maxnimpos] double arrays
            of per-frame parameters.

    Example:
    -------

    Modification History:
    --------------------
    Written by:	Joseph Harrington, Cornell.
    2005-09-16  jh@oobleck.astro.cornell.edu
    2005-10-26  jh        Fixed frame times.
    2005-10-27	jh        Moved calculation of some constants out of the
		          routine.  Filled in header.  Corrected some
		          array datatypes.
    2005-11-25	jh        Converted to using FP array.
    2006-01-04	jh        Header tweak.
    2006-03-20  jh        Added zodi, ism, cib header values.
    2007-03-07  khorning  Adapted program to use for non-subarray data
    2007-07-15  jh        Made nimpos be a long, not integer, array.
    2010-08-24  patricio  Converted to python.
    2010-10-27  patricio  Comments added.

  """
  # General variables
  dpref     = event.dpref          # data directory preffix
  expadj    = event.expadj         # id number of fisrt image
  ndcenum   = event.ndcenum        # number of dcenum
  npos      = event.npos           # number of position
  nnod      = event.nnod           # number of nodding positions
  #fpref     = event.fpref          # file names prefix
  pipev     = event.pipev          # spitzer pipeline version
  bcddir    = event.inst.bcddir    # directory containing bcd files
  bcdsuf    = event.inst.bcdsuf    # bcd file suffix
  buncsuf   = event.inst.buncsuf   # uncertainities file preffix
  #bdmsksuf  = event.inst.bdmsksuf  # badpixelmask file suffix
  brmsksuf  = event.inst.brmsksuf  # badpixelmask file suffix
  masksuf   = event.masksuf        # badpixelmask file suffix
  nx        = event.nx             #
  ny        = event.ny             #
  nz        = event.nz             # number of subarrays in datafile
  nh        = event.nh             #
  framtime  = event.framtime       #
  bcdlist   = event.bcdfiles       # lists of files

  # AORs/cal AORs variables
  aorname   = event.aorname[np.where(event.aortype==type)]
  if type !=0:
    naor      = event.calnaor
    nexpid    = event.calnexpid
    maxnimpos = event.calmaxnimpos
    nmcyc     = event.calnmcyc
  else:
    naor      = event.naor
    nexpid    = event.nexpid
    maxnimpos = event.maxnimpos
    nmcyc     = event.nmcyc
    nscyc     = event.nscyc

  # Allocate space for returned arrays
  headerdtype = 'a'+str(nh*81)
  #head   = np.zeros( (maxnimpos / nz, npos), dtype=headerdtype)
  print(maxnimpos/nz)
  maxnimpos = int(maxnimpos)
  print(maxnimpos/nz)
  print(npos)
  print(nz)

  #nz = int(nz)
  #npos = int(npos)
  #print(headerdtype)
  tmp = int(maxnimpos / nz)
  head   = np.zeros( (tmp, npos), dtype=headerdtype)
  data   = np.zeros( (maxnimpos, ny, nx,  npos), dtype=float)
  uncd   = np.zeros( (maxnimpos, ny, nx,  npos), dtype=float)
  bdmskd = np.zeros( (maxnimpos, ny, nx,  npos), dtype=int)
  brmskd = np.zeros( (maxnimpos, ny, nx,  npos), dtype=int)


  # Allocate space for the frame paramaters
  fp = FrameParameters()
  fpsize = np.zeros((npos, maxnimpos))
  fp.frmobs   = np.copy(fpsize)  # sequential frame number
  fp.pos      = np.copy(fpsize)  # position number
  fp.aor      = np.copy(fpsize)  # sequential AOR number
  fp.expid    = np.copy(fpsize)  # EXPosure ID
  fp.dce      = np.copy(fpsize)  # Data Collection Event
  fp.subarn   = np.copy(fpsize)  # subarray frame number
  fp.time     = np.copy(fpsize)  # frame mid-time, seconds J2000.0
  fp.bmjd     = np.copy(fpsize)  # frame mid-time, BJD-2400000.5
  fp.zodi     = np.copy(fpsize)  # zodiacal light estimate, see header comment
  fp.ism      = np.copy(fpsize)  # interstellar medium estimate,see head comment
  fp.cib      = np.copy(fpsize)  # cosmic infrared background,see header comment
  fp.afpat2b  = np.copy(fpsize)  # temperatures, K, see header comment
  fp.afpat2e  = np.copy(fpsize)
  fp.ashtempe = np.copy(fpsize)
  fp.atctempe = np.copy(fpsize)
  fp.acetempe = np.copy(fpsize)
  fp.apdtempe = np.copy(fpsize)
  fp.acatmp1e = np.copy(fpsize)
  fp.acatmp2e = np.copy(fpsize)
  fp.acatmp3e = np.copy(fpsize)
  fp.acatmp4e = np.copy(fpsize)
  fp.acatmp5e = np.copy(fpsize)
  fp.acatmp6e = np.copy(fpsize)
  fp.acatmp7e = np.copy(fpsize)
  fp.acatmp8e = np.copy(fpsize)
  # mips frame parameters
  fp.cmd_t_24 = np.copy(fpsize)
  fp.ad24tmpa = np.copy(fpsize)
  fp.ad24tmpb = np.copy(fpsize)
  fp.acsmmtmp = np.copy(fpsize)
  fp.aceboxtm = np.copy(fpsize)
  fp.pxscl2   = np.copy(fpsize)
  fp.pxscl1   = np.copy(fpsize)

  fp.heady    = np.copy(fpsize)
  fp.headx    = np.copy(fpsize)
  fp.filename = np.zeros((npos, maxnimpos), dtype='S150')

  nimpos = np.zeros(npos, np.long)

  # conveniences
  salist  = np.arange(nz)
  sadind  = np.arange(nz, dtype=np.double)

  # position of the star
  sky = [[event.ra*180/np.pi, event.dec*180/np.pi]]

  # dictionary to get position in MIPS
  mirind = {1929.:0, 2149.5:1, 1907.5:2, 2128.:3,
            1886.:4, 2106.5:5, 1864.5:6}

  # Write to log first line
  if log != None:
    log.writelog('  aor  expid  dcenum   pos')
  else:
    print('  aor  expid  dcenum   pos')

  # pattern to find     expid      dcenum
  pattern = re.compile("_([0-9]{4})_([0-9]{4})_")


  # Obtain data
  for aor in np.arange(naor):
    direc   = dpref + aorname[aor] + bcddir
    bcd   = bcdlist[aor]

    for i in np.arange(len(bcd)):
      # Read data
      try:
        dataf   = pf.getdata(direc + bcd[i])
        bcdhead   = pf.getheader(direc + bcd[i])
      except: # If a file doesn't exist, skip to next file.
        log.writelog(direc + bcd[i] + " File not found!")
        continue

      try: # Read uncertainity and mask files
        # Replace suffix in bcd file to get the corresponding file.
          uncfile = re.sub(bcdsuf, buncsuf, direc + bcd[i])
          uncf    = pf.getdata(uncfile)
          mskfile = re.sub(bcdsuf, masksuf, direc + bcd[i])
          bdmskf  = pf.getdata(mskfile)
      except:
        bdmskf    = np.ones((nz, ny, nx), np.long)

      try: # Mips
        brmskfile = re.sub(bcdsuf, brmsksuf, direc + bcd[i])
        brmskf    = pf.getdata(brmskfile)
      except:
        brmskf    = -np.ones((nz, ny, nx), np.long)

      # Obtain expid and dcenum
      index = pattern.search(bcd[i])
      expid  = int(index.group(1))
      dcenum = int(index.group(2))

      # Do I really need this?
      if np.size(bdmskf) == 1:
        bdmskf = -np.ones((nz, ny, nx), np.long)
      if np.size(brmskf) == 1:
        brmskf = -np.ones((nz, ny, nx), np.long)

      # Find dither position
      try:
        pos = bcdhead['DITHPOS'] - 1
      except:
        pos = 0  # No dither position in stare data
      if event.inst.name == 'irs':
        pos = expid % npos
      elif event.inst.name == 'mips':
        nod = expid % nnod
        pos = nod * nscyc + mirind[bcdhead['CSM_PRED']]

      be = nimpos[pos]      # begining
      en = nimpos[pos] + nz # end

      # Store data
      data  [be:en, :, :, pos] = dataf.reshape( (nz,ny,nx))
      uncd  [be:en, :, :, pos] = uncf.reshape(  (nz,ny,nx))
      bdmskd[be:en, :, :, pos] = bdmskf.reshape((nz,ny,nx))
      brmskd[be:en, :, :, pos] = brmskf.reshape((nz,ny,nx))
      # All the single numbers per frame that we care about
      fp.frmobs[pos, be:en] = np.sum(nimpos) + salist
      fp.pos   [pos, be:en] = pos
      fp.aor   [pos, be:en] = aor
      fp.expid [pos, be:en] = expid
      fp.dce   [pos, be:en] = dcenum
      fp.subarn[pos, be:en] = salist
      fp.bmjd  [pos, be:en] = bcdhead['BMJD_OBS'] + framtime*(sadind+0.5)/86400.
      # ccampo 2011/3/18: changed to UTC from SCLK to avoid timing inconsistencies
      try:
        fp.time[pos, be:en] = bcdhead['UTCS_OBS'] + framtime*(sadind+0.5)
      except:
        pass
      try:
        fp.zodi    [pos, be:en] = bcdhead['ZODY_EST']
        fp.ism     [pos, be:en] = bcdhead['ISM_EST']
        fp.cib     [pos, be:en] = bcdhead['CIB_EST']
        fp.afpat2b [pos, be:en] = bcdhead['AFPAT2B']
        fp.afpat2e [pos, be:en] = bcdhead['AFPAT2E']
        fp.ashtempe[pos, be:en] = bcdhead['ASHTEMPE'] + 273.0
        fp.atctempe[pos, be:en] = bcdhead['ATCTEMPE'] + 273.0
        fp.acetempe[pos, be:en] = bcdhead['ACETEMPE'] + 273.0
        fp.apdtempe[pos, be:en] = bcdhead['APDTEMPE'] + 273.0
        fp.acatmp1e[pos, be:en] = bcdhead['ACATMP1E']
        fp.acatmp2e[pos, be:en] = bcdhead['ACATMP2E']
        fp.acatmp3e[pos, be:en] = bcdhead['ACATMP3E']
        fp.acatmp4e[pos, be:en] = bcdhead['ACATMP4E']
        fp.acatmp5e[pos, be:en] = bcdhead['ACATMP5E']
        fp.acatmp6e[pos, be:en] = bcdhead['ACATMP6E']
        fp.acatmp7e[pos, be:en] = bcdhead['ACATMP7E']
        fp.acatmp8e[pos, be:en] = bcdhead['ACATMP8E']
      except:
        pass

      try:
        fp.pxscl2[pos, be:en]   = np.abs(bcdhead['PXSCAL2'])
        fp.pxscl1[pos, be:en]   = np.abs(bcdhead['PXSCAL1'])
        fp.acatmp5e[pos, be:en] = bcdhead['CMD_T_24']
        fp.acatmp6e[pos, be:en] = bcdhead['AD24TMPA']
        fp.acatmp6e[pos, be:en] = bcdhead['AD24TMPB']
        fp.acatmp5e[pos, be:en] = bcdhead['ACSMMTMP']
        fp.acatmp6e[pos, be:en] = bcdhead['ACEBOXTM'] + 273.0
      except:
        pass

      # Store filename
      fp.filename[pos, be:en] = direc + bcd[i]

      # Store header
      head[np.int(nimpos[pos] / nz), pos] = np.str(bcdhead)

      # Header position of the star:
      bcdhead["NAXIS"] = 2

      wcs = pw.WCS(bcdhead, naxis=2)
      pix = wcs.wcs_world2pix(sky,0)
      fp.headx[pos, be:en] = pix[0][0]
      fp.heady[pos, be:en] = pix[0][1]

      # Print to log and screen:
      if log != None:
        log.writelog('%4d'%aor + '%7d'%expid + '%7d'%dcenum + '%7d'%pos)
      else:
        print('%4d'%aor + '%7d'%expid + '%7d'%dcenum + '%7d'%pos)

      nimpos[pos] += nz

  # frame tags in fp

  # where there exist data
  fp.exist = np.zeros((npos, maxnimpos), np.long)
  for pos in np.arange(npos):
    fp.exist[pos, 0:nimpos[pos]] = 1

  fp.im = np.copy(fpsize)  # Frame within position
  for pos in np.arange(npos):
    fp.im[pos, 0:nimpos[pos]] = np.arange(nimpos[pos], dtype=np.double)

  if event.inst.name != 'mips':
    fp.cycpos = np.trunc(fp.frmobs / (npos * nmcyc * nz)) # Cycle number
    fp.visobs = np.trunc(fp.frmobs / (nmcyc * nz))# Visit number within obs. set
    fp.frmvis  = fp.im % (nmcyc * nz)             # Frame within visit

  else:
    fp.cycpos = np.trunc(fp.frmobs / (2*ndcenum)) # Cycle number
    fp.visobs = np.trunc(fp.frmobs / ndcenum)     # Visit number within obs. set
    fp.frmvis = np.trunc(fp.frmobs % ndcenum)     # Frame within visit

    # Image scale:
    for pos in np.arange(npos):
      last = nimpos[pos]
      if np.all(fp.pxscl1[pos, 0:last] == fp.pxscl1[pos, 0]):
        event.posscl[1, pos] = np.abs(fp.pxscl1[pos, 0])
      if np.all(fp.pxscl2[pos, 0:last] == fp.pxscl2[pos, 0]):
        event.posscl[0, pos] = np.abs(fp.pxscl2[pos, 0])

  # Update event
  event.data   = data
  event.uncd   = uncd
  event.bdmskd = bdmskd
  event.brmskd = brmskd
  event.head   = head
  event.fp     = fp
  event.nimpos = nimpos

  return
