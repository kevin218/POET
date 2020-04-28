# $Author: patricio $
# $Revision: 291 $
# $Date: 2010-07-01 13:31:46 -0400 (Thu, 01 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_3center.py $
# $Id: poet_3center.py 291 2010-07-01 17:31:46Z patricio $

import numpy  as np
from astropy.io import fits as pf
import sys, time, os, shutil, copy
import reader3      as rd
import logedit      as le
import manageevent  as me
import centerdriver as cd
import timer        as t
from multiprocessing import Process, Array

"""
POET_P3CENTERING WORKFLOW:
--------------------------
This beautiful piece of code consists of three sections:
- run_centering
- centering
- do_center

"""


def centering(event, pcf, centerdir):
  tini = time.time()

  # Create centering log
  logname = event.logname
  log = le.Logedit(centerdir + "/" + logname, logname)
  log.writelog("\nStart " + centerdir + " centering: " + time.ctime())

  os.chdir(centerdir)

  # copy center.pcf in centerdir
  pcf.make_file("center.pcf")

  # Parse the attributes from the control file to the event:
  attrib = vars(pcf)
  keys = attrib.keys()
  for key in keys:
    setattr(event, key, attrib.get(key).get())

  # Check least asym parameters work:
  if event.method in ['lac', 'lag']:
    if event.ctrim < (event.cradius + event.csize) and event.ctrim is not 0:
      event.ctrim = event.cradius + event.csize + 1
      log.writelog('Trim radius is too small, changed to: %i'%event.ctrim)
    if event.psfctrim < (event.psfcrad + event.psfcsize) and event.psfctrim is not 0:
      event.psfctrim = event.psfcrad + event.psfcsize + 1
      log.writelog('PSF Trim radius is too small, changed to: %i'
                   %event.psfctrim)

  # Centering bad pixel mask:
  centermask = np.ones((event.ny, event.nx))
  if event.ymask is not None:
    ymask = np.asarray(event.ymask, int)
    xmask = np.asarray(event.xmask, int)
    for i in np.arange(len(ymask)):
      centermask[ymask[i], xmask[i]] = 0

  # PSF:
  # Re-evaluate if a PSF has been redefined:
  if event.newpsf is not None:
    event.ispsf = os.path.isfile(event.newpsf)
    if event.ispsf:
      event.psffile = event.newpsf
      log.writelog('The PSF file has been redefined!')
      log.writelog("PSF:     " + event.psffile)

  # PSF Centering:
  if event.ispsf:
    event.psfim = pf.getdata(event.psffile)
    # Guess of the center of the PSF (center of psfim)
    psfctrguess = np.asarray(np.shape(event.psfim))/2
    # Do not find center of PSF:
    if event.nopsfctr:
      event.psfctr = psfctrguess
    # Find center of PSF:
    else:
      '''
      if event.method == "bpf" or event.method == "ipf":
        method = "fgc"
      else:
        method = event.method
      event.psfctr, extra = cd.centerdriver(method, event.psfim, psfctrguess,
                                 event.psfctrim, event.psfcrad, event.psfcsize)
      '''
      # Always use 'fgc' on PSF, for testing
      event.psfctr, extra = cd.centerdriver("fgc", event.psfim, psfctrguess,
                                 event.psfctrim, event.psfcrad, event.psfcsize) #FINDME
    log.writelog('PSF center found.')
    print(event.psfctr) #FINDME
  else:
    event.psfim  = None
    event.psfctr = None
    log.writelog('No PSF supplied.')

  # Find center of the mean Image:
  event.targpos = np.zeros((2, event.npos))
  for pos in np.arange(event.npos):
    meanim = event.meanim[:,:,pos]
    guess  = event.srcest[:, pos]
    targpos, extra = cd.centerdriver(event.method, meanim,
                                     guess, event.ctrim,
                                     event.cradius, event.csize,
                                     fitbg=event.fitbg, psf=event.psfim,
                                     psfctr=event.psfctr, expand=event.expand)
    event.targpos[:,pos] = targpos
  log.writelog("Center position(s) of the mean Image(s):\n" +
               str(np.transpose(event.targpos)))
  
  # Inclusion ::::::::
  # Multy Process set up:
  # Shared memory arrays allow only 1D Arrays :(
  event.maxnimpos = int(event.maxnimpos)
  event.npos = int(event.npos)
  x       = Array("d", np.zeros(event.npos * event.maxnimpos))
  y       = Array("d", np.zeros(event.npos * event.maxnimpos))
  sx      = Array("d", np.zeros(event.npos * event.maxnimpos))
  sy      = Array("d", np.zeros(event.npos * event.maxnimpos))
  flux    = Array("d", np.zeros(event.npos * event.maxnimpos))
  sky     = Array("d", np.zeros(event.npos * event.maxnimpos))
  goodfit = Array("d", np.zeros(event.npos * event.maxnimpos))

  # Size of chunk of data each core will process:
  chunksize = event.maxnimpos/event.ccores + 1
  print("Number of cores: " + str(event.ccores))

  # Start Muti Procecess: ::::::::::::::::::::::::::::::::::::::
  processes = []
  for nc in np.arange(event.ccores):
    start =  nc    * chunksize # Starting index to process
    end   = (nc+1) * chunksize # Ending   index to process
    proc = Process(target=do_center, args=(start, end, event, centermask, log,
                                      x, y, sx, sy, flux, sky, goodfit))
    processes.append(proc)
    proc.start()
  # Make sure all processes finish their work:
  for nc in np.arange(event.ccores):
    processes[nc].join()

  # Put the results in the event. I need to reshape them:
  event.fp.x        = np.asarray(x      ).reshape(event.npos,event.maxnimpos)
  event.fp.y        = np.asarray(y      ).reshape(event.npos,event.maxnimpos)
  event.fp.sx       = np.asarray(sx     ).reshape(event.npos,event.maxnimpos)
  event.fp.sy       = np.asarray(sy     ).reshape(event.npos,event.maxnimpos)
  # If PSF fit:
  if event.method in ["ipf", "bpf"]:
    event.fp.flux    = np.asarray(flux   ).reshape(event.npos,event.maxnimpos)
    event.fp.psfsky  = np.asarray(sky    ).reshape(event.npos,event.maxnimpos)
    event.fp.goodfit = np.asarray(goodfit).reshape(event.npos,event.maxnimpos)
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # Pixel R position:
  event.fp.r = np.sqrt((event.fp.x % 1.0 - 0.5)**2.0 +
                       (event.fp.y % 1.0 - 0.5)**2.0 )

  log.writelog("End frames centering.")

  # Save
  print("\nSaving")
  if event.denoised:
    me.saveevent(event, event.eventname + "_ctr", save=['dendata', 'data',
                                                        'uncd', 'mask'])
  else:
    me.saveevent(event, event.eventname + "_ctr", save=['data', 'uncd', 'mask'])

  # Print time elapsed and close log:
  cwd = os.getcwd() + "/"
  log.writelog("Output files (" + event.centerdir + "):")
  log.writelog("Data:")
  log.writelog(" " + cwd + event.eventname + "_ctr.dat")
  log.writelog(" " + cwd + event.eventname + "_ctr.h5")
  log.writelog("Log:")
  log.writelog(" " + cwd + event.logname)

  dt = t.hms_time(time.time()-tini)
  log.writeclose("\nEnd Centering. Time (h:m:s):  %s"%dt  +
                 "  (" + event.centerdir + ")")
  print("-------------  ------------\n")

  if hasattr(event, 'runp4') and event.runp4 == True:
    os.chdir(event.eventdir)
    os.system("poet.py p4 %s/"%event.centerdir)

def run_centering(eventname, control, cwd=None):
  """
  Read the control file.
  Load the event.
  Launch a thread for each centering run.
  """

  if cwd is None:
    cwd = os.getcwd()
  os.chdir(cwd)
  pcf = rd.read_pcf(control)
  nruns = len(pcf)

  if nruns == 1: #, I may be in the center dir, to re-run:
    # Get name of centering dir:
    centerdir = pcf[0].method.get()
    if pcf[0].pcfname.get() is not "":
      centerdir += "_" + pcf[0].pcfname.get()

    if cwd[-len(centerdir):] == centerdir:
      # Go to dir where poet2 files were saved.
      cwd = cwd[:-len(centerdir)]
      os.chdir(cwd)

  # Load the event:
  try:
    event = me.loadevent(eventname, load=['dendata', 'data','uncd','mask'])
    print("Performing centering on denoised data")
  except:
    event = me.loadevent(eventname, load=['data','uncd','mask'])
    event.denoised = False

  # Loop over each run:
  for run in np.arange(nruns):
    os.chdir(cwd)

    # Make a copy of the event:
    this_event = copy.copy(event)

    # Name of the directory to put the results:
    centerdir = pcf[run].method.get()
    if pcf[run].pcfname.get() is not "":
      centerdir += "_" + pcf[run].pcfname.get()
    this_event.centerdir = centerdir

    # Create the centering directory if it doesn't exist:
    if not os.path.exists(centerdir):
      os.mkdir(centerdir)

    # copy the photometry control file to centerdir
    shutil.copy('photom.pcf', centerdir + '/photom.pcf')

    # Launch the thread:
    p = Process(target=centering, args=(this_event, pcf[run], centerdir))
    p.start()


def do_center(start, end, event, centermask, log, x, y, sx, sy, flux, sky, goodfit):

  # Initialize a Timer to report progress:
  if start == 0:  # Only for the fisrt chunk
    clock = t.Timer(event.npos*end,
                    progress=np.array([0.05, 0.1, 0.2, 0.3, 0.4,  0.5,
                                       0.6,  0.7, 0.8, 0.9, 0.99, 1.1]))

  # Use denoised data if exists:
  if event.denoised:
    data = event.dendata
  else:
    data = event.data

  # Finally, do the centering:
  for  pos in np.arange(event.npos):
    # Recalculate star/end, Care not to go out of bounds:
    end   = np.amin([end,   event.nimpos[pos]])
    start = np.amin([start, event.nimpos[pos]]) # is this necessary?
    if event.noctr:   # Just use the mean x,y in this case
      y[pos*event.npos+start:pos*event.npos+end] = event.targpos[0, pos]
      x[pos*event.npos+start:pos*event.npos+end] = event.targpos[1, pos]
    else:
      for im in np.arange(start, end):
        # Index in the share memory arrays:
        ind = int(pos * event.npos + im)
        try:
          if event.weights:   # weight by uncertainties in fitting?
            uncd = event.uncd[im,:,:,pos]
          else:
            uncd = None
          # Do the centering:
          im = int(im)
          pos = int(pos)
          ind = int(ind)
          position, extra = cd.centerdriver(event.method, data[im,:,:,pos],
                                 event.targpos[:,pos], event.ctrim,
                                 event.cradius, event.csize,
                                 event.mask[im,:,:,pos]*centermask,
                                 uncd, fitbg=event.fitbg,
                                 expand=event.expand,
                                 psf=event.psfim, psfctr=event.psfctr)

          y[ind], x[ind] = position

          if event.method == "fgc":
            sy[ind] = extra[0]
            sx[ind] = extra[1]
          if event.method == "ipf" or event.method == "bpf":
            flux[ind] = extra[0]
            sky [ind] = extra[1]
            # FINDME: define some criterion for good/bad fit.
            goodfit[ind] = 1
        except:
          y[ind], x[ind] = event.targpos[:, pos]
          sy[ind], sx[ind] = 0.0, 0.0
          flux[ind], sky[ind] = 0.0, 0.0
          goodfit[ind] = 0
          log.writelog("Centering failed in im, pos: %5i"%im + ", %2i"%pos)

        if start == 0:
          # Report progress:
          clock.check(pos*end + im, name=event.centerdir)
