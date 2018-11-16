# $Author: patricio $
# $Revision: 301 $
# $Date: 2010-07-10 01:03:18 -0400 (Sat, 10 Jul 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/poet_4photom.py $
# $Id: poet_4photom.py 301 2010-07-10 05:03:18Z patricio $

import numpy as np
import time, os, copy
import apphot         as ap
import reader3        as rd
import poet_dooptphot as do
import logedit        as le
import timer          as t
import manageevent    as me
import pixlevdecorr as pld
from multiprocessing import Process, Array


def photometry(event, pcf, photdir, mute):
  tini = time.time()

  # Create photometry log
  logname = event.logname
  log = le.Logedit(photdir + "/" + logname, logname)
  log.writelog("\nStart " + photdir + " photometry: " + time.ctime())

  parentdir = os.getcwd() + "/"
  os.chdir(photdir)

  # copy photom.pcf in photdir
  pcf.make_file("photom.pcf")

  # Parse the attributes from the control file to the event:
  attrib = vars(pcf)
  keys = attrib.keys()
  for key in keys:
    setattr(event, key, attrib.get(key).get())

  maxnimpos, npos = event.maxnimpos, event.npos
  # allocating frame parameters:
  event.fp.aplev     = np.zeros((npos, maxnimpos)) # aperture flux
  event.fp.aperr     = np.zeros((npos, maxnimpos)) # aperture error
  event.fp.nappix    = np.zeros((npos, maxnimpos)) # number of aperture  pixels
  event.fp.skylev    = np.zeros((npos, maxnimpos)) # background sky flux level
  event.fp.skyerr    = np.zeros((npos, maxnimpos)) # sky error
  event.fp.nskypix   = np.zeros((npos, maxnimpos)) # number of sky pixels
  event.fp.nskyideal = np.zeros((npos, maxnimpos)) # ideal number of sky pixels
  event.fp.status    = np.zeros((npos, maxnimpos)) # apphot return status
  event.fp.good      = np.zeros((npos, maxnimpos)) # good flag
  event.fp.betaper   = np.zeros((npos, maxnimpos)) # beta aperture


  # Aperture photometry:
  if not event.dooptimal or event.from_aper == None:

    # Multy Process set up:
    # Shared memory arrays allow only 1D Arrays :(
    aplev     = Array("d", np.zeros(npos * maxnimpos))
    aperr     = Array("d", np.zeros(npos * maxnimpos))
    nappix    = Array("d", np.zeros(npos * maxnimpos))
    skylev    = Array("d", np.zeros(npos * maxnimpos))
    skyerr    = Array("d", np.zeros(npos * maxnimpos))
    nskypix   = Array("d", np.zeros(npos * maxnimpos))
    nskyideal = Array("d", np.zeros(npos * maxnimpos))
    status    = Array("d", np.zeros(npos * maxnimpos))
    good      = Array("d", np.zeros(npos * maxnimpos))
    betaper   = Array("d", np.zeros(npos * maxnimpos))
    # Size of chunk of data each core will process:
    chunksize = maxnimpos/event.ncores + 1

    print("Number of cores: " + str(event.ncores))
    # Start Muti Procecess:
    processes = []
    for nc in np.arange(event.ncores):
      start =  nc    * chunksize # Starting index to process
      end   = (nc+1) * chunksize # Ending   index to process
      proc = Process(target=do_aphot, args=(start, end, event, log, mute,
                                            aplev, aperr,
                                            nappix, skylev, skyerr, nskypix,
                                            nskyideal, status, good, betaper))
      processes.append(proc)
      proc.start()

    # Make sure all processes finish their work:
    for nc in np.arange(event.ncores):
      processes[nc].join()


    # Put the results in the event. I need to reshape them:
    event.fp.aplev     = np.asarray(aplev    ).reshape(npos,maxnimpos)
    event.fp.aperr     = np.asarray(aperr    ).reshape(npos,maxnimpos)
    event.fp.nappix    = np.asarray(nappix   ).reshape(npos,maxnimpos)
    event.fp.skylev    = np.asarray(skylev   ).reshape(npos,maxnimpos)
    event.fp.skyerr    = np.asarray(skyerr   ).reshape(npos,maxnimpos)
    event.fp.nskypix   = np.asarray(nskypix  ).reshape(npos,maxnimpos)
    event.fp.nskyideal = np.asarray(nskyideal).reshape(npos,maxnimpos)
    event.fp.status    = np.asarray(status   ).reshape(npos,maxnimpos)
    event.fp.good      = np.asarray(good     ).reshape(npos,maxnimpos)
    event.fp.betaper   = np.asarray(betaper  ).reshape(npos,maxnimpos)

    # raw photometry (no sky subtraction):
    event.fp.apraw = ( event.fp.aplev + ( event.fp.skylev * event.fp.nappix ) )

    # Print results into the log if it wans't done before:
    for  pos in np.arange(npos):
      for i in np.arange(event.nimpos[pos]):
        log.writelog('\nframe =%7d       '%i                 +
                       'pos   =%5d       '%pos               +
                       'y =%7.3f       '  %event.fp.y[pos,i] +
                       'x =%7.3f'         %event.fp.x[pos,i] + '\n' +
                       'aplev =%11.3f   ' %event.fp.aplev    [pos,i] +
                       'aperr =%9.3f   '  %event.fp.aperr    [pos,i] +
                       'nappix =%6.2f'    %event.fp.nappix   [pos,i] + '\n' +
                       'skylev=%11.3f   ' %event.fp.skylev   [pos,i] +
                       'skyerr=%9.3f   '  %event.fp.skyerr   [pos,i] +
                       'nskypix=%6.2f   ' %event.fp.nskypix  [pos,i] +
                       'nskyideal=%6.2f'  %event.fp.nskyideal[pos,i] + '\n' +
                       'status=%7d       '%event.fp.status   [pos,i] +
                       'good  =%5d'       %event.fp.good     [pos,i] + '\n' +
                       'betaper  =%5.3f'  %event.fp.betaper     [pos,i], mute=True)

  elif event.from_aper != None:
    # Load previous aperture photometry if required for optimal:
    evt = me.loadevent(parentdir + event.from_aper + "/" +
                       event.eventname + "_pht"          )
    event.fp.aplev     = evt.fp.aplev
    event.fp.aperr     = evt.fp.aperr
    event.fp.nappix    = evt.fp.nappix
    event.fp.skylev    = evt.fp.skylev
    event.fp.skyerr    = evt.fp.skyerr
    event.fp.nskypix   = evt.fp.nskypix
    event.fp.nskyideal = evt.fp.nskyideal
    event.fp.status    = evt.fp.status
    event.fp.good      = evt.fp.good
    event.fp.apraw     = evt.fp.apraw


  if event.dooptimal:
    ofp, psf = do.dooptphot(event.data, event.uncd, event.mask, event.fp,
                            event.srcest, event.nimpos, rejlim=[10.45,1000,1.5],
                            order=1, resize=event.oresize, norm=1,
                            trim=event.otrim, log=log )
    event.fp  = ofp
    event.psf = psf

  elif event.ispsf:
    # PSF aperture correction:
    log.writelog('Calculating PSF aperture:')

    event.aperfrac,    event.psfnappix,    event.psfskylev, \
    event.psfnskypix, event.psfnskyideal, event.psfstatus \
     = ap.apphot(image = event.psfim, ctr = event.psfctr,
                 photap = event.photap * event.psfexpand,
                 skyin = event.skyin  * event.psfexpand,
                 skyout = event.skyout * event.psfexpand,
                 betahw = event.betahw * event.psfexpand, 
                 targpos = event.targpos,
                 med    = event.skymed,
                 expand = event.apscale,
                 isbeta = event.isbeta,
                 nappix  = True, skylev    = True,
                 nskypix = True, nskyideal = True,
                 status  = True, betaper = False)

    event.aperfrac += event.psfskylev * event.psfnappix

    event.fp.aplev /= event.aperfrac
    event.fp.aperr /= event.aperfrac
    log.writelog('Aperture contains %f of PSF.'%event.aperfrac)

  # For running pixel-level decorrelation (pld)
  if event.ispld and event.npos == 1:
    event.apdata = pld.pld_box(event.data,
                               event.targpos,
                               event.pldhw,
                               event.fp.skylev)
    log.writelog("Created " + str(event.pldhw*2+1) + "x" + str(event.pldhw*2+1) + " box around centroid for pixel-level decorrelation and normalized it in time.")
  elif event.ispld and event.npos != 1:
    log.writelog("Could not perform pixel-level decorrelation because there is more than 1 nod position.")


  # save
  print("\nSaving ...")
  me.saveevent(event, event.eventname + "_pht", delete=['data', 'uncd', 'mask'])

  # Print time elapsed and close log:
  cwd = os.getcwd() + "/"
  log.writelog("Output files (" + event.photdir + "):")
  log.writelog("Data:")
  log.writelog(" " + cwd + event.eventname + "_pht.dat")
  log.writelog("Log:")
  log.writelog(" " + cwd + logname)

  dt = t.hms_time(time.time()-tini)
  log.writeclose("\nEnd Photometry. Time (h:m:s):  %s "%dt  +
                 "  (" + photdir + ")")
  print("--------------  ------------\n")


def run_photometry(eventname, control, cwd=None):
  """
  Load the event.
  Read the control file.
  Launch a thread for each centering run.
  """

  if cwd == None:
    cwd = os.getcwd()
  os.chdir(cwd)
  pcf = rd.read_pcf(control)
  nruns = len(pcf)


  if nruns == 1: #, I may be in photdir to re-run:
    # Get name of photometry dir:
    if pcf[0].dooptimal.get():
      photdir = 'optimal' + '%02d'%pcf[0].oresize.get()

    if pcf[0].isbeta.get():
          photdir = ('bap%03d'%(pcf[0].photap.get()*100) +
                     '%02d'%pcf[0].skyin.get() + '%02d'%pcf[0].skyout.get() )

    else:
      photdir = ('ap%03d'%(pcf[0].photap.get()*100) +
                 '%02d'%pcf[0].skyin.get() + '%02d'%pcf[0].skyout.get() )
    if pcf[0].pcfname.get() != "":
      photdir += "_" + pcf[0].pcfname.get()

    # If I am in the photometry dir already:
    if cwd[-len(photdir):] == photdir:
      # Go to dir where poet3 files were saved.
      cwd = cwd[:-len(photdir)]
      os.chdir(cwd)

    mute = False
  else:
    mute = True

  # Load the event:
  event = me.loadevent(eventname, load=['data','uncd','mask'])

  # Loop over each run:
  for run in np.arange(nruns):
    os.chdir(cwd)

    # Make a copy of the event:
    this_event = copy.copy(event)

    # Get name of photometry dir:
    if pcf[run].dooptimal.get():
      photdir = 'optimal' + '%02d'%pcf[run].oresize.get()

    if pcf[run].isbeta.get():
          photdir = ('bap%03d'%(pcf[run].photap.get()*100) +
                     '%02d'%pcf[run].skyin.get() + '%02d'%pcf[run].skyout.get() )
    else:
      photdir = ('ap%03d'%(pcf[run].photap.get()*100) +
                 '%02d'%pcf[run].skyin.get() + '%02d'%pcf[run].skyout.get() )
    if pcf[run].pcfname.get() != "":
      photdir += "_" + pcf[run].pcfname.get()
    this_event.photdir = photdir

    # Create the photometry directory if it doesn't exist:
    if not os.path.exists(photdir):
      os.mkdir(photdir)

    # Launch the thread:
    p = Process(target=photometry, args=(this_event, pcf[run], photdir, mute))
    p.start()



def do_aphot(start, end, event, log, mute, aplev, aperr, nappix, skylev, skyerr,
             nskypix, nskyideal, status, good, betaper):
  """
    Notes:
    ------
    Medium level routine that performs aperture photometry.
    Each thread from the main routine (photometry) will run do_aphot once.
    do_aphot stores the values in the shared memory arrays.
  """
  # Initialize a Timer to report progress (use first Process):
  if start == 0:
    clock = t.Timer(event.npos*end,
                    progress=np.array([0.05, 0.1, 0.25, 0.5, 0.75, 1.1]))

  for pos in np.arange(event.npos):
    # Recalculate star and end indexes. Care not to go out og bounds:
    end   = np.amin([end,   event.nimpos[pos]])
    start = np.amin([start, event.nimpos[pos]])

    for i in np.arange(start, end):
      # Index in the share memory arrays:
      loc = pos * event.npos + i
      # Calculate aperture photometry:
      pos = int(pos)
      i = int(i)
      aphot = ap.apphot(image = event.data[i, :, :, pos],
                        ctr = (event.fp.y[pos,i], event.fp.x[pos,i]),
                        photap = event.photap, skyin = event.skyin, skyout = event.skyout,
                        betahw = event.betahw, targpos = event.targpos,
                        mask  = event.mask[i, :, :, pos],
                        imerr = event.uncd[i, :, :, pos],
                        skyfrac = event.skyfrac, med = event.skymed,
                        expand  = event.apscale, isbeta = event.isbeta,
                        nochecks  = True, aperr  = True, nappix  = True,
                        skylev    = True, skyerr = True, nskypix = True,
                        nskyideal = True, status = True, betaper = True)

      # Store values:
      loc = int(loc)
      aplev  [loc], aperr  [loc], nappix   [loc], skylev[loc], \
      skyerr[loc], nskypix[loc], nskyideal[loc], status[loc], betaper[loc] = aphot

      if status[loc] == 0:
        good[loc] = 1 # good flag

      # Print to screen only if one core:
      if event.ncores == 1 and (not mute): # (end - start) == event.maxnimpos:
        if event.verbose:
          print('\nframe =%7d       '%i                 +
                         'pos   =%5d       '%pos               +
                         'y =%7.3f       '  %event.fp.y[pos,i] +
                         'x =%7.3f'         %event.fp.x[pos,i] + '\n' +
                         'aplev =%11.3f   ' %aplev     [loc] +
                         'aperr =%9.3f   '  %aperr     [loc] +
                         'nappix =%6.2f'    %nappix    [loc] + '\n' +
                         'skylev=%11.3f   ' %skylev    [loc] +
                         'skyerr=%9.3f   '  %skyerr    [loc] +
                         'nskypix=%6.2f   ' %nskypix   [loc] +
                         'nskyideal=%6.2f'  %nskyideal [loc] + '\n' +
                         'status=%7d       '%status    [loc] +
                         'good  =%5d'       %good      [loc] + '\n' +
                         'betaper  =%5.3f'  %betaper   [loc])


          perc = 100.0*(np.sum(event.nimpos[:pos])+i+1)/np.sum(event.nimpos)
          hms = clock.hms_left(np.sum(event.nimpos[0:pos]) + i)
          print("progress: %6.2f"%perc + "%  -  Remaining time (h:m:s):" + hms)

      # First Process when ncores > 1:
      if start == 0:
        if mute or event.ncores > 1:
          clock.check(pos*end + i, name=event.photdir)
