#! /usr/bin/env python

# $Author: kevin $
# $Revision: 540 $
# $Date: 2011-08-13 00:08:23 -0400 (Sat, 13 Aug 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/pipeline/trunk/run.py $
# $Id: run.py 540 2011-08-13 04:08:23Z kevin $

"""
#TO EXECUTE DIRECTLY FROM BASH, USE THE FOLLOWING SYNTAX:
./run.py p5idl filedir
./run.py p5 filedir topdir clip
./run.py p6 newdir filedir topdir clip idl
./run.py p7 filedir topdir idl islak
./run.py p8 filedir topdir idl eclphase
./run.py p9 filedir topdir idl
#Arguments after p# are optional but must remain in order. 
#   eg. To specify topdir with p7, you MUST include filedir; 
#       however, specifying filedir does NOT require topdir or idl.
#newdir:    String to append to the default directory name.
            A new directory is ALWAYS created in BASH mode.
#filedir:   Location of the savefile to restore, default is '..' for p6 and '.' for others
#topdir:    Specify POET top directory if it's not a parent of cwd.
#           To specify default directory, type None.
#clip:      Clip data set to model and plot only a portion of the entire light curve.
#           Good for modeling a transit or eclipse in an around-the-orbit data set.
#           Syntax is start:end (eg. 0:70000 or -60500:None), defult is None for full data set.
#idl:       False (default), set True when using IDL photometry.
#islak:     True (default), set False to not load 'allknots' into memory and skip Figs 701 & 702.
#Figures will NOT display in BASH mode, but they will still save.
"""

from __future__ import print_function
import os, sys
import readeventhdf
os.environ['OMP_NUM_THREADS']='6'
sys.path.append('../')
from run import *
ancildir      = 'ancil/'
modeldir      = 'ancil/modelparams/'
isinteractive = True       #Set False and don't use '-pylab' to turn off displaying plots in interactive mode.

#USE THE COMMANDS BELOW WHEN IN INTERACTIVE MODE
#ONLY COPY AND PASTE THE LINES THAT YOU NEED
def interactive():
    '''IDL'''
    #Run p5checks or restore after p5checks
    #filedir:   Location of the savefile to restore.
    event = p5idl(filedir='.')
    event = p5idlRestore(filedir='.')
    '''POET'''
    #Restore after p5checks
    #filedir:   Location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #clip:      Clip data set to model and plot only a portion of the entire light curve.
    #           Good for modeling a transit or eclipse in an around-the-orbit data set.
    #           Syntax is [start,end] (eg. [0,70000] or [-60500,None] to use last 60500 points).
    event = poetRestore(filedir='..')

    #Run p6model
    #event:     If not specified, event will be restored using poetRestore or p5idlRestore.
    #newdir:    Creates a new model directory when True (default).
    #           Can set to False or a string that appends to the default directory name.
    #filedir:   If event is not specified, location of the savefile to restore.
    #idl:       False (default), set True when using IDL photometry.
    p6model(event)
    p6model(event, newdir=False)
    p6model(event, newdir='test')

    #Run p7anal, p8tables, and p9figs
    #event:     If not specified, event will be restored using p6Restore.
    #filedir:   If event is not specified, location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #idl:       False (default), set True when using IDL photometry.
    #islak:     True (default), set False to not load 'allknots' into memory and skip Figs 701 & 702 (p7anal only).
    #eclphase:  Nominal eclipse phase, default is 0.5 (p8tables only).
    p7anal(event)
    p8tables(event)
    p9figs(event)
    
    #Restore after p6model or p7anal, if necessary
    #You can use os.listdir('.') to find the desired model directory.
    #filedir:   Location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #idl:       False (default), set True when using IDL photometry.
    event = p6Restore(filedir = '2012-08-09_13:52-2e6')
    event = p6Restore(filedir='.')
    event = p7Restore(filedir='.')
    
    #Overwrite the savefile for a single event after running poetRestore, p6model or p7anal.
    #Do not use these commands unless you are sure you want to overwrite an existing savefile.
    #filedir:   Location of the savefile to restore (poetSave only).
    #topdir:    Specify POET top directory if it's not a parent of cwd  (poetSave only).
    poetSave(event[0], filedir='..')
    p6Save(event[0])
    p7Save(event[0])

    return

########################################################
#                                                      #
#   DO NOT RUN CODE INTERACTIVELY BEYOND THIS POINT!   #
#                                                      #
########################################################

#RUN p5checks AFTER IDL PHOTOMETRY
def p5idl(filedir='.'):
    import p5checks as p5
    global numevents
    #Read in HDF5 file
    hdf5file  = []
    for fname in os.listdir(filedir):
        if (fname.endswith(".h5")):
            hdf5file.append(fname)
    hdf5file.sort()
    numevents = len(hdf5file)
    print('Using the following HDF5 files:', hdf5file)
    #Run p5checks
    event     = []
    filename  = ''
    for i in range(numevents):
        event.append(p5.checks(hdf5file[i], i))
        event[i].ancildir   = ancildir
        exec 'import ' + event[i].eventname + 'params as op' + str(i)
        filename = ''.join((filename,event[i].eventname))
    #p5checks - save
    for i in range(numevents):
        event[i].filename = filename
        savefile  = "d-" + event[i].eventname + "-5checks.dat"
        handle    = open(savefile, 'w')
        pickle.dump(event[i], handle)
        handle.close()
    return event

#RESTORE SAVEFILE AFTER p5checks USING IDL PHOTOMETRY
def p5idlRestore(filedir='.'):
    global numevents
    loadfile = []
    for fname in os.listdir(filedir):
        if (fname.endswith("5checks.dat")):
            loadfile.append(fname)
    loadfile.sort()
    numevents = len(loadfile)
    event    = []
    for lfile in loadfile:
        print("Loading " + lfile)
        handle      = open(lfile, 'r')
        try:
            event.append(pickle.load(handle))
        except:
            event.append(unpic.unpickle_old_pyfits(lfile))
        handle.close()
    return event

#RESTORE SAVEFILE FROM POET
def poetRestore(filedir='..', topdir=None, clip=None):
    import shutil
    global numevents
    #Append system path
    if topdir == None or topdir == 'None':
        r = os.getcwd().split("/")
        topdir = "/".join(r[:r.index("run")])
    sys.path.append(topdir + '/lib/')
    import manageevent as me
    
    files    = []
    event    = []
    filename = ''
    for fname in os.listdir(filedir):
        if (fname.endswith("_p5c.dat")):
            files.append(fname[:-4])
    files.sort()
    numevents = len(files)
    if numevents == 0:
        print('Cannot find any files to restore.')
        return []
    for i in np.arange(numevents):
        #Load event
        event.append(me.loadevent(filedir+'/'+files[i]))
        print('Finished loading: ' + event[i].eventname)
        filename = ''.join((filename,event[i].eventname))
        event[i].ancildir   = ancildir
        #Clip data set to model and plot only a portion of the entire light curve
        #Good for modeling a transit or eclipse in an around-the-orbit data set
        if clip != None and clip != 'None':
            if type(clip) == str:
                #Convert from string to 2 ints
                start, end = clip.split(':',1)
                try:
                    start = int(start)
                except:
                    print("Error with format of optional variable clip.")
                    return []
                try:
                    end   = int(end)
                except:
                    end   = None
            else:
                if len(clip) == 2:
                    start, end = clip
                else:
                    start = clip[0]
                    end   = None
            #Use only data points from 'start' to 'end'
            event[i].phase      = event[i].phase[:,start:end]
            event[i].aplev      = event[i].aplev[:,start:end]
            event[i].aperr      = event[i].aperr[:,start:end]
            event[i].good       = event[i].good[:,start:end]
            event[i].time       = event[i].time[:,start:end]
            event[i].y          = event[i].y[:,start:end]
            event[i].x          = event[i].x[:,start:end]
            event[i].juldat     = event[i].juldat[:,start:end]
            event[i].bjdutc     = event[i].bjdutc[:,start:end]
            event[i].bjdtdb     = event[i].bjdtdb[:,start:end]
    #Create and populate ancil directory, if it doesn't already exist
    if os.path.isdir(ancildir) == False:
        os.mkdir(ancildir, 0775)
        sys.path.append(ancildir)
        for i in np.arange(numevents):
            #Copy params file into new ancil dir.
            paramsfile = event[i].topdir + modeldir + event[i].eventname + 'params.py'
            event[i].paramsfile = ancildir + event[i].eventname + 'params.py'
            if os.path.isfile(event[i].paramsfile) == False:
                if os.path.isfile(paramsfile):
                    shutil.copy(paramsfile, event[i].paramsfile)
                else:
                    shutil.copy(event[i].topdir + modeldir + 'eg00params.py', event[i].paramsfile)
            #Copy initial parameters file into new ancil dir
            initpfile = []
            for f in os.listdir(event[i].topdir + modeldir):
                if f.startswith(event[i].eventname) and f.endswith('.txt'):
                    initpfile.append(f)
            if len(initpfile) == 0:
                shutil.copy(event[i].topdir + modeldir + 'eg00-initvals.txt', ancildir)
            for f in initpfile:
                if os.path.isfile(ancildir + f) == False:
                    shutil.copy(event[i].topdir + modeldir + f, ancildir)
    else:
        #On initial setup, rename eg00params and eg00-initvals
        for i in np.arange(numevents):
            event[i].paramsfile   = ancildir + event[i].eventname + 'params.py'
            event[i].initvalsfile = ancildir + event[i].eventname + '-initvals.txt'
            # Copy eg00params
            if os.path.isfile(event[i].paramsfile) == False:
                shutil.copy(ancildir + 'eg00params.py', event[i].paramsfile)
            # Copy eg00-initvals
            if os.path.isfile(event[i].initvalsfile) == False:
                shutil.copy(ancildir + 'eg00-initvals.txt', event[i].initvalsfile)
    return event

#OVERWRITE event SAVEFILE AFTER RUNNING p5checks in POET
def poetSave(event, filedir='..', topdir=None):
    #Append system path
    if topdir == None:
        r = os.getcwd().split("/")
        topdir = "/".join(r[:r.index("run")])
    sys.path.append(topdir + '/lib/')
    import manageevent as me
    me.saveevent(event, filedir+"/"+event.eventname + "_p5c")

#RUN p6model
def p6model(event=None, newdir=True, filedir='..', topdir=None, clip=None, idl=False):
    import datetime
    import p6model as p6
    reload(p6)
    global numevents, nummodels, isinteractive
    if event == None:
        if idl == False:
            event = poetRestore(filedir, topdir, clip)
        else:
            event = p5idlRestore(filedir)
    if len(event) == 0:
        print("Event object is empty.")
        return

    if newdir == True:
        #CREATE NEW DIRECTORY FOR MODEL RUN
        modeldir = datetime.datetime.strftime(datetime.datetime.today(), "%Y-%m-%d_%H:%M")
        os.mkdir(modeldir)
    elif type(newdir) == str:
        #CREATE NEW DIRECTORY FOR MODEL RUN WITH CUSTOM NAME
        modeldir = datetime.datetime.strftime(datetime.datetime.today(), "%Y-%m-%d_%H:%M")+"-"+newdir
        os.mkdir(modeldir)
    else:
        try:
            modeldir = event[0].modeldir
        except:
            print("Model directory has not been specified.")
            return

    #RELOAD MODEL PARAMETERS
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        exec('import ' + event[j].eventname + 'params as op' + str(j))
        exec("reload(op" + str(j) + ")")
        event[j].params = readeventhdf.initParams()
        exec("op" + str(j) + ".modelparams(event["+ str(j) + "].params)")
        #exec("event[" + str(j) + "].params = op" + str(j) + ".params")
        event[j].fit = []
        event[j].modeldir = modeldir
        if hasattr(event[j].params, 'noisewavelet') and event[j].params.noisewavelet[0] == 'all':
            event[j].params.noisewavelet = pywt.wavelist()
            nummodels[j] = len(pywt.wavelist())
            event[j].params.model = [event[j].params.model[0] for i in range(nummodels[j])]
        else:
            nummodels[j] = len(event[j].params.model)
        if j > 0 and nummodels[j] != nummodels[j-1]:
            print("WARNING: Number of models in each event does not match.")

    #INITIALIZE OUTPUT TYPE: stdout, A FILE OBJECT OR A FILE
    printout = printoutput.init(event[0].params.printout, event)
    
    #Execute rundmc
    for i in range(nummodels.min()):
        p6.rundmc(event, i, printout, isinteractive=isinteractive)
        if hasattr(event[0].params, 'savedata') and event[0].params.savedata == False:
            pass
        else:
            for j in range(numevents):
                p6Save(event[j])

    #PRINT PARAMETERS USED FOR COMPARISON
    print("\nBest-fit eclipse depths or transit radius ratios with errors:", file=printout)
    for j in range(numevents):
        event[j].minbic = np.inf     #Minimum BIC value of all fits for one event.
        print(event[j].eventname, file=printout)
        for i in range(len(event[j].fit)):
            if hasattr(event[j].fit[i].i,'depth'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.depth  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.depth,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'depth2'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.depth2  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.depth2,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'depth3'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.depth3  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.depth3,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'trqrprs'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.trqrprs  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.trqrprs,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'trq2rprs'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.trq2rprs  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.trq2rprs,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'trrprs'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.trrprs  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.trrprs,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'trrprs2'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.trrprs2  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.trrprs2,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'rprs'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.rprs  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.rprs,1],
                      event[j].fit[i].saveext, file=printout)
            if hasattr(event[j].fit[i].i,'rprs2'):
                print(event[j].fit[i].bestp  [event[j].fit[i].i.rprs2  ], 
                      event[j].fit[i].medianp[event[j].fit[i].i.rprs2,1],
                      event[j].fit[i].saveext, file=printout)
            event[j].minbic = np.min((event[j].minbic,event[j].fit[i].bic))

    print("\n     S/N      SDNR     \xce\x94BIC       MODEL   NUMIT  BIN_SZ(y,x)   MinNumPts     Wavelet", file=printout)
    #Delta = '\xce\x94' in utf-8
    for j in range(numevents):
        print(event[j].eventname, file=printout)
        minbic = event[j].minbic
        for i in range(len(event[j].fit)):
            try:
                sdnr  = event[j].fit[i].sdnr
                bic   = event[j].fit[i].bic - minbic
                model = event[j].fit[i].saveext
                numit = event[j].params.numit[1]
                if len(event[j].params.ystep) == len(event[j].params.model):
                    ystep, xstep = event[j].params.ystep[i], event[j].params.xstep[i]
                else:
                    ystep, xstep = event[j].params.ystep[0], event[j].params.xstep[0] 
                if len(event[j].params.minnumpts) == len(event[j].params.model):
                    minnumpts = event[j].params.minnumpts[i] 
                else:
                    minnumpts = event[j].params.minnumpts[0]
                if event[j].params.model[i].__contains__('rednoise'):
                    if len(event[j].params.noisewavelet) == len(event[j].params.model):
                        wavelet = event[j].params.noisewavelet[i]
                    else:
                        wavelet = event[j].params.noisewavelet[0]
                else:
                    wavelet = 'None'
                if hasattr(event[j].fit[i].i,'depth'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.depth  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.depth,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'depth2'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.depth2  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.depth2,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'depth3'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.depth3  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.depth3,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'trrprs'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.trrprs  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.trrprs,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'trrprs2'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.trrprs2  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.trrprs2,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'trqrprs'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.trqrprs  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.trqrprs,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'trq2rprs'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.trq2rprs  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.trq2rprs,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'rprs'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.rprs  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.rprs,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
                if hasattr(event[j].fit[i].i,'rprs2'):
                    snr   = event[j].fit[i].bestp  [event[j].fit[i].i.rprs2  ] / \
                            event[j].fit[i].medianp[event[j].fit[i].i.rprs2,1]
                    print('%8.4f %9.7f %8.1f %11s %7.1e %6.3f,%5.3f %4.0f %12s' % 
                         (snr, sdnr, bic, model, numit, ystep, xstep, minnumpts, wavelet), file=printout)
            except:
                print("Error calculating values. %13s" % event[j].fit[i].saveext, file=printout)
    printoutput.close(printout)
    if isinteractive == False:
        plt.close('all')
    else:
        plt.show()
    return

#SAVE event AFTER p6model
def p6Save(event):
    cwd = os.getcwd().split("/")
    if cwd[-1] == event.modeldir:
        savefile  = "d-" + event.eventname + "-6model.dat"
    else:
        savefile  = event.modeldir + "/d-" + event.eventname + "-6model.dat"
    handle    = open(savefile, 'w')
    pickle.dump(event, handle)
    handle.close()
    return

#RESTORE SAVEFILE AFTER p6model
def p6Restore(filedir='.', topdir=None, idl=False):
    global numevents, nummodels
    '''
    if idl == False:
        #Append system path
        if topdir == None:
            r = os.getcwd().split("/")
            topdir = "/".join(r[:r.index("run")])
        sys.path.append(topdir + '/lib/')
    '''
    loadfile = []
    for fname in os.listdir(filedir):
        if (fname.endswith("6model.dat")):
            loadfile.append(filedir+"/"+fname)
    loadfile.sort()
    #print("Loading files:", loadfile)
    #loadfile = [loadfile[-1]]   #***EDIT THIS LINE MANUALLY***
    numevents = len(loadfile)
    if numevents == 0:
        print('Cannot find any files to restore.')
        return []
    event    = []
    for lfile in loadfile:
        print("Loading " + lfile)
        handle      = open(lfile, 'r')
        try:
            event.append(pickle.load(handle))
        except:
            event.append(unpic.unpickle_old_pyfits(lfile))
        handle.close()
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        nummodels[j] = len(event[j].params.model)
    return event

#RUN 7anal
def p7anal(event=None, filedir='.', topdir=None, idl=False, islak=True):
    import p7anal as p7
    reload(p7)
    global numevents, nummodels, isinteractive
    if event == None:
        event = p6Restore(filedir, topdir, idl)
    cwd = os.getcwd().split("/")
    if cwd[-1] == event[0].modeldir:
        os.chdir('..')
    printout = printoutput.init(event[0].params.printout, event)
    for j in range(numevents):
        print("\n" + event[j].eventname, file=printout)
        for i in range(nummodels.min()):
            print("\nCurrent model = " + str(event[j].params.model[i]), file=printout)
            p7.stdanal(event[j], event[j].fit[i], j*nummodels.min()+i, printout, islak)
            p7Save(event[j])
    printoutput.close(printout)
    if isinteractive == False:
        plt.close('all')
    else:
        plt.show()
    return

#SAVE event AFTER p7anal
def p7Save(event):
    cwd = os.getcwd().split("/")
    if cwd[-1] == event.modeldir:
        savefile  = "d-" + event.eventname + "-7anal.dat"
    else:
        savefile  = event.modeldir + "/d-" + event.eventname + "-7anal.dat"
    handle    = open(savefile, 'w')
    pickle.dump(event, handle)
    handle.close()
    return

#RESTORE SAVEFILE AFTER p7anal
def p7Restore(filedir='.', topdir=None, idl=False):
    global numevents, nummodels
    if idl == False:
        #Append system path
        if topdir == None:
            r = os.getcwd().split("/")
            topdir = "/".join(r[:r.index("run")])
        sys.path.append(topdir + '/lib/')
    loadfile = []
    for fname in os.listdir(filedir):
        if (fname.endswith("7anal.dat")):
            loadfile.append(filedir+"/"+fname)
    loadfile.sort()
    numevents = len(loadfile)
    if numevents == 0:
        print('Cannot find any files to restore.')
        return []
    event    = []
    for lfile in loadfile:
        print("Loading " + lfile)
        handle      = open(lfile, 'r')
        try:
            event.append(pickle.load(handle))
        except:
            event.append(unpic.unpickle_old_pyfits(lfile))
        handle.close()
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        nummodels[j] = len(event[j].params.model)
    return event

#RUN p8tables
def p8tables(event=None, filedir='.', topdir=None, idl=False, eclphase=0.5):
    import p8tables as p8
    reload(p8)
    global numevents, nummodels
    if event == None:
        event = p7Restore(filedir, topdir, idl)
    cwd = os.getcwd().split("/")
    if cwd[-1] == event[0].modeldir:
        os.chdir('..')
    printout = printoutput.init(event[0].params.printout, event)
    for j in range(numevents):
        print("\n" + event[j].eventname, file=printout)
        event[j].meanphase = eclphase
        for i in range(nummodels.min()):
            print("\nCurrent model = " + str(event[j].params.model[i]), file=printout)
            p8.tables(event[j], i, printout)
    printoutput.close(printout)
    return

#RUN p9figs
def p9figs(event=None, filedir='.', topdir=None, idl=False):
    import p9figs as p9
    reload(p9)
    global numevents, nummodels, isinteractive
    if event == None:
        event = p7Restore(filedir, topdir, idl)
    cwd = os.getcwd().split("/")
    if cwd[-1] == event[0].modeldir:
        os.chdir('..')
    printout = printoutput.init(event[0].params.printout, event)
    for j in range(numevents):
        print("\n" + event[j].eventname, file=printout)
        for i in range(nummodels.min()):
            print("\nCurrent model = " + str(event[j].params.model[i]), file=printout)
            p9.figs  (event[j], i, j*nummodels.min()+i)
    printoutput.close(printout)
    if isinteractive == False:
        plt.close('all')
    else:
        plt.show()
    return

import cPickle as pickle
import matplotlib.pyplot as plt
import np_unpickle as unpic
import numpy as np
import printoutput, readeventhdf, pywt
sys.path.append(ancildir)
sys.path.append('../'+ancildir)

#MAIN IS CALLED WHEN EXECUTING DIRECTLY FROM BASH
def main(args):
    length = len(args)
    #Run p5idl
    if   args[1] == 'p5idl':
        p5idl(*args[2:])
    #Run p5
    elif args[1] == 'p5':
        event = poetRestore(*args[2:])
    #Run p6model
    elif args[1] == 'p6':
        p6model(None, *args[2:])
    #Run p7anal
    elif args[1] == 'p7':
        p7anal(None, *args[2:])
    #Run p8tables
    elif args[1] == 'p8':
        p8tables(None, *args[2:])
    #Run p9figs
    elif args[1] == 'p9':
        p9figs(None, *args[2:])
    else:
        print("Unrecognized function.")
    return 0

#CALLS main IN BASH MODE THEN EXITS CLEANLY
if __name__ == '__main__':
    isinteractive = False
    sys.exit(main(sys.argv))

