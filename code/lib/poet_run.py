
import numpy as np
import matplotlib.pyplot as plt
import os, sys, shutil, datetime, pickle
import printoutput, readeventhdf, pywt
import manageevent as me
from importlib import reload

ancildir      = 'ancil/'
sys.path.append('../')
sys.path.append(ancildir)
sys.path.append('../'+ancildir)

#RESTORE SAVEFILE FROM POET
def poetRestore(filedir='..', topdir=None, clip=None):
    #global numevents
    #Append system path
    if topdir == None or topdir == 'None':
        r = os.getcwd().split("/")
        topdir = "/".join(r[:r.index("run")])
    sys.path.append(topdir + '/lib/')
    
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
            event[i].sy         = event[i].sy[:,start:end]
            event[i].sx         = event[i].sx[:,start:end]
            event[i].juldat     = event[i].juldat[:,start:end]
            event[i].bjdutc     = event[i].bjdutc[:,start:end]
            event[i].bjdtdb     = event[i].bjdtdb[:,start:end]
    #Create and populate ancil directory, if it doesn't already exist
    if os.path.isdir(ancildir) == False:
        os.mkdir(ancildir, 775)
        sys.path.append(ancildir)
        for i in np.arange(numevents):
            #Copy params file into new ancil dir.
            paramsfile = event[i].topdir + modeldir + event[i].eventname + '_params.py'
            event[i].paramsfile = ancildir + event[i].eventname + '_params.py'
            if os.path.isfile(event[i].paramsfile) == False:
                if os.path.isfile(paramsfile):
                    shutil.copy(paramsfile, event[i].paramsfile)
                else:
                    shutil.copy(event[i].topdir + modeldir + 'eg00_params.py', event[i].paramsfile)
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
            event[i].paramsfile   = ancildir + event[i].eventname + '_params.py'
            event[i].initvalsfile = ancildir + event[i].eventname + '-initvals.txt'
            # Copy eg00params
            if os.path.isfile(event[i].paramsfile) == False:
                print("Missing file: "+event[i].paramsfile)
                try:
                    shutil.copy(ancildir + 'eg00_params.py', event[i].paramsfile)
                except:
                    print("Missing file: "+ancildir + 'eg00_params.py')
            # Copy eg00-initvals
            if os.path.isfile(event[i].initvalsfile) == False:
                print("Missing file: "+event[i].initvalsfile)
                try:
                    shutil.copy(ancildir + 'eg00-initvals.txt', event[i].initvalsfile)
                except:
                    print("Missing file: "+ancildir + 'eg00-initvals.txt')
    return event

#OVERWRITE event SAVEFILE AFTER RUNNING p5checks in POET
def poetSave(event, filedir='..', topdir=None):
    #Append system path
    if topdir == None:
        r = os.getcwd().split("/")
        topdir = "/".join(r[:r.index("run")])
    sys.path.append(topdir + '/lib/')
    me.saveevent(event, filedir+"/"+event.eventname + "_p5c")

#RUN p6model
def p6model(event=None, newdir=True, filedir='..', topdir=None, clip=None, idl=False, isinteractive=True):
    """
    
    """
    import p6model as p6
    reload(p6)
    #global numevents, nummodels, isinteractive
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
    numevents   = len(event)
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        exec('import ' + event[j].eventname + '_params as op' + str(j))
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
    '''
    
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
        event.append(pickle.load(handle))
        handle.close()
    nummodels = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        nummodels[j] = len(event[j].params.model)
    return event

#RUN 7anal
def p7anal(event=None, filedir='.', topdir=None, idl=False, islak=True):
    import p7anal as p7
    #reload(p7)
    #global numevents, nummodels, isinteractive
    if event == None:
        event = p6Restore(filedir, topdir, idl)
    cwd = os.getcwd().split("/")
    if cwd[-1] == event[0].modeldir:
        os.chdir('..')
    printout    = printoutput.init(event[0].params.printout, event)
    numevents   = len(event)
    nummodels   = np.zeros(numevents,dtype=int)
    for j in range(numevents):
        nummodels[j] = len(event[j].params.model)
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
    #global numevents, nummodels
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
        event.append(pickle.load(handle))
        handle.close()
    #nummodels = np.zeros(numevents,dtype=int)
    #for j in range(numevents):
    #    nummodels[j] = len(event[j].params.model)
    return event

#RUN p8tables
def p8tables(event=None, filedir='.', topdir=None, idl=False, eclphase=0.5):
    import p8tables as p8
    #reload(p8)
    #global numevents, nummodels
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
    #reload(p9)
    #global numevents, nummodels, isinteractive
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

def binparams(ev, binsize=64):
    '''
    Reduce the number of data points by binning the flux, time, uncertainties, positions, etc.
    '''
    for i in range(len(ev)):
        try:
            del(ev[i].bjdutc,ev[i].apdata)
        except:
            pass
        nobj    = len(ev[i].phase[0])
        nbins   = int(np.ceil(nobj/binsize))
        phase   = np.zeros(nbins)
        bjdtdb  = np.zeros(nbins)
        time    = np.zeros(nbins)
        bjdcor  = np.zeros(nbins)
        juldat  = np.zeros(nbins)
        aplev   = np.zeros(nbins)
        aperr   = np.zeros(nbins)
        x       = np.zeros(nbins)
        y       = np.zeros(nbins)
        sx      = np.zeros(nbins)
        sy      = np.zeros(nbins)
        good    = np.ones(nbins, dtype=int)
        pos     = np.zeros(nbins)
        frmvis  = np.zeros(nbins)
        for j in range(nbins):
            start       = int(1.*j*nobj/nbins)
            end         = int(1.*(j+1)*nobj/nbins)
            isgood      = np.where(ev[i].good[0,start:end])[0] + start
            if len(isgood) == 0:
                good[j] = 0
            else:
                phase[j]    = np.mean(ev[i].phase[0,isgood])
                bjdtdb[j]   = np.mean(ev[i].bjdtdb[0,isgood])
                time[j]     = np.mean(ev[i].time[0,isgood])
                bjdcor[j]   = np.mean(ev[i].bjdcor[0,isgood])
                juldat[j]   = np.mean(ev[i].juldat[0,isgood])
                aplev[j]    = sum(ev[i].aplev[0,isgood]/ev[i].aperr[0,isgood]**2)/ \
                              sum(1/ev[i].aperr[0,isgood]**2)
                aperr[j]    = np.sqrt(1 / sum(1/ev[i].aperr[0,isgood]**2))
                x[j]        = np.mean(ev[i].x[0,isgood])
                y[j]        = np.mean(ev[i].y[0,isgood])
                sx[j]       = np.mean(ev[i].sx[0,isgood])
                sy[j]       = np.mean(ev[i].sy[0,isgood])
                #good[j]     = np.min(ev[i].good[0,isgood])
                pos[j]      = np.mean(ev[i].pos[0,isgood])
                frmvis[j]   = np.mean(ev[i].frmvis[0,isgood])
        #Make copy of parameters
        #ev[i].phase_old = np.copy(ev[i].phase)
        #Replace parameters
        ev[i].phase     = np.copy(phase)
        ev[i].bjdtdb    = np.copy(bjdtdb)
        ev[i].time      = np.copy(time)
        ev[i].bjdcor    = np.copy(bjdcor)
        ev[i].juldat    = np.copy(juldat)
        ev[i].aplev     = np.copy(aplev)
        ev[i].aperr     = np.copy(aperr)
        ev[i].x         = np.copy(x)
        ev[i].y         = np.copy(y)
        ev[i].sx        = np.copy(sx)
        ev[i].sy        = np.copy(sy)
        ev[i].good      = np.copy(good)
        ev[i].pos       = np.copy(pos)
        ev[i].frmvis    = np.copy(frmvis)
    return ev

