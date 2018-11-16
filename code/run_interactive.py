#! /usr/bin/env python

"""
#TO EXECUTE DIRECTLY FROM BASH, USE THE FOLLOWING SYNTAX:
./run_interactive.py p5idl filedir
./run_interactive.py p5 filedir topdir clip
./run_interactive.py p6 newdir filedir topdir clip idl
./run_interactive.py p7 filedir topdir idl islak
./run_interactive.py p8 filedir topdir idl eclphase
./run_interactive.py p9 filedir topdir idl
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

import os, sys
import poet_run as run

os.environ['OMP_NUM_THREADS']='4'
isinteractive = True       #Set False and don't use '-pylab' to turn off displaying plots in interactive mode.

#USE THE COMMANDS BELOW WHEN IN INTERACTIVE MODE
#ONLY COPY AND PASTE THE LINES THAT YOU NEED
def interactive():
    '''POET'''
    #Restore after p5checks
    #filedir:   Location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #clip:      Clip data set to model and plot only a portion of the entire light curve.
    #           Good for modeling a transit or eclipse in an around-the-orbit data set.
    #           Syntax is [start,end] (eg. [0,70000] or [-60500,None] to use last 60500 points).
    ev = run.poetRestore(filedir='..')
    
    #Optional: reduce the number of data points by binning
    #ev = run.binparams(ev, binsize=64)

    #Run p6model
    #ev:        If not specified, ev will be restored using poetRestore or p5idlRestore.
    #newdir:    Creates a new model directory when True (default).
    #           Can set to False or a string that appends to the default directory name.
    #filedir:   If event is not specified, location of the savefile to restore.
    #idl:       False (default), set True when using IDL photometry.
    run.p6model(ev)
    run.p6model(ev, newdir=False)
    run.p6model(ev, newdir='test')

    #Run p7anal, p8tables, and p9figs
    #ev:        If not specified, event will be restored using p6Restore.
    #filedir:   If event is not specified, location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #idl:       False (default), set True when using IDL photometry.
    #islak:     True (default), set False to not load 'allknots' into memory and skip Figs 701 & 702 (p7anal only).
    #eclphase:  Nominal eclipse phase, default is 0.5 (p8tables only).
    run.p7anal(ev)
    run.p8tables(ev)
    run.p9figs(ev)
    
    #Restore after p6model or p7anal, if necessary
    #You can use os.listdir('.') to find the desired model directory.
    #filedir:   Location of the savefile to restore.
    #topdir:    Specify POET top directory if it's not a parent of cwd.
    #idl:       False (default), set True when using IDL photometry.
    ev = run.p6Restore(filedir = '2012-08-09_13:52-2e6')
    ev = run.p6Restore(filedir='.')
    ev = run.p7Restore(filedir='.')
    
    #Overwrite the savefile for a single event after running poetRestore, p6model or p7anal.
    #Do not use these commands unless you are sure you want to overwrite an existing savefile.
    #filedir:   Location of the savefile to restore (poetSave only).
    #topdir:    Specify POET top directory if it's not a parent of cwd  (poetSave only).
    run.poetSave(ev[0], filedir='..')
    run.p6Save(ev[0])
    run.p7Save(ev[0])

    return


#MAIN IS CALLED WHEN EXECUTING DIRECTLY FROM BASH
def main(args):
    length = len(args)
    #Run p5
    elif args[1] == 'p5':
        ev = run.poetRestore(*args[2:])
    #Run p6model
    elif args[1] == 'p6':
        run.p6model(None, *args[2:])
    #Run p7anal
    elif args[1] == 'p7':
        run.p7anal(None, *args[2:])
    #Run p8tables
    elif args[1] == 'p8':
        run.p8tables(None, *args[2:])
    #Run p9figs
    elif args[1] == 'p9':
        run.p9figs(None, *args[2:])
    else:
        print("Unrecognized function.")
    return 0

#CALLS main IN BASH MODE THEN EXITS CLEANLY
if __name__ == '__main__':
    isinteractive = False
    sys.exit(main(sys.argv))

