#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import os, sys, copy

os.chdir('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/')

import manageevent as me
os.environ['OMP_NUM_THREADS']='6'
sys.path.append('../')
sys.path.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/')
sys.path.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/models_c/ext_func/')
sys.path.append('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/models_c/py_func/')
from run import *
ancildir      = 'ancil/'
modeldir      = 'ancil/modelparams/'
isinteractive = True       #Set False and don't use '-pylab' to turn off displaying plots in interactive mode.

os.chdir('/Users/kristin/Documents/STScI/KBS_WFC3-IR_Algorithm/2017-07-15-w1_spec_width_20/fitbghw_20/whiteLC_07202017/')
ev = w4Restore(filedir='..', fname='d-wasp79b_1125_1650-w3.dat')
#ev = w4Restore(filedir='..', fname='d-wa062bph1_1150_1650-w3-WHITE.dat') # change file name
ev = ancil(ev) #look at this

w6model(ev, newdir='test')