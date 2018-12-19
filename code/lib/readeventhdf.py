# Read event's HDF5 using h5py, returning a python object
# All arrays, variables are "dot" attributes, like filed in C structs
# Returned object is an instance of AnEvent.
# HDF groups, which were originally structs in IDL, are made 
# into instances of DataGob with members as attributes, nested 
# as deep as needed.

# 06-28-2009	Kevin Stevenson
#				Updated h5py command to list members in savefile


import h5py
import numpy

class initParams:
   def __init__(self):
      pass

class indices: 
	def __init__(self):
		self.size = 0

class fits:
	def __init__(self):
		self.i     = indices()

class AnEvent:
   def __init__(self):
      pass




class DataGob:
   def __init__(self):
      pass
   


def _undress(e):
   #print("----\ndigesting item of class %s named %s"  %(e.__class__, e.name))
   if e.__class__ == h5py.Group:
      g=DataGob()
      for n in e.names:
         g.__dict__[n]=_undress(e[n])
      return g
   else:
      return e.value



def ReadEventHDF(fname):
   hdf=h5py.File(fname, 'r')
   E=AnEvent()
   #for n in hdf.names:
   for n in list(hdf):
      E.__dict__[n]=_undress(hdf[n])
   return E
   
