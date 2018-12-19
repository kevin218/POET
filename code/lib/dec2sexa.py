"""

   
"""

import numpy    # to allow input being an numpy array of floats


def dec2sexa1(dec, ndecdig=6):
   '''Convert a single scalar decimal values (degrees or hours) 
      to sexagesimal
      Returns tuple with (units,min,sec,usec)
      Units,min,sec are integer 
      fracsec is floating point of same type as dec'''
   
   if dec<0:
      dec=-dec
      thesign="-"
   else:
      thesign=" "
   
   if ndecdig<=0:
      ndecdig=0
      secw=2
   else:
      secw=3+ndecdig
   
   h= int(dec)   
   f = 60.*(dec-h)
   m = int(f)
   f = 60.*(f-m)
   s=round(f+ .5*(10**(-ndecdig-1)), ndecdig)  
   parts=(thesign,h,m,s)
   format="%%s%%2d:%%02d:%%0%d.%df"  % (secw, ndecdig)
   #print("TEST: format="+format)
   return format  % parts





def dec2sexa(dec, ndecdig=6):
   '''Convert decimal to sexagesimal, accepting scalars, lists or
      arrays'''
   d=dec
   if type(d)==type(1):
      return dec2sexa1(float(dec))
   elif type(d)==type(1.0):
      return dec2sexa1(dec,ndecdig=ndecdig)
   elif type(d)==type([]):
      return [dec2sexa1(x,ndecdig=ndecdig) for x in dec]
   elif type(d)==type(numpy.array([])):
      return [dec2sexa1(x,ndecdig=ndecdig) for x in dec]
   else:
      return None

