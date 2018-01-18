"""
 NAME:
	ECLIPSETIMING

 PURPOSE:
	This function  algorithm

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = 

 INPUTS:
    	e:	eccentricity
	w:	arguement of pariapsis
	T:	Orbital period

 OUTPUTS:
	This function returns an array of the best fitting parameters,
	an array of all parameters over all iterations, and numaccept.

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-06-30
			kbstevenson@gmail.com
"""

def eclipsetiming(e, w, T):

   import numpy as np

#GJ436
e  = 0.15
w  = 351*np.pi/180       #radians
RA = 11.703056*np.pi/180
dec = 26.706389*np.pi/180
#phase = 0.58

#HAT-P-2
e = 0.5163
w = 189.92*np.pi/180       #radians
#Phase = 0.18

#HD17156
e = 0.6717
w = 121.23*np.pi/180
#Phase = ?? 0.33 ??

#HD80606
e = 0.93
w = 300*np.pi/180
#Phase = 

#COMPUTE PRIMARY AND SECONDARY TRUE ANOMALLIES
#FINDME: Incorrect!!!
fp = np.pi/2 - w  #-2*np.pi
fs = w + np.pi/2
#COMPUTE ECCENTRIC ANOMALLIES
Ep = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(fp/2))
Es = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(fs/2))
#COMPUTE MEAN ANOMOLIES
Mp = Ep-e*np.sin(Ep)
Ms = Es-e*np.sin(Es)
#COMPUTE PHASE W.R.T. SEMI-MAJOR AXIS
tp = Mp/np.pi
ts = Ms/np.pi
tp
ts
1 - abs(ts) - abs(tp)

   return ts - tp

