
import numpy as np
'''
;+
; NAME:
;	TIME2PHASE
;
; PURPOSE:
;	This function converts a time to a phase in uniform, circular motion.
;
; CATEGORY:
;	Physics.
;
; CALLING SEQUENCE:
;
;	Result = TIME2PHASE(Times, Tzero, Period)
;
; INPUTS:
;	Times:	Array of times to convert.
;	Tzero:	A time, in the same system and units as that used for
;		Times, when the motion was at zero phase.
;	Period:	The duration of one cycle of the motion, same units
;		and Times.
;
; OUTPUTS:
;	This function returns a double-precision array or scalar of
;	the same shape as Times, giving the corresponding phases.  All
;	phases are in the range [0.d,1.].  -0.0 is sometimes seen.
;
; EXAMPLE:
;	print, time2phase([0.,1.2,3.3], 1.1, 1.)
;      0.89999998      0.10000002      0.19999981
;
; MODIFICATION HISTORY:
; 	Written by:	Joseph Harrington, Cornell.  2005 Oct 14
;			jh@oobleck.astro.cornell.edu
;   		Kevin Stevenson					2009 Mar 9
;			Converted to Python
;-
function time2phase, time, tzero, period

phase = ( (time - tzero) / period ) mod 1d
negs = where(phase lt 0d)
if negs[0] ne -1 then $
  phase[negs] += 1d

return, phase
end
'''
def time2phase(time, t0, period):
   return ((time - t0) / period) % 1.
   




