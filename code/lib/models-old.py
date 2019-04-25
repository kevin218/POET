import numpy as np
import scipy.interpolate as spi
import numexpr as ne
import models
import orbit

def setupmodel(model, ind):
   """

 NAME:
	SETUPMODEL

 PURPOSE:
	This function sets up the indeces and parameter names for the given model, also defined in this file.

 FUNCTYPE = 0: eclipse model function
            1: ramp model function
            2: intra-pixel model function
            3: position function
            4: position sensitivity function
            5: flat field correction
            6: Interpolated intra-pixel model function

   """
   
   parname  = []
   myfuncs  = []
   saveext = []
   functype = np.zeros(len(model))
   for i in range(len(model)):
   	if   model[i] == 'mandelecl':
		#DEFINE INDICES
		ind.midpt = ind.size
		ind.width = ind.size + 1
		ind.depth = ind.size + 2
		ind.t12   = ind.size + 3
		ind.t34   = ind.size + 4
		ind.flux  = ind.size + 5
		ind.size += 6
		#DEFINE NAMES
		parname.insert(ind.midpt,'Eclipse Phase')
		parname.insert(ind.width,'Eclipse Duration')
		parname.insert(ind.depth,'Eclipse Flux Ratio')
		parname.insert(ind.t12,  'Ingress Time')
		parname.insert(ind.t34,  'Egress Time')
		parname.insert(ind.flux,r'System Flux ($\mu Jy$)')
		#DEFINE ECLIPSE MODEL
		myfuncs.append(mandelecl)
		saveext.append('m1')
		functype[i] = 0
   	elif model[i] == 'mandelecl2':
		#DEFINE INDICES
		ind.midpt2 = ind.size
		ind.width2 = ind.size + 1
		ind.depth2 = ind.size + 2
		ind.t122   = ind.size + 3
		ind.t342   = ind.size + 4
		ind.flux2  = ind.size + 5
		ind.size  += 6
		#DEFINE NAMES
		parname.insert(ind.midpt2,'Eclipse Phase')
		parname.insert(ind.width2,'Eclipse Duration')
		parname.insert(ind.depth2,'Eclipse Flux Ratio')
		parname.insert(ind.t122,  'Ingress Time')
		parname.insert(ind.t342,  'Egress Time')
		parname.insert(ind.flux2,r'System Flux ($\mu Jy$)')
		#DEFINE ECLIPSE MODEL
		myfuncs.append(mandelecl)
		saveext.append('m2')
		functype[i] = 0
	elif model[i] == 'mandelgeom':
		#DEFINE INDICES
		ind.midpt3 = ind.size
		ind.width3 = ind.size + 1
		ind.rp_rs = ind.size + 2
		ind.b3   = ind.size + 3
		ind.flux  = ind.size + 4
		ind.size  += 5
		#DEFINE NAMES
		parname.insert(ind.midpt3,'Center')
		parname.insert(ind.width3,'Duration')
		parname.insert(ind.rp_rs,r'$R_p/R_s$')
		parname.insert(ind.b3,  'Impact Parameter')
		parname.insert(ind.flux3,r'System Flux ($\mu Jy$)')
		#DEFINE ECLIPSE MODEL
		myfuncs.append(mandel_geom)
		saveext.append('mg')
		functype[i] = 0	
	elif model[i] == 'mandelorbit':
		#DEFINE INDICES
		ind.e = ind.size
		ind.omega = ind.size + 1
		ind.i = ind.size + 2
		ind.rplanet   = ind.size + 3
		ind.rstar  = ind.size + 4
		ind.mstar = ind.size + 5
		ind.depth = ind.size + 6
		ind.flux4 = ind.size + 7
		ind.size  += 7
		#DEFINE NAMES
		parname.insert(ind.e,'Eccentricity')
		parname.insert(ind.omega,'Argument of Periapsis')
		parname.insert(ind.i,'Inclination')
		parname.insert(ind.rplanet,  r'Planet Radius ($R_J$)')
		parname.insert(ind.rstar, r'Star Radius ($R_\odot$)')
		parname.insert(ind.mstar, r'Star Mass ($M_\odot$)')
		parname.insert(ind.depth, 'Eclipse Depth')
		parname.insert(ind.flux4, r'System Flux ($\mu Jy$)')
		#DEFINE ECLIPSE MODEL
		myfuncs.append(mandelecl_orbit)
		saveext.append('mo')
		functype[i] = 0	
	elif model[i] == 'risingexp':
		#DEFINE INDICES
		ind.regoal  = ind.size
		ind.rem     = ind.size + 1
		ind.ret0    = ind.size + 2
		ind.size   += 3
		#DEFINE NAMES
		parname.insert(ind.regoal, 'Ramp, Ideal Flux')
		parname.insert(ind.rem,    'Ramp, Curvature')
		parname.insert(ind.ret0,   'Ramp, Phase Offset')
		#DEFINE RAMP MODEL
		myfuncs.append(risingexp)
		saveext.append('re')
		functype[i] = 1
   	elif model[i] == 'fallingexp':
		#DEFINE INDICES
		ind.fegoal  = ind.size
		ind.fem     = ind.size + 1
		ind.fet0    = ind.size + 2
		ind.size   += 3
		#DEFINE NAMES
		parname.insert(ind.fegoal, 'Ramp, Ideal Flux')
		parname.insert(ind.fem,    'Ramp, Curvature')
		parname.insert(ind.fet0,   'Ramp, Phase Offset')
		#DEFINE RAMP MODEL
		myfuncs.append(fallingexp)
		saveext.append('fe')
		functype[i] = 1
   	elif model[i] == 'quadramp':
		#DEFINE INDICES
		ind.qra     = ind.size
		ind.qrb     = ind.size + 1
		ind.qrc     = ind.size + 2
		ind.qrx0    = ind.size + 3
		ind.size   += 4
		#DEFINE NAMES
		parname.insert(ind.qra,    'Ramp, Quadratic Term')
		parname.insert(ind.qrb,    'Ramp, Linear Term')
		parname.insert(ind.qrc,    'Ramp, Constant Term')
		parname.insert(ind.qrx0,   'Ramp, Phase Offset')
   		#DEFINE RAMP MODEL
		myfuncs.append(quadramp)
		saveext.append('qd')
		functype[i] = 1
   	elif model[i] == 'linramp':
		#DEFINE INDICES
		ind.lina     = ind.size
		ind.linb     = ind.size + 1
		ind.linx0    = ind.size + 2
		ind.size    += 3
		#DEFINE NAMES
		parname.insert(ind.lina,    'Ramp, Linear Term')
		parname.insert(ind.linb,    'Ramp, Constant Term')
		parname.insert(ind.linb,    'Ramp, Phase Offset')
   		#DEFINE RAMP MODEL
		myfuncs.append(linramp)
		saveext.append('ln')
		functype[i] = 1
   	elif model[i] == 'logramp':
		#DEFINE INDICES
		ind.logt0    = ind.size
		ind.loga     = ind.size + 1
		ind.logb     = ind.size + 2
		ind.logc     = ind.size + 3
		ind.logd     = ind.size + 4
		ind.loge     = ind.size + 5
		ind.size    += 6
		#DEFINE NAMES
		parname.insert(ind.logt0,   'Ramp, Phase Offset')
		parname.insert(ind.loga,    'Ramp, Quartic Term')
		parname.insert(ind.logb,    'Ramp, Cubic Term')
		parname.insert(ind.logc,    'Ramp, Quadratic Term')
		parname.insert(ind.logd,    'Ramp, Linear Term')
		parname.insert(ind.loge,    'Ramp, Constant Term')
   		#DEFINE RAMP MODEL
		myfuncs.append(logramp)
		saveext.append('lg')
		functype[i] = 1
   	elif model[i] == 'log4qramp':
		#DEFINE INDICES
		ind.log4dt0    = ind.size
		ind.log4da     = ind.size + 1
		ind.log4db     = ind.size + 2
		ind.log4dc     = ind.size + 3
		ind.log4dd     = ind.size + 4
		ind.log4de     = ind.size + 5
		ind.log4df     = ind.size + 6
		ind.log4dg     = ind.size + 7
		ind.size    += 8
		#DEFINE NAMES
		parname.insert(ind.log4dt0,   'Ramp, Phase Offset')
		parname.insert(ind.log4da,    'Ramp, Quartic Log')
		parname.insert(ind.log4db,    'Ramp, Cubic Log')
		parname.insert(ind.log4dc,    'Ramp, Quadratic Log')
		parname.insert(ind.log4dd,    'Ramp, Linear Log')
		parname.insert(ind.log4de,    'Ramp, Quadratic Term')
		parname.insert(ind.log4df,    'Ramp, Linear Term')
		parname.insert(ind.log4dg,    'Ramp, Constant Term')
   		#DEFINE RAMP MODEL
		myfuncs.append(log4qramp)
		saveext.append('l4q')
		functype[i] = 1
   	elif model[i] == 'llramp':
		#DEFINE INDICES
		ind.llx0    = ind.size
		ind.lla     = ind.size + 1
		ind.llb     = ind.size + 2
		ind.llc     = ind.size + 3
		ind.size   += 4
		#DEFINE NAMES
		parname.insert(ind.llx0,   'Ramp, Phase Offset')
		parname.insert(ind.lla,    'Ramp, Log Term')
		parname.insert(ind.llb,    'Ramp, Linear Term')
		parname.insert(ind.llc,    'Ramp, Constant Term')
   		#DEFINE RAMP MODEL
		myfuncs.append(llramp)
		saveext.append('ll')
		functype[i] = 1
   	elif model[i] == 'lqramp':
		#DEFINE INDICES
		ind.llx0    = ind.size
		ind.lla     = ind.size + 1
		ind.llb     = ind.size + 2
		ind.llc     = ind.size + 3
		ind.lld     = ind.size + 4
		ind.size   += 5
		#DEFINE NAMES
		parname.insert(ind.llx0,   'Ramp, Phase Offset')
		parname.insert(ind.lla,    'Ramp, Log Term')
		parname.insert(ind.llb,    'Ramp, Quadratic Term')
		parname.insert(ind.llc,    'Ramp, Linear Term')
		parname.insert(ind.lld,    'Ramp, Constant Term')
   		#DEFINE RAMP MODEL
		myfuncs.append(lqramp)
		saveext.append('lq')
		functype[i] = 1
   	elif model[i] == 'sindecay':
		#DEFINE INDICES
		ind.sdx0    = ind.size
		ind.sda     = ind.size + 1
		ind.sdb     = ind.size + 2
		ind.sdc     = ind.size + 3
		ind.sdd     = ind.size + 4
		ind.size   += 5
		#DEFINE NAMES
		parname.insert(ind.sdx0,   'Decay, Phase Offset')
		parname.insert(ind.sda,    'Decay, Amplitude')
		parname.insert(ind.sdb,    'Decay, Exp Term')
		parname.insert(ind.sdc,    'Decay, Period')
		parname.insert(ind.sdd,    'Decay, Constant Term')
   		#DEFINE RAMP MODEL
		myfuncs.append(sindecay)
		saveext.append('sd')
		functype[i] = 1
   	elif model[i] == 'quadip':
		#DEFINE INDICES
		ind.qy2    = ind.size
		ind.qx2    = ind.size + 1
		ind.qxy    = ind.size + 2
		ind.qy     = ind.size + 3
		ind.qx     = ind.size + 4
		ind.qipc   = ind.size + 5
		ind.size  += 6
		#DEFINE NAMES
		parname.insert(ind.qy2,    'Intra-pixel, Quadratic Term in y')
		parname.insert(ind.qx2,    'Intra-pixel, Quadratic Term in x')
		parname.insert(ind.qxy,    'Intra-pixel, Cross Term')
		parname.insert(ind.qy,     'Intra-pixel, Linear Term in y')
		parname.insert(ind.qx,     'Intra-pixel, Linear Term in x')
		parname.insert(ind.qipc,   'Intra-pixel, Constant Term')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(quadip)
		saveext.append('qip')
		functype[i] = 2
   	elif model[i] == 'quadip4':
		#DEFINE INDICES
		ind.q0y2   = ind.size
		ind.q0x2   = ind.size + 1
		ind.q0xy   = ind.size + 2
		ind.q0y    = ind.size + 3
		ind.q0x    = ind.size + 4
		ind.q0c    = ind.size + 5
		ind.q1y2   = ind.size + 6
		ind.q1x2   = ind.size + 7
		ind.q1xy   = ind.size + 8
		ind.q1y    = ind.size + 9
		ind.q1x    = ind.size + 10
		ind.q1c    = ind.size + 11
		ind.q2y2   = ind.size + 12
		ind.q2x2   = ind.size + 13
		ind.q2xy   = ind.size + 14
		ind.q2y    = ind.size + 15
		ind.q2x    = ind.size + 16
		ind.q2c    = ind.size + 17
		ind.q3y2   = ind.size + 18
		ind.q3x2   = ind.size + 19
		ind.q3xy   = ind.size + 20
		ind.q3y    = ind.size + 21
		ind.q3x    = ind.size + 22
		ind.q3c    = ind.size + 23
		ind.size  += 24
		#DEFINE NAMES
		parname.insert(ind.q0y2,    'Intra-pixel, Q0, Quad Term in y')
		parname.insert(ind.q0x2,    'Intra-pixel, Q0, Quad Term in x')
		parname.insert(ind.q0xy,    'Intra-pixel, Q0, Cross Term')
		parname.insert(ind.q0y,     'Intra-pixel, Q0, Linear Term in y')
		parname.insert(ind.q0x,     'Intra-pixel, Q0, Linear Term in x')
		parname.insert(ind.q0c,     'Intra-pixel, Q0, Constant Term')
		parname.insert(ind.q1y2,    'Intra-pixel, Q1, Quad Term in y')
		parname.insert(ind.q1x2,    'Intra-pixel, Q1, Quad Term in x')
		parname.insert(ind.q1xy,    'Intra-pixel, Q1, Cross Term')
		parname.insert(ind.q1y,     'Intra-pixel, Q1, Linear Term in y')
		parname.insert(ind.q1x,     'Intra-pixel, Q1, Linear Term in x')
		parname.insert(ind.q1c,     'Intra-pixel, Q1, Constant Term')
		parname.insert(ind.q2y2,    'Intra-pixel, Q2, Quad Term in y')
		parname.insert(ind.q2x2,    'Intra-pixel, Q2, Quad Term in x')
		parname.insert(ind.q2xy,    'Intra-pixel, Q2, Cross Term')
		parname.insert(ind.q2y,     'Intra-pixel, Q2, Linear Term in y')
		parname.insert(ind.q2x,     'Intra-pixel, Q2, Linear Term in x')
		parname.insert(ind.q2c,     'Intra-pixel, Q2, Constant Term')
		parname.insert(ind.q3y2,    'Intra-pixel, Q3, Quad Term in y')
		parname.insert(ind.q3x2,    'Intra-pixel, Q3, Quad Term in x')
		parname.insert(ind.q3xy,    'Intra-pixel, Q3, Cross Term')
		parname.insert(ind.q3y,     'Intra-pixel, Q3, Linear Term in y')
		parname.insert(ind.q3x,     'Intra-pixel, Q3, Linear Term in x')
		parname.insert(ind.q3c,     'Intra-pixel, Q3, Constant Term')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(quadip4)
		saveext.append('qip4')
		functype[i] = 2
   	elif model[i] == 'cubicip':
		#DEFINE INDICES
		ind.cy3    = ind.size
		ind.cx3    = ind.size + 1
		ind.cy2x   = ind.size + 2
		ind.cyx2   = ind.size + 3
		ind.cy2    = ind.size + 4
		ind.cy2    = ind.size + 5
		ind.cx2    = ind.size + 6
		ind.cxy    = ind.size + 7
		ind.cy     = ind.size + 8
		ind.cx     = ind.size + 9
		ind.cipc   = ind.size + 10
		ind.size  += 11
		#DEFINE NAMES
		parname.insert(ind.cy3,    'Intra-pixel, Cubic Term in y')
		parname.insert(ind.cx3,    'Intra-pixel, Cubic Term in x')
		parname.insert(ind.cy2x,   'Intra-pixel, y^2*x Cross Term')
		parname.insert(ind.cyx2,   'Intra-pixel, y*x^2 Cross Term')
		parname.insert(ind.cy2,    'Intra-pixel, Quadratic Term in y')
		parname.insert(ind.cx2,    'Intra-pixel, Quadratic Term in x')
		parname.insert(ind.cxy,    'Intra-pixel, y*x Cross Term')
		parname.insert(ind.cy,     'Intra-pixel, Linear Term in y')
		parname.insert(ind.cx,     'Intra-pixel, Linear Term in x')
		parname.insert(ind.cipc,   'Intra-pixel, Constant Term')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(cubicip)
		saveext.append('cip')
		functype[i] = 2
   	elif model[i] == 'sexticip':
		#DEFINE INDICES
		ind.sy6    = ind.size
		ind.sx6    = ind.size + 1
		ind.sy5    = ind.size + 2
		ind.sx5    = ind.size + 3
		ind.sy4    = ind.size + 4
		ind.sx4    = ind.size + 5
		ind.sy3    = ind.size + 6
		ind.sx3    = ind.size + 7
		ind.sy2    = ind.size + 8
		ind.sx2    = ind.size + 9
		ind.sy     = ind.size + 10
		ind.sx     = ind.size + 11
		ind.sipc   = ind.size + 12
		ind.size  += 13
		#DEFINE NAMES
		parname.insert(ind.sy6,    'Intra-pixel, Sextic Term in y')
		parname.insert(ind.sx6,    'Intra-pixel, Sextic Term in x')
		parname.insert(ind.sy5,    'Intra-pixel, Quintic Term in y')
		parname.insert(ind.sx5,    'Intra-pixel, Quintic Term in x')
		parname.insert(ind.sy4,    'Intra-pixel, Quartic Term in y')
		parname.insert(ind.sx4,    'Intra-pixel, Quartic Term in x')
		parname.insert(ind.sy3,    'Intra-pixel, Cubic Term in y')
		parname.insert(ind.sx3,    'Intra-pixel, Cubic Term in x')
		parname.insert(ind.sy2,    'Intra-pixel, Quadratic Term in y')
		parname.insert(ind.sx2,    'Intra-pixel, Quadratic Term in x')
		parname.insert(ind.sy,     'Intra-pixel, Linear Term in y')
		parname.insert(ind.sx,     'Intra-pixel, Linear Term in x')
		parname.insert(ind.sipc,   'Intra-pixel, Constant Term')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(sexticip)
		saveext.append('6ip')
		functype[i] = 2
   	elif model[i] == 'sexticipc':
		#DEFINE INDICES
		ind.sy6    = ind.size
		ind.sx6    = ind.size + 1
		ind.sy5    = ind.size + 2
		ind.sx5    = ind.size + 3
		ind.sy4    = ind.size + 4
		ind.sx4    = ind.size + 5
		ind.sy3    = ind.size + 6
		ind.sx3    = ind.size + 7
		ind.sy2x   = ind.size + 8
		ind.sx2y   = ind.size + 9
		ind.sy2    = ind.size + 10
		ind.sx2    = ind.size + 11
		ind.sxy    = ind.size + 12
		ind.sy     = ind.size + 13
		ind.sx     = ind.size + 14
		ind.sipc   = ind.size + 15
		ind.size  += 16
		#DEFINE NAMES
		parname.insert(ind.sy6,    'Intra-pixel, Sextic Term in y')
		parname.insert(ind.sx6,    'Intra-pixel, Sextic Term in x')
		parname.insert(ind.sy5,    'Intra-pixel, Quintic Term in y')
		parname.insert(ind.sx5,    'Intra-pixel, Quintic Term in x')
		parname.insert(ind.sy4,    'Intra-pixel, Quartic Term in y')
		parname.insert(ind.sx4,    'Intra-pixel, Quartic Term in x')
		parname.insert(ind.sy3,    'Intra-pixel, Cubic Term in y')
		parname.insert(ind.sx3,    'Intra-pixel, Cubic Term in x')
		parname.insert(ind.sy2x,   'Intra-pixel, y^2*x Cross Term')
		parname.insert(ind.sx2y,   'Intra-pixel, y*x^2 Cross Term')
		parname.insert(ind.sy2,    'Intra-pixel, Quadratic Term in y')
		parname.insert(ind.sx2,    'Intra-pixel, Quadratic Term in x')
		parname.insert(ind.sxy,    'Intra-pixel, y*x Cross Term')
		parname.insert(ind.sy,     'Intra-pixel, Linear Term in y')
		parname.insert(ind.sx,     'Intra-pixel, Linear Term in x')
		parname.insert(ind.sipc,   'Intra-pixel, Constant Term')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(sexticipc)
		saveext.append('6ipc')
		functype[i] = 2
   	elif model[i] == 'medianip':
		#DEFINE INDICES
		ind.rad    = ind.size
		ind.size  += 1
		#DEFINE NAMES
		parname.insert(ind.rad,    'Intra-pixel, Radius')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(medianip)
		saveext.append('mip')
		functype[i] = 6
   	elif model[i] == 'nnint':
		#DEFINE INDICES
		ind.nnip   = ind.size
		ind.size  += 1
		#DEFINE NAMES
		parname.insert(ind.nnip,   'Interpolation, min # pts')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(nnint)
		saveext.append('nni')
		functype[i] = 6
   	elif model[i] == 'bilinint':
		#DEFINE INDICES
		ind.blip   = ind.size
		ind.size  += 1
		#DEFINE NAMES
		parname.insert(ind.blip,   'Interpolation, min # pts')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(bilinint)
		saveext.append('bli')
		functype[i] = 6
   	elif model[i] == 'ipspline':
		#DEFINE INDICES
		ind.ipsx0   = ind.size
		ind.ipsy0   = ind.size + 1
		ind.size  += 2
		#DEFINE NAMES
		parname.insert(ind.ipsx0,    'Intra-pixel, Knot x0')
		parname.insert(ind.ipsy0,    'Intra-pixel, Knot y0')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(ipspline)
		saveext.append('ipspl')
		functype[i] = 2
   	elif model[i] == 'posflux':
		#DEFINE INDICES
		ind.p0    = ind.size
		ind.p1    = ind.size + 1
		ind.p2    = ind.size + 2
		ind.p3    = ind.size + 3
		ind.p4    = ind.size + 4
		ind.p5    = ind.size + 5
		ind.p6    = ind.size + 6
		ind.p7    = ind.size + 7
		ind.p8    = ind.size + 8
		ind.size  += 9
		#DEFINE NAMES
		parname.insert(ind.p0,    'Position 0')
		parname.insert(ind.p1,    'Position 1')
		parname.insert(ind.p2,    'Position 2')
		parname.insert(ind.p3,    'Position 3')
		parname.insert(ind.p4,    'Position 4')
		parname.insert(ind.p5,    'Position 5')
		parname.insert(ind.p6,    'Position 6')
		parname.insert(ind.p7,    'Position 7')
		parname.insert(ind.p8,    'Position 8')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(posflux)
		saveext.append('pf')
		functype[i] = 3
   	elif model[i] == 'vissen':
		#DEFINE INDICES
		ind.vsx0    = ind.size
		ind.vsa    = ind.size + 1
		ind.vsb    = ind.size + 2
		ind.vsc    = ind.size + 3
		ind.size  += 4
		#DEFINE NAMES
		parname.insert(ind.vsx0,   'Visit Sensitivity, Phase Offset')
		parname.insert(ind.vsa,    'Visit Sensitivity, Log Term')
		parname.insert(ind.vsb,    'Visit Sensitivity, Linear Term')
		parname.insert(ind.vsc,    'Visit Sensitivity, Constant Term')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(vissen)
		saveext.append('vs')
		functype[i] = 4
   	elif model[i] == 'vsspline':
		#DEFINE INDICES
		ind.vss0    = ind.size
		ind.vss1    = ind.size + 1
		ind.vss2    = ind.size + 2
		ind.vss3    = ind.size + 3
		ind.vss4    = ind.size + 4
		ind.vss5    = ind.size + 5
		ind.vss6    = ind.size + 6
		ind.vss7    = ind.size + 7
		ind.vss8    = ind.size + 8
		ind.vss9    = ind.size + 9
		ind.vss10   = ind.size + 10
		ind.vss11   = ind.size + 11
		ind.size  += 12
		#DEFINE NAMES
		parname.insert(ind.vss0,   'Visit Sensitivity, Knot 0')
		parname.insert(ind.vss1,   'Visit Sensitivity, Knot 1')
		parname.insert(ind.vss2,   'Visit Sensitivity, Knot 2')
		parname.insert(ind.vss3,   'Visit Sensitivity, Knot 3')
		parname.insert(ind.vss4,   'Visit Sensitivity, Knot 4')
		parname.insert(ind.vss5,   'Visit Sensitivity, Knot 5')
		parname.insert(ind.vss6,   'Visit Sensitivity, Knot 6')
		parname.insert(ind.vss7,   'Visit Sensitivity, Knot 7')
		parname.insert(ind.vss8,   'Visit Sensitivity, Knot 8')
		parname.insert(ind.vss9,   'Visit Sensitivity, Knot 9')
		parname.insert(ind.vss10,  'Visit Sensitivity, Knot 10')
		parname.insert(ind.vss11,  'Visit Sensitivity, Knot 11')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(vsspline)
		saveext.append('vss')
		functype[i] = 4
   	elif model[i] == 'flatfield3':
		#DEFINE INDICES
		ind.ff30    = ind.size
		ind.ff31    = ind.size + 1
		ind.ff32    = ind.size + 2
		ind.ff33    = ind.size + 3
		ind.ff34    = ind.size + 4
		ind.ff35    = ind.size + 5
		ind.ff36    = ind.size + 6
		ind.ff37    = ind.size + 7
		ind.ff38    = ind.size + 8
		ind.flux4   = ind.size + 9
		ind.size   += 10
		#DEFINE NAMES
		parname.insert(ind.ff30,   'Flat Field, Pixel 0')
		parname.insert(ind.ff31,   'Flat Field, Pixel 1')
		parname.insert(ind.ff32,   'Flat Field, Pixel 2')
		parname.insert(ind.ff33,   'Flat Field, Pixel 3')
		parname.insert(ind.ff34,   'Flat Field, Pixel 4')
		parname.insert(ind.ff35,   'Flat Field, Pixel 5')
		parname.insert(ind.ff36,   'Flat Field, Pixel 6')
		parname.insert(ind.ff37,   'Flat Field, Pixel 7')
		parname.insert(ind.ff38,   'Flat Field, Pixel 8')
		parname.insert(ind.flux4, r'System Flux ($\mu Jy$)')
		#DEFINE INTRA-PIXEL MODEL
		myfuncs.append(flatfield)
		saveext.append('ff3')
		functype[i] = 5
   	else:
		print(str(model[i]) + " model not found.")
		myfuncs.append(-1)
		functype[i] = -1

   return myfuncs, functype, parname, ind, "".join(saveext)

def mandelecl(eclparams, x, etc = []):
   """
 NAME:
	MANDELECLIPSE

 PURPOSE:
	This function computes the secondary eclipse shape using equations provided by Mandel & Agol (2002)

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = MANDELECLIPSE([midpt, width, depth, x12, x34, out], x)

 INPUTS:
	midpt:	Center of eclipse
	width:	Eclipse duration from contacts 1 to 4
    	depth:	Eclipse depth
	x12:	Eclipse duration from contacts 1 to 2
	x34:	Eclipse duration from contacts 3 to 4
	out:	Flux offset from 0
	x:	Array of phase points

 OUTPUTS:
	This function returns the flux for each point in x.

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-05-08
			kevin218@knights.ucf.edu
   """
   #print("you made it in the function")
   midpt        = eclparams[0]
   width        = eclparams[1]
   depth        = eclparams[2]
   t12          = eclparams[3]
   t34          = eclparams[4]
   flux         = eclparams[5]
   
   #COMPUTE TIME OF CONTACT POINTS
   t1           = midpt - width / 2
   t2           = t1 + t12
   t4           = midpt + width / 2
   t3           = max(t4 - t34, t2)
   ieclipse     = np.where((x >= t2) & (x <= t3))
   iingress     = np.where((x >  t1) & (x <  t2))
   iegress      = np.where((x >  t3) & (x <  t4))

   p            = np.sqrt(abs(depth))*np.sign(depth)
   y            = np.ones(len(x))
   y[ieclipse]  = 1.0 - depth
   if p != 0:
   	  #Use Mandel & Agol (2002) for ingress of eclipse
      z           = -2 * p * (x[iingress] - t1) / t12 + 1 + p
      k0          = np.arccos((p ** 2 + z ** 2 - 1) / 2 / p / z)
      k1   	      = np.arccos((1 - p ** 2 + z ** 2) / 2 / z)
      y[iingress] = 1.0 - np.sign(depth) / np.pi * (p ** 2 * k0 + k1 	     \
                  - np.sqrt((4 * z ** 2 - (1 + z ** 2 - p ** 2) ** 2) / 4))
      #Use Mandel & Agol (2002) for egress of eclipse
      z    	      = 2 * p * (x[iegress] - t3) / t34 + 1 - p
      k0   	      = np.arccos((p ** 2 + z ** 2 - 1) / 2 / p / z)
      k1   	      = np.arccos((1 - p ** 2 + z ** 2) / 2 / z)
      y[iegress]  = 1.0 - np.sign(depth) / np.pi * (p ** 2 * k0 + k1 	     \
                  - np.sqrt((4 * z ** 2 - (1 + z ** 2 - p ** 2) ** 2) / 4))
   return y*flux

def mandel_geom(params, x, etc = []):
	"""
	NAME:
	MANDELGEOM

	PURPOSE:
	This function computes a transit shape using equations provided by Mandel & Agol (2002).

	CATEGORY:
	Astronomy.

	CALLING SEQUENCE:

	Result = mandelgeom([midpt, width, rp_rs, b, flux], x)

	INPUTS:
	midpt:	Center of eclipse
	width:	Eclipse duration from contacts 1 to 4
	rp_rs:	Planet-star radius ratio
	b:	Impact parameter
	flux:	Stellar flux
	x:	Array of phase points

	OUTPUTS:
	This function returns the flux for each point in x.

	PROCEDURE:

	EXAMPLE:



	MODIFICATION HISTORY:
	Written by:	Ryan A. Hardy, UCF  	2009-12-10
		hardy.r@gmail.com
	"""
	import orbit
	midpt, width, rp_rs, b, flux = params
	ingress = orbit.limbtime(b, width, 1, rp_rs)[0]
	trpars = np.array([midpt, width, rp_rs**2, ingress, ingress, flux])
	return mandelecl(trpars, x)
def mandelecl_orbit(params, x, etc = []):
	"""
	NAME:
	MANDELORBIT

	PURPOSE:
	This function computes the transit shape using equations provided by Mandel & Agol (2002)

	CATEGORY:
	Astronomy.

	CALLING SEQUENCE:

	Result = mandelgeom([midpt, width, rp_rs, b, flux], x)

	INPUTS:
	midpt:	Center of eclipse
	width:	Eclipse duration from contacts 1 to 4
	rp_rs:	Planet-star radius ratio
	b:	Impact parameter
	flux:	Stellar flux
	x:	Array of phase points

	OUTPUTS:
	This function returns the flux for each point in x.

	PROCEDURE:

	EXAMPLE:



	MODIFICATION HISTORY:
	Written by:	Ryan A. Hardy, UCF  	2009-12-10
		hardy.r@gmail.com
	"""
	G = 6.674e-11
	msun = 1.98892e30
	rsun = 6.955e8
	rjupiter = 71492000
	e, omega, i, period, rplanet, rstar, mstar, ecldepth, flux = params
	rplanet *= rjupiter
	rstar *= rsun
	mstar *= msun
	rp_rs = rplanet/rstar
	period*=86400
	a = (G*mstar*(period/(2*np.pi))**2)**(1/3.0)
	r = a*(1-e**2)/(1+e*np.cos(np.pi*(90-omega)/180.0))
	btr = r*np.cos(np.pi/180.0*i)/(rstar)
	trdur = orbit.duration(e, period/86400., omega, mstar/msun, rstar, rplanet, i*np.pi/180.)/(1440.*period/86400)
	#Fold data
	x = np.abs(x % 1)
	trlimb = orbit.limbtime(btr, trdur, 1, rp_rs)[0]
	ecldur, becl, ecllimb, midpt = orbit.scaled_eclipse(e, omega, trdur, rp_rs, btr)
	trpars = np.array([0, trdur, rp_rs**2, trlimb, trlimb, flux])
	eclpars = np.array([midpt, ecldur, ecldepth, ecllimb, ecllimb, flux])
	eclipse = mandelecl(eclpars, x)
	transit = mandelecl(trpars, x)
	return eclipse*transit/flux

def fallingexp(rampparams, x, etc = []):
#+
# NAME:
#	FALLINGEXP
#
# PURPOSE:
#	This function creates a model that fits a falling exponential with eclipse
#
# CATEGORY:
#	Astronomy.
#
# CALLING SEQUENCE:
#
#	Result = FALLINGEXP([midpt,width,depth,goal,m,x0], x)
#
# INPUTS:
#    	midpt:	Midpoint of eclipse
#	width:	Eclipse durations
#	depth:	Depth of eclipse
#	goal:	goal as x -> inf
#	m:	rise exp
#	x0:	time offset
#	x:	Array of time/phase points
#
# OUTPUTS:
#	This function returns an array of y values by combining an eclipse and a rising exponential
#
# SIDE EFFECTS:
#
# RESTRICTIONS:
#
# PROCEDURE:
#
# EXAMPLE:
#
#
#
# MODIFICATION HISTORY:
# 	Written by:	Kevin Stevenson, UCF  	2008-06-16
#			kevin218@knights.ucf.edu
#

   goal  = rampparams[0]
   m     = rampparams[1]
   x0    = rampparams[2]
   
   return ne.evaluate('goal * (1 + exp(-m * (x - x0)))')
   #return goal * (1 + np.exp(-m * (x - x0)))

def risingexp(rampparams, x, etc = []):
#+
# NAME:
#	RISINGEXP
#
# PURPOSE:
#	This function creates a model that fits a rising exponential with eclipse
#
# CATEGORY:
#	Astronomy.
#
# CALLING SEQUENCE:
#
#	Result = RISINGEXP([midpt,width,depth,x12,x34,goal,m,x0], x)
#
# INPUTS:
#    	midpt:	Midpoint of eclipse
#	width:	Eclipse durations
#	depth:	Depth of eclipse
#	x12:	Ingress time
#	x34:	Egress time
#	goal:	goal as x -> inf
#	m:	rise exp
#	x0:	time offset
#	x:	Array of time/phase points
#
# OUTPUTS:
#	This function returns an array of y values by combining an eclipse and a rising exponential
#
# SIDE EFFECTS:
#
# RESTRICTIONS:
#
# PROCEDURE:
#
# EXAMPLE:
#
#
#
# MODIFICATION HISTORY:
# 	Written by:	Kevin Stevenson, UCF  	2008-06-24
#			kevin218@knights.ucf.edu
#
   
   goal  = rampparams[0]
   m     = rampparams[1]
   x0    = rampparams[2]
   
   return ne.evaluate('goal*(1 - exp(-m*(x - x0)))')
   #return goal * (1 - np.exp(-m * (x - x0)))

def quadramp(rampparams, x, etc = []):
   """
 NAME:
	QUADRAMP

 PURPOSE:
	This function creates a model that fits a quadratically ramped eclipse

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = QUADRAMP([midpt,width,depth,a,b,c],x)

 INPUTS:
    	midpt:	Midpoint of eclipse
	width:	Eclipse durations
	depth:	Depth of eclipse
	a:	x^2 constant
	b:	x constant
	c:	x=0 offset
	x0: time/phase offset (constant)
	x:	Array of time/phase points

 OUTPUTS:
	This function returns an array of y values by combining an eclipse and a quadratic

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-06-22
			kevin218@knights.ucf.edu
   """

   a     = rampparams[0]
   b     = rampparams[1]
   c     = rampparams[2]
   x0    = rampparams[3]

   return ne.evaluate('a*(x-x0)**2 + b*(x-x0) + c')

def linramp(rampparams, x, etc = []):
   """
 NAME:
	LINRAMP

 PURPOSE:
	This function creates a model that fits a linear ramped eclipse

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = LINRAMP([a,b],x)

 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-07-07
			kevin218@knights.ucf.edu
   """

   a     = rampparams[0]
   b     = rampparams[1]
   x0    = rampparams[2]

   return ne.evaluate('a*(x-x0) + b')

def logramp(rampparams, x, etc = []):
   """
 NAME:
	LOGRAMP

 PURPOSE:
	This function creates a model that fits a natural log + linear ramped eclipse

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = LOGRAMP([midpt,width,depth,x12,x34,x0,b,c],x)

 INPUTS:
    
	x0:	time offset
	b:	x constant
	c:	x=0 offset
	x:	Array of time/phase points

 OUTPUTS:
	This function returns an array of y values by combining an eclipse and the ramp model

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-06-26
			kevin218@knights.ucf.edu
   """

   x0    = rampparams[0]
   a     = rampparams[1]
   b     = rampparams[2]
   c     = rampparams[3]
   d     = rampparams[4]
   e     = rampparams[5]
   
   return ne.evaluate('a*log(x-x0)**4 + b*log(x-x0)**3 + c*log(x-x0)**2 + d*log(x-x0) + e')

def llramp(rampparams, x, etc = []):
   """
 NAME:
	LLRAMP

 PURPOSE:
	This function creates a model that fits a log + linear ramped eclipse

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = LLRAMP([a,b,c],x)

 INPUTS:
	a:	log(x) constant
	b:	x constant
	c:	x=0 offset
	x:	Array of time/phase points

 OUTPUTS:
	This function returns an array of y values...

 EXAMPLE:


 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-08-31
			kevin218@knights.ucf.edu
   """

   x0    = rampparams[0]
   a     = rampparams[1]
   b     = rampparams[2]
   c     = rampparams[3]

   #Account for case where x < x0
   #Cannot have log of a negative number
   xnew = np.copy(x)
   xnew[np.where(x0 > x)] = x0 + 1e-15

   return ne.evaluate('a*log(xnew-x0) + b*(xnew-x0) + c')
   #return a*np.log(xnew-x0) + b*(xnew-x0) + c

def lqramp(rampparams, x, etc = []):
   """
 NAME:
	LQRAMP

 PURPOSE:
	This function creates a model that fits a log + quadratic ramped eclipse

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = LQRAMP([a,b,c,d],x)

 INPUTS:
    x0: phase offset
	a:	log(x) term
	b:	quadratic term
	c:	linear term
    d:  constant term
	x:	Array of time/phase points

 OUTPUTS:
	This function returns an array of y values...

 EXAMPLE:


 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2009-11-28
			kevin218@knights.ucf.edu
   """

   x0    = rampparams[0]
   a     = rampparams[1]
   b     = rampparams[2]
   c     = rampparams[3]
   d     = rampparams[4]

   #Account for case where x < x0
   #Cannot have log of a negative number
   xnew = np.copy(x)
   xnew[np.where(x0 > x)] = x0 + 1e-15

   return ne.evaluate('a*log(xnew-x0) + b*(xnew-x0)**2 + c*(xnew-x0) + d')
   #return a*np.log(xnew-x0) + b*(xnew-x0) + c


def log4qramp(rampparams, x, etc = []):
   """
 NAME:
	LOG4QRAMP

 PURPOSE:
	This function creates a model that fits a quartic-log + quadratic ramped eclipse

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = LQRAMP([a,b,c,d],x)

 INPUTS:
    x0: phase offset
	a:	log(x)^4 term
	b:	log(x)^3 term
	c:	log(x)^2 term
	d:	log(x) term
	e:	quadratic term
	f:	linear term
    g:  constant term
	x:	Array of time/phase points

 OUTPUTS:
	This function returns an array of y values...

 EXAMPLE:


 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2009-11-28
			kevin218@knights.ucf.edu
   """

   x0    = rampparams[0]
   a     = rampparams[1]
   b     = rampparams[2]
   c     = rampparams[3]
   d     = rampparams[4]
   e     = rampparams[5]
   f     = rampparams[6]
   g     = rampparams[7]

   #Account for case where x < x0
   #Cannot have log of a negative number
   xnew = np.copy(x)
   xnew[np.where(x0 > x)] = x0 + 1e-15

   return ne.evaluate('a*log(xnew-x0)**4 + b*log(xnew-x0)**3 + c*log(xnew-x0)**2 + d*log(xnew-x0) + e*(xnew-x0)**2 + f*(xnew-x0) + g')
   #return a*np.log(xnew-x0) + b*(xnew-x0) + c

def quadip(ipparams, position, etc = []):
   """
 NAME:
	QUADIP

 PURPOSE:
	This function fits the intra-pixel sensitivity effect using a 2-D quadratic.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = QUADIP(ipparams, position)

 INPUTS:
	midpt:	
	

 OUTPUTS:
	
 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-07-05
			kevin218@knights.ucf.edu
   """

   a       = ipparams[0]
   b       = ipparams[1]
   c       = ipparams[2]
   d       = ipparams[3]
   e       = ipparams[4]
   f       = ipparams[5]
   y, x, q = position

   return ne.evaluate('a*y**2 + b*x**2 + c*y*x + d*y + e*x + f')

def quadip4(ipparams, position, etc = []):
   """
 NAME:
	QUADIP4

 PURPOSE:
	This function fits the intra-pixel sensitivity effect using a 2-D quadratic in each pixel quadrant.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = QUADIP(ipparams, position)

 INPUTS:
	midpt:	
	

 OUTPUTS:
	

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-08-18
			kevin218@knights.ucf.edu
   """

   a0      = ipparams[0]
   b0      = ipparams[1]
   c0      = ipparams[2]
   d0      = ipparams[3]
   e0      = ipparams[4]
   f0      = ipparams[5]
   a1      = ipparams[6]
   b1      = ipparams[7]
   c1      = ipparams[8]
   d1      = ipparams[9]
   e1      = ipparams[10]
   f1      = ipparams[11]
   a2      = ipparams[12]
   b2      = ipparams[13]
   c2      = ipparams[14]
   d2      = ipparams[15]
   e2      = ipparams[16]
   f2      = ipparams[17]
   a3      = ipparams[18]
   b3      = ipparams[19]
   c3      = ipparams[10]
   d3      = ipparams[21]
   e3      = ipparams[22]
   f3      = ipparams[23]
   y, x, q = position

   y0      = y[np.where(q == 0)]
   x0      = x[np.where(q == 0)]
   y1      = y[np.where(q == 1)]
   x1      = x[np.where(q == 1)]
   y2      = y[np.where(q == 2)]
   x2      = x[np.where(q == 2)]
   y3      = y[np.where(q == 3)]
   x3      = x[np.where(q == 3)]

   output  = np.zeros(y.size)

   output[np.where(q == 0)] = ne.evaluate('a0*y0**2 + b0*x0**2 + c0*y0*x0 + d0*y0 + e0*x0 + f0')
   output[np.where(q == 1)] = ne.evaluate('a1*y1**2 + b1*x1**2 + c1*y1*x1 + d1*y1 + e1*x1 + f1')
   output[np.where(q == 2)] = ne.evaluate('a2*y2**2 + b2*x2**2 + c2*y2*x2 + d2*y2 + e2*x2 + f2')
   output[np.where(q == 3)] = ne.evaluate('a3*y3**2 + b3*x3**2 + c3*y3*x3 + d3*y3 + e3*x3 + f3')

   return (output)

def cubicip(ipparams, position, etc = []):
   """
 NAME:
	CUBICIP

 PURPOSE:
	This function fits the intra-pixel sensitivity effect using a 2-D cubic .

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = CUBICIP(ipparams, position)

 INPUTS:
	midpt:	
	

 OUTPUTS:
	

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2008-07-08
			kevin218@knights.ucf.edu
   """

   a       = ipparams[0]
   b       = ipparams[1]
   c       = ipparams[2]
   d       = ipparams[3]
   e       = ipparams[4]
   f       = ipparams[5]
   g       = ipparams[6]
   h       = ipparams[7]
   i       = ipparams[8]
   j       = ipparams[9]
   y, x, q = position

   return ne.evaluate('a*y**3 + b*x**3 + c*y**2*x + d*y*x**2 + e*y**2 + f*x**2 + g*y*x + h*y + i*x + j')

def sexticip(ipparams, position, etc = []):
   """
 NAME:
	SEXTIP

 PURPOSE:
	This function fits the intra-pixel sensitivity effect using a 2-D 6th-order polynomial.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = SEXTICIP(ipparams, position)

 INPUTS:
	midpt:	
	

 OUTPUTS:
	

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-01-22
			kevin218@knights.ucf.edu
   """

   y6      = ipparams[0]
   x6      = ipparams[1]
   y5      = ipparams[2]
   x5      = ipparams[3]
   y4      = ipparams[4]
   x4      = ipparams[5]
   y3      = ipparams[6]
   x3      = ipparams[7]
   y2      = ipparams[8]
   x2      = ipparams[9]
   y1      = ipparams[10]
   x1      = ipparams[11]
   c       = ipparams[12]
   y, x, q = position

   return ne.evaluate('y6*y**6 + x6*x**6 + y5*y**5 + x5*x**5 + y4*y**4 + x4*x**4 + \
                       y3*y**3 + x3*x**3 + y2*y**2 + x2*x**2 + y1*y + x1*x + c')

def sexticipc(ipparams, position, etc = []):
   """
 NAME:
	SEXTIPC

 PURPOSE:
	This function fits the intra-pixel sensitivity effect using a 2-D 6th-order polynomial,
    with cross terms.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = SEXTICIPC(ipparams, position)

 INPUTS:
	midpt:	
	

 OUTPUTS:
	

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-02-01
			kevin218@knights.ucf.edu
   """

   y6      = ipparams[0]
   x6      = ipparams[1]
   y5      = ipparams[2]
   x5      = ipparams[3]
   y4      = ipparams[4]
   x4      = ipparams[5]
   y3      = ipparams[6]
   x3      = ipparams[7]
   y2x     = ipparams[8]
   x2y     = ipparams[9]
   y2      = ipparams[10]
   x2      = ipparams[11]
   xy      = ipparams[12]
   y1      = ipparams[13]
   x1      = ipparams[14]
   c       = ipparams[15]
   y, x, q = position

   return ne.evaluate('y6*y**6 + x6*x**6 + y5*y**5 + x5*x**5 + y4*y**4 + x4*x**4 + y3*y**3 + x3*x**3 + \
                       y2x*y**2*x + x2y*x**2*y + y2*y**2 + x2*x**2 + xy*x*y + y1*y + x1*x + c')

def medianip(ipparams, posflux, etc = [], retbinflux = False):
    """
 NAME:
	MEDIANIP

 PURPOSE:
	This function fits the intra-pixel sensitivity effect using the median 
    within a given radius of the current position.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = MEDINAIP(ipparams, position)

 INPUTS:
	midpt:	
	

 OUTPUTS:
	

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-06-06
			kevin218@knights.ucf.edu
    """
    
    r          = ipparams[0]
    y, x, flux, whereltrad = posflux[0:4]
    
    fluxcorr   = np.ones(flux.size)
    print(len(flux), len(etc))
    fluxip     = flux / etc
    for i in range(flux.size):
        #if i%10000 == 0:
            #print(i)
        fluxcorr[i] = np.median(fluxip[whereltrad[i]])
    
    if retbinflux == False:
        return fluxcorr
    else:
        return [fluxcorr, -1]

def nnint(ipparams, posflux, etc = [], retbinflux = False):
    """
 NAME:
	NNINT

 PURPOSE:
	This function fits the intra-pixel sensitivity effect using the mean 
    within a given binned position (nearest-neighbor interpolation).

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = NNINT(ipparams, position)

 INPUTS:
	midpt:	
	

 OUTPUTS:
	

 SIDE EFFECTS:

 RESTRICTIONS:

 PROCEDURE:

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-06-07
			kevin218@knights.ucf.edu
    """
    #ystep, xstep = ipparams
    y, x, flux, wherebinflux = posflux[0:4]
    szwbf      = len(wherebinflux)
    fluxcorr   = np.ones(flux.size)
    fluxip     = flux / etc
    binflux    = np.zeros(szwbf)
    for i in range(szwbf):
        if len(wherebinflux[i][0]) > 0:
            fluxcorr[wherebinflux[i]] = binflux[i] = np.mean(fluxip[wherebinflux[i]])
    
    #Return fit with or without binned flux
    if retbinflux == False:
        return fluxcorr
    else:
        return [fluxcorr, binflux]
    ''' 
    ystep, xstep = ipparams
    y, x, flux, wherebinflux = position
    ####ASSIGN BINS FOR 2D BINNING
    ymin     = np.floor(10*y.min())/10
    ymax     = np.ceil (10*y.max())/10
    xmin     = np.floor(10*x.min())/10
    xmax     = np.ceil (10*x.max())/10
    ysize    = (ymax-ymin)/ystep + 1                    #Number of bins in y dimension
    xsize    = (xmax-xmin)/xstep + 1                    #Number of bins in x dimension
    binx, biny = np.meshgrid(np.linspace(xmin,xmax,xsize),np.linspace(ymin,ymax,ysize))
    ###
    meanbinflux = np.ones(len(wherebinflux))
    weight      = np.ones(len(wherebinflux))
    fluxcorr    = np.ones(flux.size)
    fluxip      = flux / etc
    for i in range(len(wherebinflux)):
        if len(wherebinflux[i][0]) > 0:
            meanbinflux[i] = np.mean(fluxip[wherebinflux[i]])
            weight[i]      = 1
    
    tck = spi.bisplrep(binx.flatten(), biny.flatten(), meanbinflux, w=weight, kx=1, ky=1, 
                       task=-1, ty=np.linspace(ymin,ymax,7), tx=np.linspace(xmin,xmax,7))
    #print(tck)
    output = np.ones(y.size)
    for i in range(y.size):
        output[i] = spi.bisplev(x[i], y[i], tck, dx=0, dy=0)
    
    return output
    '''

def bilinint(ipparams, posflux, etc = [], retbinflux = False):
    """
    Uses bilinear interpolation to fit binned flux vs position.
    
	This function fits the intra-pixel sensitivity effect using 
    bilinearly interpolated flux of mean binned flux values.

    Parameters
    ----------
	ipparams :  tuple
                unused
    y :         1D array, size = # of measurements
                Pixel position along y
    x :         1D array, size = # of measurements
                Pixel position along x
    flux :      1D array, size = # of measurements
                Observed flux at each position
    wherebinflux :  1D array, size = # of bins
                    Measurement number assigned to each bin
    gridpt :    1D array, size = # of measurements
                Bin number in which each measurement is assigned
    dy1 :       1D array, size = # of measurements
                (y - y1)/(y2 - y1)
    dy2 :       1D array, size = # of measurements
                (y2 - y)/(y2 - y1)
    dx1 :       1D array, size = # of measurements
                (x - x1)/(x2 - x1)
    dx2 :       1D array, size = # of measurements
                (x2 - x)/(x2 - x1)
    ysize :     int
                Number of bins along y direction
    xsize :     int
                Number of bins along x direction
    
    Returns
    -------
    output :    1D array, size = # of measurements
                Intra-pixel-corrected flux multiplier

    Optional
    --------
    binflux :   1D array, size = # of bins
                Binned Flux values

    Examples
    --------
    None

    Revisions
    ---------
    2010-06-11  Kevin Stevenson, UCF
			    kevin218@knights.ucf.edu
                Original version
    """
    
    y, x, flux, wherebinflux, [gridpt, dy1, dy2, dx1, dx2], [ysize, xsize] = posflux
    gridpt     = gridpt.astype(int)
    szwbf      = len(wherebinflux)
    binflux    = np.zeros(szwbf)
    ipflux     = flux / etc
    for i in range(szwbf):
        if len(wherebinflux[i][0]) > 0:
            binflux[i] = np.mean(ipflux[wherebinflux[i]])

    output = binflux[gridpt      ]*dy2*dx2 + binflux[gridpt      +1]*dy2*dx1 + \
             binflux[gridpt+xsize]*dy1*dx2 + binflux[gridpt+xsize+1]*dy1*dx1
    #Return fit with or without binned flux
    if retbinflux == False:
        return output
    else:
        return [output, binflux]

def ipspline(ipparams, position, etc = []):
    """
 NAME:
	IPSPLINE

 PURPOSE:
	This function creates a bivariate B-spline that fits the intra-pixel effect.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = IPSPLINE([k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11], position, etc)

 INPUTS:
    k#:     Knot coefficient
	x,y:    Array of x,y pixel positions and quadrant locations
    etx:    Knot locations

 OUTPUTS:
	This function returns an array of y values...

 REFERENCES:
    See SI from Harrinton et al. (2007)

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-06-08
			kevin218@knights.ucf.edu
    """
    y, x, q = position
    yknots, xknots = etc
    
    tck = spi.bisplrep(xknots.flatten(), yknots.flatten(), ipparams, kx=3, ky=3)
    #print(tck)
    #tck = [yknots, xknots, ipparams, 3, 3]
    #func = spi.interp2d(xknots, yknots, ipparams, kind='cubic'
    output = np.ones(y.size)
    for i in range(y.size):
        output[i] = spi.bisplev(x[i], y[i], tck, dx=0, dy=0)
    
    return output

def sindecay(rampparams, x, etc = []):
   """
 NAME:
	LLRAMP

 PURPOSE:
	This function creates a model that fits a sinusoidal decay

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = SINDECAY([a,b,c,d],x)

 INPUTS:
    x0: phase/time offset
	a:	amplitude
	b:	exponential constant
	c:	period
    d:  vertical offset
	x:	Array of time/phase points

 OUTPUTS:
	This function returns an array of y values...

 EXAMPLE:


 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2009-07-26
			kevin218@knights.ucf.edu
   """

   x0    = rampparams[0]
   a     = rampparams[1]
   b     = rampparams[2]
   c     = rampparams[3]
   d     = rampparams[4]

   return ne.evaluate('a*exp(b*x)*cos(2*pi*(x-x0)/c) + d')

def posflux(posparams, x, etc = []):
    """
 NAME:
	POSFLUX

 PURPOSE:
	This function creates a model that fits the median flux at each mosition

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = POSFLUX([p0,p1,p2,p3,p4,p5,p6,p7,p8],x)

 INPUTS:
    p#: position #
	x:	Array of position points

 OUTPUTS:
	This function returns an array of y values...

 EXAMPLE:


 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-01-20
			kevin218@knights.ucf.edu
    """
     
    normfactors = np.ones(x.size)
    for i in range(np.size(posparams)):
        normfactors[np.where(x == i)] = posparams[i]
    
    return normfactors

def vissen(visparams, (x, knots), etc = []):
    """
 NAME:
	VISSEN

 PURPOSE:
	This function creates a model that fits the visit sensitivity.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = VISSEN([p0,p1,p2,p3,p4,p5,p6,p7,p8],x)

 INPUTS:
    p#:    position #
	x:	   Array of frame numbers in current visit
    knots: Not required for this function

 OUTPUTS:
	This function returns an array of y values...

 REFERENCES:
    See SI from Harrinton et al. (2007)

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-01-20
			kevin218@knights.ucf.edu
    """
    
    x0    = visparams[0]
    a     = visparams[1]
    b     = visparams[2]
    c     = visparams[3]
    
    #Account for case where x < x0
    #Cannot have log of a negative number
    xnew = np.copy(x)
    xnew[np.where(x0 > x)] = x0 + 1e-15
    
    return ne.evaluate('a*log(xnew-x0) + b*(xnew-x0) + c')

def vsspline(visparams, (x, knots), etc = []):
    """
 NAME:
	VSSPLINE

 PURPOSE:
	This function creates a cubic spline that fits the visit sensitivity.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = VSSPLINE([k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11],(x, knots))

 INPUTS:
    k#:    knot coefficient
	x:	   Array of frame numbers in current visit
    knots: knot locations

 OUTPUTS:
	This function returns an array of y values...

 REFERENCES:
    See SI from Harrinton et al. (2007)

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-04-03
			kevin218@knights.ucf.edu
    """
    '''
    k0    = visparams[0]
    k1    = visparams[1]
    k2    = visparams[2]
    k3    = visparams[3]
    k4    = visparams[4]
    k5    = visparams[5]
    k6    = visparams[6]
    k7    = visparams[7]
    k8    = visparams[8]
    k9    = visparams[9]
    k10   = visparams[10]
    k11   = visparams[11]
    '''
    tck = spi.splrep(knots, visparams, k=3)
    #tck = (knots, visparams, 3)
    
    return spi.splev(x, tck, der=0)


def flatfield(ffparams, ap, etc = []):
    """
 NAME:
	VSSPLINE

 PURPOSE:
	This function perform flat field correction on a n-by-n pixel aperture.

 CATEGORY:
	Astronomy.

 CALLING SEQUENCE:

	Result = FLATFIELD([ff0,ff1,ff3,ff3,ff4,ff5,ff6,ff7,ff8],(flux, ap))

 INPUTS:
    ff#:    Flux modifier, starting with upper left pixel, reads normally (shape = n*n or n,n)
	flux:   Array of flux values (shape = #images)
    ap:     3D array of flux values within the aperture (shape = n,n,#images)

 OUTPUTS:
	This function returns an array of new flux values.

 COMMENTS:
    Apertures for each image must be the same size as ffparams. (ap[0].size=ffparams.size)

 REFERENCES:
    See SI from Harrinton et al. (2007)

 EXAMPLE:



 MODIFICATION HISTORY:
 	Written by:	Kevin Stevenson, UCF  	2010-05-10
			kevin218@knights.ucf.edu
    """
    ff   = np.array(ffparams[:-1])
    flux = ffparams[-1]
    #if ff.ndim == 1:
    #    ff = np.reshape(ff, (ap.shape[1],ap.shape[2]))
    ap = np.reshape(ap, (-1, ff.size))

    #Test which method is faster
    #deltaflux = np.sum(np.sum((ffparams*ap - ap), axis=1), axis=1)
    #ne.evaluate('flux + deltaflux')
    return (np.sum((ff*ap - ap), axis=1) + flux)


