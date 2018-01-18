
import numpy as np

# properties for GJ436b;
def gj436b(info):

   # measured stellar parameters
   #info.ra = sexi2dec('11:42:11.0941') / 12e0 * pi # SIMBAD 
   #info.dec = sexi2dec('+26:42:23.652') / 180e0 * pi # SIMBAD
   info.parallax     = 97.73e-3      # mas, SIMBAD 
   info.parallaxerr  = 2.27e-3       # ", SIMBAD 
   info.lstar        = 0.019 * info.lsun # Wikipedia
   info.lstarerr     = 0.019 * info.lsun # FINDME:UNKNOWN ERROR
   info.rstar        = 0.464e0 * info.rsun # m, Rsun, Torres 2007
   info.rstarerr     = 0.010e0 * info.rsun # m, Rsun, Torres 2007
   info.metalstar    = -0.32e0 # [Fe/ H], C. Bean et al. 2006
   info.metalstarerr = 0.12e0 # [Fe/ H], C. Bean et al. 2006
   info.tstar        = 3684e0 # K, Torres 2007
   info.tstarerr     = 71e0 # K, Torres 2007
   info.logg         = 4.80e0 # alog10(cm s-2), C. Bean et al. 2006
   info.loggerr      = 0.10e0 # alog10(cm s-2), C. Bean et al. 2006
   
   # measured planet parameters
   info.rplan        = 0.438e0 * info.rjup # m, Bean et al. 2008
   info.rplanerr     = 0.040e0 * info.rjup  # m, Bean et al. 2008
   info.semimaj      = 0.02872e0 * info.au # Bean et al. 2008
   info.semimajerr   = 0.00006e0 * info.au # Bean et al. 2008
   info.incl         = 85.8e0 * np.pi / 180. # Bean et al. 2008
   info.inclerr      = 0.25e0 * np.pi / 180. # Bean et al. 2008
   info.ephtime      = 2454222.616e0   # days, Gillon et al. 2007
   info.ephtimeerr   = 0.001e0   # days, Gillon et al. 2007
   info.period       = 2.643904e0 # Bean et al. 2008
   info.perioderr    = 0.000005e0 # Bean et al. 2008
   info.transdur     = 0.96e0 / 24.e0   # FINDME: Reference this value
   info.transdurerr  = 0.96e0 / 24.e0   # FINDME: Need to find FINDME:OLDVALUE
   info.arat         = (info.rplan / info.rstar) ** 2   # Area Ratio FINDME: tentative 
   info.araterr      = (info.arat)*2*((info.rplanerr/info.rplan)**2 + (info.rstarerr/info.rstar)**2)**0.5
   
   return info

# orbit data for HAT-P-1b
def hatp1b(info):

   # measured stellar parameters
   info.ra = 0 #sexi2dec('22:57:46.825') / 12e0 * np.pi # 
   info.dec = 0 #sexi2dec('38:40:29.83') / 180e0 * np.pi # 
   # RA & Dec from SIMBAD which references Roeser & Bastian 1988
   
   info.parallax = 7.194e-3      # mas, FINDME:CITESOURCE
   info.parallaxerr = 1.139e-3       # ", FINDME:CITESOURCE
   info.lstar = 1.51e0 * info.lsun # W, FINDME:CITESOURCE
   #info.lstarerr = x.xxd * info.lsun # W, FINDME:CITESOURCE
   info.rstar = 1.115e0 * info.rsun # m, Rsun,  Johnson et al. 2008
   info.rstarerr = 0.050e0 * info.rsun # m, Rsun,  Johnson et al. 2008
   info.metalstar = 0.13e0 # [Fe/ H], Bakos et al. 2007
   info.metalstarerr = 0.02e0 # [Fe/ H], Bakos et al. 2007
   info.tstar = 5975e0 # K, Bakos et al. 2007
   info.tstarerr = 45e0 # K, Bakos et al. 2007
   info.logg = 4.45e0 # alog10(cm s-2), Bakos et al. 2007
   info.loggerr = 0.06e0 # alog10(cm s-2), Bakos et al. 2007
   
   # measured planet parameters
   info.rplan = 1.225e0 * info.rjup # m, Johnson et al. 2008
   info.rplanerr = 0.056e0 * info.rjup # m, Johnson et al. 2008
   info.semimaj = 0.0553e0 * info.au   # Johnson et al. 2008
   info.semimajerr = 0.0014e0 * info.au   # Johnson et al. 2008
   info.incl = 86.28e0 * np.pi / 180e0 # Johnson et al. 2008
   info.inclerr = 0.20e0 * np.pi / 180e0 # Johnson et al. 2008
   info.ephtime = 2454363.94656e0 # days, Johnson et al. 2008
   info.ephtimeerr = 0.00072e0       # days, Johnson et al. 2008
   info.period = 4.4652934e0     # days, Johnson et al. 2008
   info.perioderr = 0.000093e0      # days, Johnson et al. 2008
   info.transdur = 2.798e0         # hours, Johnson et al. 2008
   info.transdurerr = 0.019e0         # hours, Johnson et al. 2008
   info.arat = 0.11295 ** 2      # Johnson et al. 2008
   info.araterr = 0.00073 * 1.414  # Johnson et al. 2008
   
   return info


# set universal constants in info structure
def univ(info):
   import numpy as np

   # conversions
   # steradians per sq. arcsecond, from MDH 2.1, p. 25., sect. 3.7.1: 2.35045e-11
   info.srperas = 4.e0 * np.pi / ((360.e0 * 60.e0 * 60.e0) ** 2 / np.pi)
   # 1d-6  converts uJy to Jy, 1d-26 converts Jy to W m^-2 Hz^-1
   info.ujy2mks = 1e-6 * 1e-26
   
   # time
   info.mjdoff = 2400000.5e0
   info.j2kjd = 2451545.e0 # Julian date of J2000.0=1.5 Jan 2000 (see ESAA ch 27)
   
   # km per parsec
   info.mppc = 3.0856776e16   # m, AllenII, p.12
   
   # units
   info.lsun = 3.827e26        # W, Wikipedia
   info.rsun = 695508.e3       # m, +- 26 km, AllenII
   info.rjup = 71492.e3        # m, AllenII
   info.au = 1.4959787066e11 # m, AllenII
   info.stefboltz = 5.67040e-8      # J m^-2 s^-1 K^-4, Wikipedia
   # see also AllenII, p.11 (less accurate)
   info.c = 2.99792458e8    # m/s, speed of light
   info.h = 6.6260693e-34   # J s, Planck's constant, Wikipedia
   info.k = 1.3806503e-23   # J/K, Boltzmann's constant, Google
   
   return info


