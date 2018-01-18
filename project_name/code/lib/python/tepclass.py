# $Author: ccampo $
# $Date: 2011-10-21 14:19:01 -0400 (Fri, 21 Oct 2011) $
# $Revision: 573 $
# $HeadURL: file:///home/esp01/svn/code/libpython/trunk/tepclass.py $
# $Id: tepclass.py 573 2011-10-21 18:19:01Z ccampo $

import numpy as np
import sexa2dec as s2d

# constants
rsun   = 6.95508e8   # solar radius, m,    Wikipedia
msun   = 1.98892e30  # solar mass, kg,     Wikipedia
rjup   = 7.1492e7    # jupiter radius, m,  Wikipedia
mjup   = 1.8689e27   # jupiter mass, kg,   Wikipedia
rearth = 6378.1*1e3  # earth radius, m,    Wikipedia
au2m   = 1.49598e11  # 1 AU in meters,     Wikipedia
d2s    = 86400.      # 1 day in sec,       Wikipedia
yr2s   = 31556926.   # 1 year in sec,      Wikipedia
hr2s   = 3600.       # 1 hour in sec,      Wikipedia
mearth = 5.9722e24

# custom defined exception; used below.
class UnitError(Exception):
    def __init__(self, val):
        self.val = val
    def __str__(self):
        return repr(self.val)

# each parameter is an instance of this class
class param:
    # constructor provides 4 variables; value, uncertainty, unit, and reference
    def __init__(self, val, uncert, unit, ref):
        self.val     = val
        self.uncert  = uncert
        self.unit    = unit
        self.ref     = ref

# class of the tepfile
class tepfile:
    # generate variable names dynamically to the names in the tepfile
    def __init__(self, filename, conv=True):
        self.fname = filename  # filename of the tepfile

        # tepfile version number; needed to understand how to parse
        tepin  = open(filename, 'r')
        lines  = tepin.readlines()
        verstr = lines[1]
        verstr = verstr[17:20] # indicies of version num (can be 3 digits long)
        tepin.close()

        # get the version number
        if verstr.endswith(', '):
            verstr = verstr.strip(', ')
        else:
            verstr = verstr.strip(',')
        self.version = np.float(verstr) # add version to tepfile

        # array of values as strings
        valarr = np.loadtxt(filename, dtype=np.str)

        # parameter names
        parnames = valarr[:, 0]

        # values to keep as strings
        strvals = [
            'planetname',
            'startype',
            'ra',
            'dec',
            'announcedate',
            'submitter',
            'email',
            'date',
            'time',
            ]

        # load all parameters
        for i, parname in enumerate(parnames):
            parname = parname.lower()
            # versions < 4 do not have units or refs
            if self.version < 4:
                exec('self.{pname} = param("{val}", "{uncert}",' \
                     ' None, None)'.format(pname  = parname,
                                           val    = valarr[i, 1],
                                           uncert = valarr[i, 2]))
            elif self.version == 4:
                exec('self.{pname} = param("{val}", "{uncert}",' \
                     ' "{unit}", "{ref}")'.format(pname  = parname,
                                                  val    = valarr[i, 1],
                                                  uncert = valarr[i, 2],
                                                  unit   = valarr[i, 3],
                                                  ref    = valarr[i, 4]))
            # assign floats to decimal values
            if parname not in strvals:
                exec("self.{pname}.val    = " \
                     "np.float(self.{pname}.val)".format(pname = parname))
                exec("self.{pname}.uncert = " \
                     "np.float(self.{pname}.uncert)".format(pname = parname))

        # do the unit conversion to MKS
        if conv:
            self.convunits()

    ##################################
    # method to convert units to MKS #
    ##################################
    def convunits(self):
        """
        Converts units to MKS.
        """
        # convert ra and dec to radians
        self.ra.val  = s2d.sexa2dec(self.ra.val)  /  12.0 * np.pi
        self.dec.val = s2d.sexa2dec(self.dec.val) / 180.0 * np.pi
        # update units
        self.ra.unit  = 'rad'
        self.dec.unit = 'rad'

        if self.version >= 4:
            # check temperature
            if self.ts.val != -1:
                # celcius to kelvin
                if self.ts.unit.lower() == 'c':
                    self.ts.val     += 273.
                    if self.ts.uncert != -1:
                        self.ts.uncert  += 273.
                    self.ts.unit     = 'K'

                # already kelvin
                elif self.ts.unit.lower() == 'k':
                    pass
                else:
                    raise UnitError("unrecognized units for self.ts!")

            # check stellar radius
            if self.rs.val != -1:
                #  solar radii to meters
                if self.rs.unit.lower() == 'rsun':
                    self.rs.val    *= rsun
                    if self.rs.uncert != -1:
                        self.rs.uncert *= rsun
                    self.rs.unit    = 'm'

                #  earth radii to meters
                elif self.rs.unit.lower() == 'rearth':
                    self.rs.val    *= rearth
                    if self.rs.uncert != -1:
                        self.rs.uncert *= rearth
                    self.rs.unit    = 'm'

                # km to m
                elif self.rs.unit.lower() == 'km':
                    self.rs.val    *= 1e3
                    if self.rs.uncert != -1:
                        self.rs.uncert *= 1e3
                    self.rs.unit    = 'm'

                # already m
                elif self.rs.unit.lower() == 'm':
                    pass
                else:
                    raise UnitError("unrecognized units for self.rs!")

            # check log(g)
            if self.loggstar.val != -1:
                # cgs to mks
                if self.loggstar.unit.lower() == 'cgs':
                    g    = 10**(self.loggstar.val)
                    g   *= 0.01   # cm to meters
                    # do uncert
                    if self.loggstar.uncert != -1:
                        guc    = 10**(self.loggstar.uncert)
                        guc   *= 0.01
                    else:
                        guc = -1
                # already mks
                elif self.loggstar.unit.lower() == 'mks':
                    g   = 10**(self.loggstar.val)
                    guc = 10**(self.loggstar.uncert)
                else:
                    raise UnitError("unrecognized units for self.loggstar!")

                # g is new value in MKS; logg stays the same
                self.g = param(g, guc, 'mks', None)
            else:
                self.g = param(-1, -1, None, None)


            # check stellar mass
            if self.ms.val != -1:
                # solar mass to kg
                if self.ms.unit.lower() == 'msun':
                    self.ms.val    *= msun
                    if self.ms.uncert != -1:
                        self.ms.uncert *= msun
                    self.ms.unit    = 'kg'

                # already kg
                elif self.ms.unit.lower() == 'kg':
                    pass
                else:
                    raise UnitError("unrecognized units for self.ms!")

            # FINDME: check epoch

            # check pmRA
            if self.pmra.val != -1:
                # mas to arcsec
                if self.pmra.unit.lower() == 'marcsec/year' or \
                       self.pmra.unit.lower() == 'mas/year' or \
                       self.pmra.unit.lower() == 'mas/yr':
                    self.pmra.val    /= 1000.
                    if self.pmra.uncert != -1:
                        self.pmra.uncert /= 1000.
                    self.pmra.unit    = 'arcsec/year'

                # already arcsec
                elif self.pmra.unit.lower() == 'arcsec/year':
                    pass
                else:
                    raise UnitError("unrecognized units for self.pmra!")

            # check pmDEC
            if self.pmdec.val != -1:
                # mas to arcsec
                if self.pmdec.unit.lower() == 'marcsec/year' or \
                       self.pmdec.unit.lower() == 'mas/year' or \
                       self.pmdec.unit.lower() == 'mas/yr':
                    self.pmdec.val    /= 1000.
                    if self.pmdec.uncert != -1:
                        self.pmdec.uncert /= 1000.
                    self.pmdec.unit    = 'arcsec/year'

                # already arcsec
                elif self.pmdec.unit.lower() == 'arcsec/year':
                    pass
                else:
                    raise UnitError("unrecognized units for self.pmdec!")

            # check planet radius
            if self.rp.val != -1:
                # jupiter radii to m
                if self.rp.unit.lower() == 'rjup':
                    self.rp.val    *= rjup
                    if self.rp.uncert != -1:
                        self.rp.uncert *= rjup
                    self.rp.unit    = 'm'

                # earth radii to m
                elif self.rp.unit.lower() == 'rearth':
                    self.rp.val    *= rearth
                    if self.rp.uncert != -1:
                        self.rp.uncert *= rearth
                    self.rp.unit    = 'm'

                # already m
                elif self.rp.unit.lower() == 'm':
                    pass
                else:
                    raise UnitError("unrecognized units for self.rp!")

            # check planet mass
            if self.mp.val != -1:
                # jupiter mass to kg
                if self.mp.unit.lower() == 'mjup':
                    self.mp.val    *= mjup
                    if self.mp.uncert != -1:
                        self.mp.uncert *= mjup
                    self.mp.unit    = 'kg'

                # earth mass to kg
                elif self.mp.unit.lower() == 'mearth':
                    self.mp.val    *= mearth
                    if self.mp.uncert != -1:
                        self.mp.uncert *= mearth
                    self.mp.unit    = 'kg'

                # already kg
                elif self.mp.unit.lower() == 'kg':
                    pass
                else:
                    raise UnitError("unrecognized units for self.mp!")

            # check period
            if self.period.val != -1:
                # days to sec
                if self.period.unit.lower() == 'days':
                    self.period.val    *= d2s
                    if self.period.uncert != -1:
                        self.period.uncert *= d2s
                    self.period.unit    = 'sec'

                # hrs to sec
                elif self.period.unit.lower() == 'hr':
                    self.period.val    *= hr2s   # seconds in an hour
                    if self.period.uncert != -1:
                        self.period.uncert *= hr2s
                    self.period.unit    = 'sec'

                # already sec
                elif self.period.unit.lower() == 'sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.period!")

            # check transdur
            if self.transdur.val != -1:
                # days to sec
                if self.transdur.unit.lower() == 'days':
                    self.transdur.val    *= d2s
                    if self.transdur.uncert != -1:
                        self.transdur.uncert *= d2s
                    self.transdur.unit    = 'sec'

                # hr to sec
                elif self.transdur.unit.lower() == 'hr':
                    self.transdur.val    *= hr2s
                    if self.transdur.uncert != -1:
                        self.transdur.uncert *= hr2s
                    self.transdur.unit    = 'sec'

                # min to sec
                elif self.transdur.unit.lower() == 'min':
                    self.transdur.val    *= 60.   # seconds per minute
                    if self.transdur.uncert != -1:
                        self.transdur.uncert *= 60.
                    self.transdur.unit    = 'sec'

                # already sec
                elif self.transdur.unit.lower() == 'sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.transdur!")

            # check limbtime
            if self.translimbtime.val != -1:
                # days to sec
                if self.translimbtime.unit.lower() == 'days':
                    self.translimbtime.val    *= d2s
                    if self.translimbtime.uncert != -1:
                        self.translimbtime.uncert *= d2s
                    self.translimbtime.unit    = 'sec'

                # hrs to sec
                elif self.translimbtime.unit.lower() == 'hr':
                    self.translimbtime.val    *= hr2s
                    if self.translimbtime.uncert != -1:
                        self.translimbtime.uncert *= hr2s
                    self.translimbtime.unit    = 'sec'

                # mins to sec
                elif self.translimbtime.unit.lower() == 'min':
                    self.translimbtime.val    *= 60.
                    if self.translimbtime.uncert != -1:
                        self.translimbtime.uncert *= 60.
                    self.translimbtime.unit    = 'sec'

                # already sec
                elif self.translimbtime.unit.lower() == 'sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.translimbtime!")

            # check semimaj axis
            if self.a.val != -1:
                # au to m
                if self.a.unit.lower() == 'au':
                    self.a.val    *= au2m
                    if self.a.uncert != -1:
                        self.a.uncert *= au2m
                    self.a.unit    = 'm'

                # km to m
                elif self.a.unit.lower() == 'km':
                    self.a.val    *= 1000.
                    if self.a.uncert != -1:
                        self.a.uncert *= 1000.
                    self.a.unit    = 'm'

                # rs to m
                elif self.a.unit.lower() == 'rs':
                    self.a.val    *= self.rs.val
                    if self.a.uncert != -1:
                        self.a.uncert *= self.rs.val
                    self.a.unit    = 'm'

                # already m
                elif self.a.unit.lower() == 'm':
                    pass
                else:
                    raise UnitError("unrecognized units for self.a!")

            # check ecldur
            if self.ecldur.val != -1:
                # days to sec
                if self.ecldur.unit.lower() == 'days':
                    self.ecldur.val    *= d2s
                    if self.ecldur.uncert != -1:
                        self.ecldur.uncert *= d2s
                    self.ecldur.unit    = 'sec'

                # hr to sec
                elif self.ecldur.unit.lower() == 'hr':
                    self.ecldur.val    *= hr2s
                    if self.ecldur.uncert != -1:
                        self.ecldur.uncert *= hr2s
                    self.ecldur.unit    = 'sec'

                # already sec
                elif self.ecldur.unit.lower() == 'sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.ecldur!")

            # check ecllimbtime
            if self.ecllimbtime.val != -1:
                # days to sec
                if self.ecllimbtime.unit.lower() == 'days':
                    self.ecllimbtime.val    *= d2s
                    if self.ecllimbtime.uncert != -1:
                        self.ecllimbtime.uncert *= d2s
                    self.ecllimbtime.unit    = 'sec'

                # hr to sec
                elif self.ecllimbtime.unit.lower() == 'hr':
                    self.ecllimbtime.val    *= hr2s
                    if self.ecllimbtime.uncert != -1:
                        self.ecllimbtime.uncert *= hr2s
                    self.ecllimbtime.unit    = 'sec'

                # already sec
                elif self.ecllimbtime.unit.lower() == 'sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.ecllimbtime!")

            # check inclination
            if self.i.val != -1:
                # deg to rag
                if self.i.unit.lower() == 'deg':
                    self.i.val = self.i.val * np.pi / 180.0
                    if self.i.uncert != -1:
                        self.i.uncert = self.i.uncert * np.pi / 180.0
                    self.i.unit = 'rad'

                # already rad
                elif self.i.unit.lower() == 'rad':
                    pass
                else:
                    raise UnitError("unrecognized units for self.i!")


            # check RVK
            if self.rvk.val != -1:
                # km/sec to m/sec
                if self.rvk.unit.lower() == 'km/sec':
                    self.rvk.val    *= 1000.
                    if self.rvk.uncert != -1:
                        self.rvk.uncert *= 1000.
                    self.rvk.unit    = 'm/sec'

                # already m/sec
                elif self.rvk.unit.lower() == 'm/sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.rvk!")

            # check RVgamma
            if self.rvgamma.val != -1:
                # km/sec to m/sec
                if self.rvgamma.unit.lower() == 'km/sec':
                    self.rvgamma.val    *= 1000.
                    if self.rvgamma.uncert != -1:
                        self.rvgamma.uncert *= 1000.
                    self.rvgamma.unit    = 'm/sec'

                # already m/sec
                elif self.rvgamma.unit.lower() == 'm/sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.rvgamma!")

            # check RVgammadot
            if self.rvgammadot.val != -1:
                # m/sec/yr to m/sec/sec
                if self.rvgammadot.unit.lower() == 'm/sec/yr':
                    self.rvgammadot.val /= 31556926.   # seconds in a year, wikipedia
                    if self.rvgammadot.uncert != -1:
                        self.rvgammadot.uncert /= 31556926.
                    self.rvgammadot.unit    = 'm/sec/sec'

                # m/sec/day to m/sec/sec
                elif self.rvgammadot.unit.lower() == 'm/sec/day':
                    self.rvgammadot.val /= d2s
                    if self.rvgammadot.uncert != -1:
                        self.rvgammadot.uncert /= d2s
                    self.rvgammadot.unit    = 'm/sec/sec'

                # km/sec/yr to m/sec/sec
                elif self.rvgammadot.unit.lower() == 'km/sec/yr':
                    self.rvgammadot.val *= 1000.
                    self.rvgammadot.val /= 31556926.
                    if self.rvgammadot.uncert != -1:
                        self.rvgammadot.uncert *= 1000.
                        self.rvgammadot.uncert /= 31556926.
                    self.rvgammadot.unit    = 'm/sec/sec'

                # km/sec/day to m/sec/sec
                elif self.rvgammadot.unit.lower() == 'km/sec/day':
                    self.rvgammadot.val *= 1000.
                    self.rvgammadot.val /= d2s
                    if self.rvgammadot.uncert != -1:
                        self.rvgammadot.uncert *= 1000.
                        self.rvgammadot.uncert /= d2s
                    self.rvgammadot.unit    = 'm/sec/sec'

                # already m/sec/sec
                elif self.rvgammadot.unit.lower() == 'm/sec/sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.rvgammadot!")

            # check RVvsinI
            if self.rvvsini.val != -1:
                # km/sec to m/sec
                if self.rvvsini.unit.lower() == 'km/sec':
                    self.rvvsini.val    *= 1000.
                    if self.rvvsini.uncert != -1:
                        self.rvvsini.uncert *= 1000.
                    self.rvvsini.unit    = 'm/sec'

                # already m/sec
                elif self.rvvsini.unit.lower() == 'm/sec':
                    pass
                else:
                    raise UnitError("unrecognized units for self.rvvsini!")

            # check OMEGA
            if self.omega.val != -1:
                # deg to rad
                if self.omega.unit.lower() == 'deg':
                    self.omega.val *= np.pi/180.
                    if self.omega.uncert != -1:
                        self.omega.uncert *= np.pi/180.
                    self.omega.unit    = 'rad'

                # already rad
                elif self.omega.unit.lower() == 'rad' or \
                         self.omega.unit.lower() == 'radians':
                    if self.omega.unit.lower() == 'radians':
                        self.omega.unit = 'rad'
                    pass
                else:
                    raise UnitError("unrecognized units for self.omega!")


            # check pmRA
            if self.pmra.val != -1:
                # mas/year
                if self.pmra.unit.lower() == 'mas/year' or \
                       self.pmra.unit.lower() == 'mas/yr' or \
                       self.pmra.unit.lower() == 'marcsec/year' or \
                       self.pmra.unit.lower() == 'marcsec/yr':
                    self.pmra.val /= 1000.
                    if self.pmra.uncert != -1:
                        self.pmra.uncert /= 1000.
                    self.pmra.unit    = 'arcsec/year'

                # already arcsec/year
                elif self.pmra.unit.lower() == 'arcsec/year':
                    pass
                else:
                    raise UnitError("unrecognized units for self.pmra!")

            # check pmDEC
            if self.pmdec.val != -1:
                # mas/year
                if self.pmdec.unit.lower() == 'mas/year' or \
                       self.pmdec.unit.lower() == 'mas/yr' or \
                       self.pmdec.unit.lower() == 'marcsec/year' or \
                       self.pmdec.unit.lower() == 'marcsec/yr':
                    self.pmdec.val /= 1000.
                    if self.pmdec.uncert != -1:
                        self.pmdec.uncert /= 1000.
                    self.pmdec.unit    = 'arcsec/year'

                # already arcsec/year
                elif self.pmdec.unit.lower() == 'arcsec/year':
                    pass
                else:
                    raise UnitError("unrecognized units for self.pmdec!")

        # do conversions for versions < 4
        elif self.version < 4:
            # convert stellar radius to meters
            if self.rs.val != -1:
                self.rs.val    *= rsun
                if self.rs.uncert != -1:
                    self.rs.uncert *= rsun
            # convert planet radius to meters
            if self.rp.val != -1:
                self.rp.val    *= rjup
                if self.rp.uncert != -1:
                    self.rp.uncert *= rjup
            # convert planet mass to kg
            if self.mass.val != -1:
                self.mass.val    *= mjup
                if self.mass.uncert != -1:
                    self.mass.uncert *= mjup
            # convert period to seconds
            if self.period.val != -1:
                self.period.val    *= d2s
                if self.period.uncert != -1:
                    self.period.uncert *= d2s
            # convert semimajor axis to meters
            if self.a.val != -1:
                self.a.val    *= au2m
                if self.a.uncert != -1:
                    self.a.uncert *= au2m
            # convert omega to radians
            if self.omega.val != -1:
                self.omega.val    *= np.pi/180.
                if self.omega.uncert != -1:
                    self.omega.uncert *= np.pi/180.
            # convert inclination to radians
            if self.i.val != -1:
                self.i.val    *= np.pi/180.
                if self.i.uncert != -1:
                    self.i.uncert *= np.pi/180.
