import numpy as np
import pylab as plt
import pyfits as pf

ch1hdu = pf.open('WASP-18b_secondary_2008-12-20_Spitzer_IRAC_3.6_microns.fits')
ch2hdu = pf.open('WASP-18b_secondary_2008-12-20_Spitzer_IRAC_5.8_microns.fits')
ch3hdu = pf.open('WASP-18b_secondary_2008-12-24_Spitzer_IRAC_4.5_microns.fits')
ch4hdu = pf.open('WASP-18b_secondary_2008-12-24_Spitzer_IRAC_8.0_microns.fits')

ch1data = ch1hdu[1].data
ch2data = ch2hdu[1].data
ch3data = ch3hdu[1].data
ch4data = ch4hdu[1].data

# The top comment for the table. It should describe the data.
topstring = 'This file contains the lightcurve from the paper:\n \
             Analysis of Exoplanet HD 149026b Using BLISS Mapping\n \
             by Kevin B. Stevenson, Joseph Harrington, Jonathan Fortney,\n \
             Ryan A. Hardy, Sarah Nymeyer, William C. Bowman, Patricio Cubillos,\n \
             M. Oliver Bowman, Matthew Hardin, and Thomas J. Loredo\n \
             which was published in 2011 in ApJ.\n \
             The data are from the Infrared Array Camera and the Infrared\n \
             Spectrograph''s photometric blue peak-up array on the US National\n \
             Aeronautics and Space Administrations Spitzer Space Telescope,\n \
             programs 254, 40135, 50517, and are available from the public\n \
             Spitzer archive (http://spitzer.caltech.edu).\n \
             The paper cited above and its electronic supplement, of which this\n \
             file is a part, describe the observations and data analysis.\n \
             The data below are in IPAC table format. The keywords describe the\n \
             columns in the table. All series (pixel position, frame number,\n \
             etc.) are zero-based. Data in a table row are valid if GOOD equals 1,\n \
             invalid otherwise. ORIGIN originally read "Spitzer Science Center",\n \
             but has been updated to reflect the origin of the present file.'
    
import irsa as irsa
#CREATE IRSA TABLES
for i in range(nfit):
    irsa.old_irsa(event[i], fit[i], topstring, obj[i]+'_irsa.tbl')
    #irsa.do_irsa(event[i], fit[i], topstring, obj[i]+'_irsa.tbl')
    

