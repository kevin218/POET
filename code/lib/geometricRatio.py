
#Estimate the geometric ratio (R1/R2*d2/d1) of two stellar sources
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import centerdriver as cd
import apphot as ap
import kurucz_inten as kur
import imageedit   as ie
import gaussian    as g
from tqdm import tqdm
#from importlib import reload
#reload(kur)


#Modified center driver
def cdriver(method, data, guess, trim, radius, size,
                 mask=None, uncd=None, fitbg=1, maskstar=True,
                 expand=5.0, psf=None, psfctr=None):
  # Default mask: all good
  if mask is None:
    mask = np.ones(np.shape(data))
  # Default uncertainties: flat image
  if uncd is None:
    uncd = np.ones(np.shape(data))
  # Trim the image if requested
  if trim != 0:
    # Integer part of center
    cen = np.rint(guess)
    # Center in the trimed image
    loc = (trim, trim)
    # Do the trim:
    img, msk, err = ie.trimimage(data, cen, loc, mask=mask, uncd=uncd)
  else:
    cen = np.array([0,0])
    loc = np.rint(guess)
    img, msk, err = data, mask, uncd
  # If all data is bad:
  if not np.any(msk):
    raise Exception('Bad Frame Exception!')
  weights = 1.0/np.abs(err)
  # Get the center with one of the methods:
  if   method == 'fgc':
    foo, bar = g.fitgaussian(img, yxguess=loc, mask=msk, weights=weights,
                         fitbg=fitbg, maskg=maskstar)
    #print(foo, bar)
    y, x = foo[2:4]
    yerr, xerr = bar[2:4]

  # Make trimming correction and return
  return ((y, x) + cen - trim), (yerr, xerr)

def compute_df(ev, guess_sec, isplots=1):
    '''
    Compute dilution factor of mean image
    '''
    
    #filterwave  = 13187.5*1e-4  #wavelength in microns 
    #filterhw    = 80.5*1e-4     #half-width


    # Calculate RA and DEC
    # foo     = hdulist['WCSCORR'].data
    #FINDME

    im      = np.squeeze(ev.meanim)
    psfim   = ev.psfim
    
    # subtract background
    # im     -= np.median(im)

    # Guess
    guess_pri = (ev.targpos[0][0], ev.targpos[1][0])
    #guess_sec = (154.0, 31.0)
    guess_psf = ev.psfctr
    
    # Create mask
    mask_sec = np.ones(im.shape)

    # Centering using 2D Gaussian
    #ctr_psf2,extra = cd.centerdriver('fgc', psfim2, guess_psf, trim=5, radius=None, size=None)
    ctr_pri, extra = cdriver('fgc', im, guess_pri, trim=3, radius=None, size=None) #, uncd=imerr
    mask_sec[int(round(ctr_pri[0]))-1:int(round(ctr_pri[0]))+2,int(round(ctr_pri[1]))-1:int(round(ctr_pri[1]))+2] = 0
    print("Masking 3x3 region around primary star for centroiding of secondary")
    ctr_sec, extra = cdriver('fgc', im, guess_sec, trim=3, mask=mask_sec, radius=None, size=None)
    ctr_psf, extra = cdriver('fgc', psfim, guess_psf, trim=5, radius=None, size=None)
    # Spitzer pixel scale ~1.2 arcsec/pixel
    separation = np.sqrt(1.2**2*(ctr_pri[0]-ctr_sec[0])**2 + 1.2**2*(ctr_pri[1]-ctr_sec[1])**2) #arcsec
    #vary = ctr_w12err[0]**2 + ctr_b6err[0]**2
    #varx = ctr_w12err[1]**2 + ctr_b6err[1]**2
    #varz = vary*4*0.121**4*(ctr_w12[0]-ctr_b6[0])**2 + varx*4*0.135**4*(ctr_w12[1]-ctr_b6[1])**2
    #seperr = 0.5*np.sqrt(varz)/separation
    print(f'Primary postion: {ctr_pri} pixels')
    print(f'Secondary postion: {ctr_sec} pixels')
    print(f'Separation: {separation} arcsec')
    print(f'PSF postion: {ctr_psf} pixels')
    
    plt.rcParams.update({'legend.fontsize':12})
    if isplots >= 3:
        #Plot image with centers
        plt.figure(3, figsize=(6,4.8))
        plt.clf()
        #vmax = im[int(ctr_sec[0]),int(ctr_sec[1])]
        vmax = im[int(round(ctr_pri[0])),int(round(ctr_pri[1]))]
        plt.imshow(im, origin='lower', aspect='auto', vmin=0, vmax=vmax, cmap=plt.cm.gray_r, interpolation='nearest')
        plt.plot([ctr_pri[1]], [ctr_pri[0]], '+', ms=10, mew=3, label='Primary')
        plt.plot([ctr_sec[1]], [ctr_sec[0]], 'x', ms=10, mew=2, label='Secondary')
        #plt.ylim(140,160)
        #plt.xlim(15,40)
        plt.xticks(size=12)
        plt.yticks(size=12)
        plt.ylabel('Pixel Number', size=14)
        plt.xlabel('Pixel Number', size=14)
        plt.legend(loc='lower right')
        plt.tight_layout()
        #plt.subplots_adjust(0.11,0.10,0.97,0.98)
        plt.savefig(f'{ev.eventname}-Spitzer-Centroids.png')


    # Perform photometry for different aperature sizes
    photap_pri      = np.arange(1.0,4.6,0.1)
    nphot_pri       = len(photap_pri)
    flux_pri     = np.zeros(nphot_pri)
    flux_psf_pri = np.zeros(nphot_pri)
    for i in tqdm(range(nphot_pri)):
        flux_pri[i]  = ap.apphot(im,  ctr_pri ,  photap_pri[i],  skyin=7,    skyout=15, expand=5, betahw=None, targpos=None)
        flux_psf_pri[i] = ap.apphot(psfim,  ctr_psf,  photap=photap_pri[i]*ev.psfexpand,  skyin=7*ev.psfexpand,    skyout=15*ev.psfexpand, expand=3, betahw=None, targpos=None)

    """
    photap2      = np.arange(1,4.1,0.1)
    nphot2       = len(photap2)
    flux_b6      = np.zeros(len(photap2))
    flux_psf_b6  = np.zeros(len(photap2))
    for i in range(nphot2):
        flux_b6[i]   = ap.apphot(im,  ctr_b6 ,  photap2[i],  skyin=12,    skyout=20, expand=5)
        flux_psf_b6[i]  = ap.apphot(psfim2,  ctr_psf2,  photap=photap2[i],  skyin=12,    skyout=20, expand=5)


    """
    '''
    plt.figure(3)
    plt.clf()
    imslice = im[150,11:36]
    plt.plot(imslice/imslice.sum(),'bo-')
    plt.plot(psfim1[13]*(imslice/imslice.sum()).max()/psfim1.max(),'ko-')
    '''

    photap_sec      = np.arange(1.0,3.1,0.1)
    nphot_sec       = len(photap_sec)
    flux_sec      = np.zeros(nphot_sec)
    flux_psf_sec  = np.zeros(nphot_sec)
    for i in tqdm(range(nphot_sec)):
        flux_sec[i]   = ap.apphot(im,  ctr_sec ,  photap_sec[i],  skyin=7,    skyout=15, expand=5, betahw=None, targpos=None)
        flux_psf_sec[i]  = ap.apphot(psfim,  ctr_psf,  photap=photap_sec[i]*ev.psfexpand,  skyin=7*ev.psfexpand,    skyout=15*ev.psfexpand, expand=3, betahw=None, targpos=None)
    #print(flux_psf_sec)

    if isplots >= 3:
        plt.figure(4)
        plt.clf()
        plt.plot(photap_pri, flux_pri/flux_psf_pri/(flux_pri[0]/flux_psf_pri[0]), 'o-', label='Primary')
        plt.plot(photap_sec, flux_sec/flux_psf_sec/(flux_sec[0]/flux_psf_sec[0]), 'o-', label='Secondary')
        #plt.plot(photap2, flux_psf_b6/flux_psf_b6[0], 'b-')
        #plt.plot(photap1, flux_psf_w12/flux_psf_w12[0]*1.18, 'k-')
        plt.legend(loc='lower right')
        plt.xlabel('Aperture Size [Pixels]', size=14)
        plt.ylabel('Normalized Flux', size=14)
        plt.savefig(f'{ev.eventname}-Spitzer-AccumulatedFlux.png')
    
    bb    = np.ones((nphot_pri,nphot_sec))*flux_sec/flux_psf_sec
    aa    = np.ones((nphot_sec,nphot_pri))*flux_pri/flux_psf_pri
    df    = bb.T/aa
    #dferr = df*np.sqrt((photerr_w12/phot_w12)**2 + (photerr_b6/phot_b6)**2)

    if isplots >= 1:
        plt.figure(1, figsize=(8,5))
        plt.clf()
        #a = plt.axes([0.1,0.06,0.93,0.93])
        plt.imshow(df, origin='lower', aspect='auto', interpolation='nearest')
        plt.xticks(np.arange(0,nphot_pri,2), np.round(photap_pri[::2],2))
        plt.yticks(np.arange(0,nphot_sec,2), np.round(photap_sec[::2],2))
        plt.xlabel('Primary Aperture Size [Pixels]', size=14)
        plt.ylabel('Secondary Aperture Size [Pixels]', size=14)
        plt.colorbar(shrink=0.86)
        plt.tight_layout()
        plt.savefig(f'{ev.eventname}-Spitzer-DilutionFactors.png')
    
    return df


##Best apertures
#phot_w12, photerr_w12 = ap.apphot(im,  ctr_w12,  photap=3.,  skyin=12,    skyout=20, expand=5, imerr=imerr, aperr=True)
#phot_b6 , photerr_b6  = ap.apphot(im,  ctr_b6 ,  photap=2.,  skyin=12,    skyout=20, expand=5, imerr=imerr, aperr=True)
#phot_psf_w12 = ap.apphot(psfim1,  ctr_psf1,  photap=3.,  skyin=12,    skyout=20, expand=5)
#phot_psf_b6  = ap.apphot(psfim2,  ctr_psf2,  photap=2.,  skyin=12,    skyout=20, expand=5)

#aa = phot_w12/phot_psf_w12
#bb = phot_b6/phot_psf_b6
#fluxRatio    = bb/aa
#fluxRatioerr = fluxRatio*np.sqrt((photerr_w12/phot_w12)**2 + (photerr_b6/phot_b6)**2)
##0.06923 +/- 0.00148

##Plot DF vs WASP-12 aperture size
#iphot2 = 10
#plt.figure(5)
#plt.clf()
#plt.plot(photap1, df[iphot2], 'bo-')
#plt.errorbar([3.0], [fluxRatio], [fluxRatioerr], color='r')
#plt.ylabel('Dilution Factor')
#plt.xlabel('WASP-12 Aperture Size [Pixels]')
#plt.savefig('DF-vs-Aperture.png')
#print(np.mean(df[iphot2,10:]),np.std(df[iphot2,10:]))

def kurucz_ratio(ev):
    '''
    
    '''
    # Kurucz model
    kuruczfile  = ev.kuruczfile
    temp_pri    = 5830  #Evans (2016)
    temp_prierr =  100  #Evans (2016)
    grav_pri    =  4.5
    temp_sec    = 4810  #Evans (2016)
    temp_secerr =  100  #Evans (2016), 
    #temp_sec    = 4570  #Evans (2018), https://arxiv.org/pdf/1709.07476.pdf
    #temp_secerr =  240  #Evans (2018), https://arxiv.org/pdf/1709.07476.pdf
    grav_sec    =  4.5
    numit       = 5000
    filterwave  = 3.6
    filterhw    = 0.4

    # Read Kurucz file
    inten, wave, grav, temp, nainten, head = kur.read(kuruczfile)
    wave   *= 1e6     #Convert to microns
    iwave   = np.where(np.bitwise_and(wave >= (filterwave-filterhw), wave <= (filterwave+filterhw)))[0]

    # Calculate best spectra
    bestinten_pri = kur.interp(inten[:,iwave], grav, temp, grav_pri, temp_pri)
    bestinten_sec = kur.interp(inten[:,iwave], grav, temp, grav_sec , temp_sec )

    ## Generate random distribution for impacting parameters (not gravity)
    #tw12    = np.random.normal(temp_w12 , temp_w12err, numit)
    #tb6     = np.random.normal(temp_b6  , temp_b6err , numit)
    ##geor    = np.random.normal(geo_ratio, gratioerr  , numit)

    ## Calculate distribution of spectra
    #inten_w12 = np.zeros((numit, len(bestinten_w12)))
    #inten_b6  = np.zeros((numit, len(bestinten_b6 )))
    #for i in range(numit):
    #    sys.stdout.write('\r'+str(i+1)+'/'+str(numit))
    #    sys.stdout.flush()
    #    inten_w12[i] = kur.interp(inten[:,iwave], grav, temp, grav_w12, tw12[i])
    #    inten_b6 [i] = kur.interp(inten[:,iwave], grav, temp, grav_b6 , tb6 [i])

    # Estimate model ratio
    kuruczRatio    = np.mean(bestinten_sec/bestinten_pri)
    print(f'Kurucz Ratio: {kuruczRatio}')
    #inten_w12err   = np.mean(np.std(inten_w12,axis=0))
    #inten_b6err    = np.mean(np.std(inten_b6,axis=0))
    #kuruczRatioerr = kuruczRatio*np.sqrt((inten_w12err/bestinten_w12.mean())**2 + (inten_b6err/bestinten_b6.mean())**2)

    ## Estimate geometric ratio
    #geoRatio    = np.sqrt(fluxRatio/kuruczRatio/2.)
    #geoRatioerr = 0.5*geoRatio*np.sqrt((fluxRatioerr/fluxRatio)**2 + (kuruczRatioerr/kuruczRatio)**2)
    ##0.3568 +/- 0.0038      Assuming no errors in Kurucz models
    ##0.3571 +/- 0.0130      Using uncertainties in stellar temperatures

    ## Effective correction factor (2*f^2)
    #corrFactor    = fluxRatio/kuruczRatio
    #corrFactorerr = corrFactor* np.sqrt((fluxRatioerr/fluxRatio)**2 + (kuruczRatioerr/kuruczRatio)**2)
    ##0.2550 +/- 0.0186
    
    return kuruczRatio

    ##Single stellar companion
    #R_b6  = 0.55
    #R_w12 = 1.57
    #distance = R_b6/R_w12/geoRatio/np.sqrt(2)
    ##0.694

    ##Binary stellar companion
    #R_b6  = 0.55
    #R_w12 = 1.57
    #R_w12err = 0.07
    #distance = R_b6/R_w12/geoRatio
    ##0.981
    #calcR_b6 = geoRatio*R_w12
    #calcR_b6err = calcR_b6*np.sqrt((R_w12err/R_w12)**2 + (geoRatioerr/geoRatio)**2)
    ##0.56 +/- 0.03

