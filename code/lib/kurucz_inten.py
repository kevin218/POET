import sys

"""
def interp3d(inten, grav, temp, wgrav, wtemp, log=None, kf=None):
    '''
    This function finds the brightness spectrum of a star by
    interpolating temperature and gravity in a grid of Kurucz
    models.  Wavelengths/frequencies are not interpolated.

    INPUTS
        Inten:    [nwavl, nmod] array of model brightnesses
            (erg cm-2 s-1 sr-1 Hz-1).

        Grav:    [nmod] array of the base-10 log of stellar surface
            gravities for the models (cm s-2).

        Temp:    [nmod] array of temperatures for the models (K).

        Wgrav:    Wanted log(gravity) (cm s-2).

        Wtemp:    Wanted temperature (K).

    KEYWORD PARAMETERS:
        LOG:    If set, interpolate in the log.
        KF:    [nwavl, ngrav, ntemp] array of Kurucz models (returned)
            (returned, erg cm-2 s-1 sr-1 Hz-1).

    OUTPUTS:
        This function returns an [nwavl] array of brightnesses
        (erg cm-2 s-1 sr-1 Hz-1).

    PROCEDURE:
        Uses cubic convolution in 2D, done explicitly at each
        wavelength.  Note that /log makes very little difference (see
        example).

        Kurucz's modeled grav and temp do not cover a rectangular
        grid, since he does not calculate low gravities for large
        stars that never have them.  The code currently just takes the
        rectangular region filled by his cooler models and makes a
        grid out of them, ignoring the rest.  If your parameters are
        too close to the edge of this grid, the code warns and exits.
        It will not be hard to recode to accomodate more models.

    EXAMPLE:
    

    MODIFICATION HISTORY:
         Written by:    Joseph Harrington, Cornell.  2006-02-02
                jh@oobleck.astro.cornell.edu
        2006-02-03 jh    Swapped wgrav and wtemp in calling order to
                match grav and temp.  Updated example for
                change of kurucz_flux_read to MKS.
        2006-02-10 jh    Fixed example, added kf keyword, changed name
                to kurucz_inten_interp.
        2006-02-11 jh    Header tweak.
        2006-02-12 jh    Converted variable names flux -> inten.
                Changed messages to printed warnings, since
                they make Monte Carlo trials bomb.
        2008-08-11 kevin
                Converted to Python
        2015-05-28 kevin
                Added 3D functionality.  Adapted from interp().
    '''
    import numpy as np
    import scipy.interpolate as spi
    
    niter  = wgrav.shape[0]
    dim    = inten.shape
    nwav   = dim[1]
    iinten = np.zeros((niter,nwav))

    # FINDME: HARDWIRED!
    ilograv = 0
    ihigrav = 10 # the full range of gravities
    ilotemp = 0
    ihitemp = 16 # the temp where Kurucz starts skipping gravities
    ngrav = ihigrav - ilograv + 1
    ntemp = ihitemp - ilotemp + 1

    kg = grav[ilograv:ihigrav+1]
    it = np.arange(ntemp) * ngrav + ilotemp
    kt = temp[it]
    
    # FINDME: this only works if the low indices are both 0 and we're not
    # at temps where Kurucz starts skipping gravities...
    if kg[ihigrav] < wgrav.min():   
          print 'kurucz_inten.interp: wgrav ' + str(wgrav.min()) + ' higher than ' + str(kg[ihigrav]) + ', see FINDME in code to fix.'
    if kt[ihitemp] < wtemp.min():   
          print 'kurucz_inten.interp: wtemp ' + str(wtemp.min()) + ' higher than ' + str(kt[ihitemp]) + ', see FINDME in code to fix.'
    kf = np.reshape(inten[0:ngrav * ntemp,:], (ntemp, ngrav, nwav))

    # calculate the indices of the parameters we want
    itemp, igrav = np.mgrid[kt[0]:kt[-1]+1:250,kg[0]:kg[-1]+0.1:0.5]

    if (log is not None):   
          kf = np.log(kf)
    
    print('Interpolating Kurucz spectrum...')
    #FINDME: meshgrid is incorrect...
    iinten = spi.interpn([itemp, igrav, range(nwav)], kf, np.meshgrid((wtemp, wgrav, range(nwav)))])


    if (log is not None):   
          iinten = exp(iinten)
    return iinten
"""

def interp2d(inten, grav, temp, wgrav, wtemp, log=None, kf=None):
    '''
    This function finds the brightness spectrum of a star by
    interpolating temperature and gravity in a grid of Kurucz
    models.  Wavelengths/frequencies are not interpolated.

    INPUTS
        Inten:    [nwavl, nmod] array of model brightnesses
            (erg cm-2 s-1 sr-1 Hz-1).

        Grav:    [nmod] array of the base-10 log of stellar surface
            gravities for the models (cm s-2).

        Temp:    [nmod] array of temperatures for the models (K).

        Wgrav:    Wanted log(gravity) (cm s-2).

        Wtemp:    Wanted temperature (K).

    KEYWORD PARAMETERS:
        LOG:    If set, interpolate in the log.
        KF:    [nwavl, ngrav, ntemp] array of Kurucz models (returned)
            (returned, erg cm-2 s-1 sr-1 Hz-1).

    OUTPUTS:
        This function returns an [nwavl] array of brightnesses
        (erg cm-2 s-1 sr-1 Hz-1).

    PROCEDURE:
        Uses cubic convolution in 2D, done explicitly at each
        wavelength.  Note that /log makes very little difference (see
        example).

        Kurucz's modeled grav and temp do not cover a rectangular
        grid, since he does not calculate low gravities for large
        stars that never have them.  The code currently just takes the
        rectangular region filled by his cooler models and makes a
        grid out of them, ignoring the rest.  If your parameters are
        too close to the edge of this grid, the code warns and exits.
        It will not be hard to recode to accomodate more models.

    EXAMPLE:
    

    MODIFICATION HISTORY:
         Written by:    Joseph Harrington, Cornell.  2006-02-02
                jh@oobleck.astro.cornell.edu
        2006-02-03 jh    Swapped wgrav and wtemp in calling order to
                match grav and temp.  Updated example for
                change of kurucz_flux_read to MKS.
        2006-02-10 jh    Fixed example, added kf keyword, changed name
                to kurucz_inten_interp.
        2006-02-11 jh    Header tweak.
        2006-02-12 jh    Converted variable names flux -> inten.
                Changed messages to printed warnings, since
                they make Monte Carlo trials bomb.
        2008-08-11 kevin
                Converted to Python
        2015-05-19 kevin
                Added 2D functionality.  Adapted from interp().
    '''
    import numpy as np
    import scipy.interpolate as spi
    
    niter  = wgrav.shape[0]
    dim    = inten.shape
    nwav   = dim[1]
    iinten = np.zeros((niter,nwav))

    # FINDME: HARDWIRED!
    ilograv = 0
    ihigrav = 10 # the full range of gravities
    ilotemp = 0
    ihitemp = 16 # the temp where Kurucz starts skipping gravities
    ngrav = ihigrav - ilograv + 1
    ntemp = ihitemp - ilotemp + 1

    kg = grav[ilograv:ihigrav+1]
    it = np.arange(ntemp) * ngrav + ilotemp
    kt = temp[it]

    #For higher temperatures where kurucz starts skipping gravities
    ilograv2 = 1 #grav 0.5-5.0
    ihigrav2 = 10
    ilotemp2 = 11
    ihitemp2 = 17
    ngrav2 = ihigrav2 - ilograv2 + 1
    ntemp2 = ihitemp2 - ilotemp2 + 1
    kg2 = grav[ilograv2:ihigrav2+1]
    it2 = np.arange(ntemp2) * ngrav2 + ilotemp2*ngrav
    kt2 = temp[it2]

    ilograv3 = 2 #grav 1.0-5.0
    ihigrav3 = 10
    ilotemp3 = 17
    ihitemp3 = 20
    ngrav3 = ihigrav3 - ilograv3 + 1
    ntemp3 = ihitemp3 - ilotemp3 + 1
    kg3 = grav[ilograv3:ihigrav3+1]
    it3 = np.arange(ntemp3)*ngrav3+(ilotemp3-ilotemp2)*ngrav2+ilotemp2*ngrav
    kt3 = temp[it3]

    ilograv4 = 3 #grav 1.5-5.0
    ihigrav4 = 10
    ilotemp4 = 20
    ihitemp4 = 23
    ngrav4 = ihigrav4 - ilograv4 + 1
    ntemp4 = ihitemp4 - ilotemp4 + 1
    kg4 = grav[ilograv4:ihigrav4+1]
    it4 = np.arange(ntemp4)*ngrav4+(ilotemp4-ilotemp3)*ngrav3+(ilotemp3-ilotemp2)*ngrav2+ilotemp2*ngrav
    kt4 = temp[it4]

    ilograv5 = 4 #grav 2.0-5.0
    ihigrav5 = 10
    ilotemp5 = 23
    ihitemp5 = 34
    ngrav5 = ihigrav5 - ilograv5 + 1
    ntemp5 = ihitemp5 - ilotemp5 + 1
    kg5 = grav[ilograv5:ihigrav5+1]
    it5 = np.arange(ntemp5)*ngrav5+(ilotemp5-ilotemp4)*ngrav4+(ilotemp4-ilotemp3)*ngrav3+(ilotemp3-ilotemp2)*ngrav2+ilotemp2*ngrav
    kt5 = temp[it5]
    
    # FINDME: this only works if the low indices are both 0 and we're not
    # at temps where Kurucz starts skipping gravities...
    if kg[ihigrav] < wgrav.min():   
          print('kurucz_inten.interp: wgrav ' + str(wgrav.min()) + ' higher than ' + str(kg[ihigrav]) + ', see FINDME in code to fix.')
    #if kt[ihitemp] < wtemp.min():   
    #      print('kurucz_inten.interp: wtemp ' + str(wtemp.min()) + ' higher than ' + str(kt[ihitemp]) + ', see FINDME in code to fix.')
    
    #isort = np.argsort(itemp)
    print('Interpolating Kurucz spectrum...')
    for i in range(nwav):
        sys.stdout.write('\r'+str(nwav)+', '+str(i)+'  ')
        sys.stdout.flush()
        #f = spi.RectBivariateSpline(itemp[isort], igrav[isort], kf[:,:,i], kx=1, ky=1, s=0)
        #iinten[:,i] = f.ev(wtemp, wgrav)
        
        #iinten[:,i] = spi.bisplev(wtemp, wgrav, tck)   
        #wtemp and wgrav need to be sorted; therefore, use for loop instead
        for j in range(niter):
            if kt[ihitemp] >= wtemp[j]:
                #print('Temp should be below 6250; real temp is'+str(wtemp[j])) #Megan debugging purposes only
                kf = np.reshape(inten[0:ngrav * ntemp,:], (ntemp, ngrav, nwav))

                # calculate the indices of the parameters we want
                itemp, igrav = np.mgrid[kt[0]:kt[-1]+1:250,kg[0]:kg[-1]+0.1:0.5]

                if (log is not None):   
                    kf = np.log(kf)

            elif kt2[-1] >= wtemp[j]:
                #print('Temp should be 6250-7750; real temp is'+str(wtemp[j]))
                kf = np.reshape(inten[ngrav*ntemp:(ngrav * ntemp+ngrav2*ntemp2),:], (ntemp2, ngrav2, nwav))

                # calculate the indices of the parameters we want
                itemp, igrav = np.mgrid[kt2[0]:kt2[-1]+1:250,kg2[0]:kg2[-1]+0.1:0.5]

                if (log is not None):   
                    kf = np.log(kf)

            elif kt3[-1] >= wtemp[j]:
                #print('Temp should be 7750-8500; real temp is'+str(wtemp[j]))
                kf = np.reshape(inten[(ngrav * ntemp+ngrav2*ntemp2):(ngrav * ntemp+ngrav2*ntemp2+ngrav3*ntemp3),:], (ntemp3, ngrav3, nwav))

                # calculate the indices of the parameters we want
                itemp, igrav = np.mgrid[kt3[0]:kt3[-1]+1:250,kg3[0]:kg3[-1]+0.1:0.5]

                if (log is not None):   
                    kf = np.log(kf)

            elif kt4[-1] >= wtemp[j]:
                #print('Temp should be 8500-9250; real temp is'+str(wtemp[j]))
                kf = np.reshape(inten[(ngrav * ntemp+ngrav2*ntemp2+ngrav3*ntemp3):(ngrav * ntemp+ngrav2*ntemp2+ngrav3*ntemp3+ngrav4*ntemp4),:], (ntemp4, ngrav4, nwav))

                # calculate the indices of the parameters we want
                itemp, igrav = np.mgrid[kt4[0]:kt4[-1]+1:250,kg4[0]:kg4[-1]+0.1:0.5]

                if (log is not None):   
                    kf = np.log(kf)

            elif kt5[-1] >= wtemp[j]:
                #print('Temp should be above 9250; real temp is'+str(wtemp[j]))
                kf = np.reshape(inten[(ngrav * ntemp+ngrav2*ntemp2+ngrav3*ntemp3+ngrav4*ntemp4):(ngrav * ntemp+ngrav2*ntemp2+ngrav3*ntemp3+ngrav4*ntemp4+ngrav5*ntemp5),:], (ntemp5, ngrav5, nwav))

                # calculate the indices of the parameters we want
                itemp, igrav = np.mgrid[kt5[0]:kt5[-1]+1:250,kg5[0]:kg5[-1]+0.1:0.5]

                if (log is not None):   
                    kf = np.log(kf)



            tck = spi.bisplrep(itemp, igrav, kf[:,:,i], kx=1, ky=1)

            iinten[j,i] = spi.bisplev(wtemp[j], wgrav[j], tck) 
        
    if (log is not None):   
          iinten = exp(iinten)

    return iinten


"""
 NAME:
    INTERP

 PURPOSE:
    This function finds the brightness spectrum of a star by
    interpolating temperature and gravity in a grid of Kurucz
    models.  Wavelengths/frequencies are not interpolated.

 CATEGORY:
    Astronomy

 CALLING SEQUENCE:

    Result = INTERP(Inten, Grav, Temp, Wgrav, Wtemp)

 INPUTS:
    Inten:    [nwavl, nmod] array of model brightnesses
        (erg cm-2 s-1 sr-1 Hz-1).

    Grav:    [nmod] array of the base-10 log of stellar surface
        gravities for the models (cm s-2).

    Temp:    [nmod] array of temperatures for the models (K).

    Wgrav:    Wanted log(gravity) (cm s-2).

    Wtemp:    Wanted temperature (K).

 KEYWORD PARAMETERS:
    LOG:    If set, interpolate in the log.
    KF:    [nwavl, ngrav, ntemp] array of Kurucz models (returned)
        (returned, erg cm-2 s-1 sr-1 Hz-1).

 OUTPUTS:
    This function returns an [nwavl] array of brightnesses
    (erg cm-2 s-1 sr-1 Hz-1).

 PROCEDURE:
    Uses cubic convolution in 2D, done explicitly at each
    wavelength.  Note that /log makes very little difference (see
    example).

    Kurucz's modeled grav and temp do not cover a rectangular
    grid, since he does not calculate low gravities for large
    stars that never have them.  The code currently just takes the
    rectangular region filled by his cooler models and makes a
    grid out of them, ignoring the rest.  If your parameters are
    too close to the edge of this grid, the code warns and exits.
    It will not be hard to recode to accomodate more models.

 EXAMPLE:
    ;Example is still in IDL:
    inten = kurucz_inten_read('fp00ak2odfnew.pck', $
                wave, grav, temp, nainten, head)
    ; find the Sun's spectrum
    wgrav = 4.44d ; log(cm s-2)
    wtemp = 5770d ; Kelvins
    kinten = kurucz_inten_interp(inten, grav, temp, wgrav, wtemp)
    c     = 2.99792458d8 ; speed of light in m/s
    freq  = reverse(c / wave)
    fsun  = reverse(kinten)
    plot, freq, fsun, /xlog, $
      xtitle='Frequency (Hz)', ytitle='I (W m-2 sr-1 Hz-1)'

    ; find the Sun's luminosity
    ; extra factor of !dpi converts brightness to flux
    print, int_tabulated(freq, fsun, /d) * !dpi * 4d * !dpi * 6.96d8^2
    ;   3.8392897e+26

    ; Compare to blackbody...
    pf = planckfreq(freq, wtemp) / !dpi * 1d-7 * 1d4 ; flux -> inten., MKS
    plot, freq, fsun, /xlog, $
      xtitle='Frequency (Hz)', ytitle='I (W m-2 sr-1 Hz-1)'
    oplot, freq, pf
    print, int_tabulated(freq, pf, /d) * !dpi * 4d * !dpi * 6.96d8^2
    ;   3.8260749e+26

    ; check the ratio
    irat = fsun / pf
    plot, freq, irat, /xlog, $
      xtitle='Frequency (Hz)', ytitle='Ikurucz / Iplanck'

    ; how much difference does log interpolation make?
    kinten2 = kurucz_inten_interp(inten, grav, temp, wgrav, wtemp, /log)
    afd = abs( (kinten-kinten2) / ((kinten+kinten2) / 2) )
    plot, wave, afd[where(kinten gt max(kinten) * 1e-2)], /xlog, /ylog

    For the sun, the fractional difference is 2e-3 or better for
    brightnesses greater than 1% of the peak, except for one spike
    to 5e-3.  Differences longward of 0.1 um are below 5e-5 and
    fall rapidly with increasing wavelength.

 MODIFICATION HISTORY:
     Written by:    Joseph Harrington, Cornell.  2006-02-02
            jh@oobleck.astro.cornell.edu
    2006-02-03 jh    Swapped wgrav and wtemp in calling order to
            match grav and temp.  Updated example for
            change of kurucz_flux_read to MKS.
    2006-02-10 jh    Fixed example, added kf keyword, changed name
            to kurucz_inten_interp.
    2006-02-11 jh    Header tweak.
    2006-02-12 jh    Converted variable names flux -> inten.
            Changed messages to printed warnings, since
            they make Monte Carlo trials bomb.
    2008-08-11 kevin
            Converted to Python
"""

def interp(inten, grav, temp, wgrav, wtemp, log=None, kf=None):

   import numpy as np
   import scipy.interpolate as interpolate

   dim    = inten.shape
   nwav   = dim[1]
   iinten = np.zeros(nwav)
   
   # FINDME: HARDWIRED!
   ilograv = 0
   ihigrav = 10 # the full range of gravities
   ilotemp = 0
   ihitemp = 11 # the temp where Kurucz starts skipping gravities
   ngrav = ihigrav - ilograv + 1
   ntemp = ihitemp - ilotemp + 1
   
   kg = grav[ilograv:ihigrav+1]
   it = np.arange(ntemp) * ngrav + ilotemp
   kt = temp[it]

   
   # FINDME: this only works if the low indices are both 0 and we're not
   # at temps where Kurucz starts skipping gravities...
   if kg[ihigrav] < wgrav:   
          print('kurucz_inten.interp: wgrav ' + str(wgrav) + ' higher than ' + str(kg[ihigrav]) + ', see FINDME in code to fix.')
   if kt[ihitemp] < wtemp:   
          print('kurucz_inten.interp: wtemp ' + str(wtemp) + ' higher than ' + str(kt[ihitemp]) + ', see FINDME in code to fix.')
   
   kf = np.reshape(inten[0:ngrav * ntemp,:], (ntemp, ngrav, nwav))
   
   # calculate the indices of the parameters we want
   temp, igrav = np.mgrid[kt[0]:kt[-1]+1:250,kg[0]:kg[-1]+0.1:0.5]
   
   if (log is not None):   
        kf = np.log(kf)
   
   for i in range(nwav):
        tck = interpolate.bisplrep(itemp, igrav, kf[:,:,i], kx=3, ky=3)
        iinten[i] = interpolate.bisplev(wtemp, wgrav, tck)
   
   if (log is not None):   
        iinten = exp(iinten)
   
   return iinten


"""
 NAME:
    READ

 PURPOSE:
    This function reads a file of stellar spectral intensity
    models from Bob Kurucz (Harvard).

 CATEGORY:
    Astronomy

 CALLING SEQUENCE:

    Result = READ(Filename)

 INPUTS:
    Filename:    Name of model file.  These come from
            http://kurucz.harvard.edu/grids.html

 KEYWORD PARAMETERS:
    FREQ:    If set, reverse first dimension of model grid and
        return frequencies, not wavelengths, in Wave.

 OUTPUTS:
    This function returns an [nwavl, nmod] array of model
        brightnesses.  The brightnesses in the file have been
        multiplied by 4, since they are Eddington fluxes
        (W m-2 sr-1 Hz-1).

    Wave:    [nwavl] array of wavelengths at which the intensities
        are calculated (meters).  Array of frequencies if
        /FREQ is set (Hz).

    Grav:    [nmod] array of the base-10 log of stellar surface
        gravities for the models (cm s-2).  NOTE: These units
        are so standard that we're leaving them as CGS!

    Temp:    [nmod] array of temperatures for the models (K).

    Nainten:    [nwavl, nmod] array of model brightnesses
        WITHOUT LINE ABSORPTION.  The intensities in the file
        have been multiplied by 4, since they are Eddington
        fluxes (W m-2 sr-1 Hz-1).

    Head:    [nmod] array of one-line headers for the models.

 EXAMPLE:
    #IDL EXAMPLE, NEEDS TO BE CONVERTED TO PYTHON
    ; read the model file

    inten = kurucz_inten_read('fp00ak2odfnew.pck', $
                 wave, grav, temp, nainten, head)

    ; plot the intensities vs. frequency in Hz
    nmod  = n_elements(head)
    loadct, 3
    waittime = 0.05
    c    = 2.99792458d8 ; speed of light in m/s
    freq = c / wave
    .run
    for i = 0, nmod - 1, 1 do begin
      plot,  freq,   inten[*, i], $
            xr=[1e14, 1e17], yr=[1e-25, 1e0], /xlog, /ylog, $
            title =   strtrim(i,       2) + '  ' $
                    + strtrim(temp[i], 2) + '  ' $
                + strtrim(grav[i], 2)
      oplot, freq, nainten[*, i], color = 180
      oplot, freq, planckfreq(freq, temp[i]) / !dpi * 1d-7 * 1d4, color=120
      wait, waittime
    endfor
    end

    ; estimate the luminosity of the sun (~4e26 W)
    tsun = 5770d
    gsun = 4.44d
    isun = (where(temp gt tsun and grav gt gsun))[0]
    tm = temp[isun]
    gm = grav[isun]
    fsun = inten[*,isun]
    freqp = reverse(freq)
    fsunp = reverse(fsun)
    ; The first !dpi converts from brightness to flux.
    print, int_tabulated(freqp, fsunp, /double) * !dpi $
                * 4 * !dpi * 6.96d8^2
    ;   4.4885570e+26
    ; The model is for a star with t=6000 and log g=4.5, so expect
    ; more than 4e26 W.

 MODIFICATION HISTORY:
     Written by:    Joseph Harrington, Cornell.  2006 Feb 1
            jh@oobleck.astro.cornell.edu
    Based on code provided by Drake Deming.
    2006 Feb 2 jh    Added solar luminosity example, fixed header,
            swapped order of grav and temp to match order
            in grid of models.
    2006 Feb 3 jh    Convert to MKS, update example.
    2006-02-10 jh    Rename to kurucz_inten_read, update header.
    2008-08-11    Kevin Stevenson, UCF
            kbstvenson@gmail.com
            Converted to Python
"""

def read(filename, freq=None):

   import numpy as np

   # read file into memory
   f = open(filename,'r')
   text = f.read()
   f.close()
   text = text.replace('\r','\n')
   filetxt = text.split('\n')

   # get, parse, and count header lines
   # Record effective temperatures and gravities
   head      = []
   temp      = np.zeros(len(filetxt))
   grav      = np.zeros(len(filetxt))-1
   header    = np.zeros(len(filetxt), int)
   startwave = 0
   for i in range(len(filetxt)):
        if filetxt[i].startswith("TEFF") == True:
          head.append(filetxt[i])
          temp[i]   = float(filetxt[i][ 5:12])
          grav[i]   = float(filetxt[i][22:29])
          header[i] = i
        elif filetxt[i].endswith("END") == True:
          startwave = i + 1

   temp   = temp  [np.where(temp   !=  0)]
   grav   = grav  [np.where(grav   != -1)]
   header = header[np.where(header !=  0)]
   nmod   = header.size
   nline  = (header[2] - header[1] - 1) / 2  #omit the header line itself

   # read and count wavelengths
   wave = np.zeros(int(header[0]*len(filetxt[startwave])/10))
   k = 0
   string = ''.join(filetxt[startwave:header[0]])
   for j in range(0,len(string),10):
    wave[k] = float(string[j:j+10])
    k += 1

   wave = wave[np.where(wave != 0)] * 1e-9  # convert nm to meters
   nwavl = wave.size

   # allocate memory for models
   inten   = np.zeros((nmod, nwavl))
   nainten = np.zeros((nmod, nwavl))

   #LOOP OVER MODELS
   for i in range(0, nmod):
       k = 0
       string1 = ''.join(filetxt[int(header[i]+1      ):int(header[i]+nline+1 ) ])
       string2 = ''.join(filetxt[int(header[i]+nline+1):int(header[i]+2*nline+1)])
       for j in range(0,len(string1),10):
           inten[i,k]   = float(string1[j:j+10])
           nainten[i,k] = float(string2[j:j+10])
           k += 1

   # convert Eddington fluxes to brightnesses, CGS (erg cm-2) -> MKS (J m-2)
   inten   *= 4. * 1e-3
   nainten *= 4. * 1e-3

   # convert to frequency if desired
   if (freq is not None):   
      c = 2.99792458e8           # speed of light in m/s
      wave = np.flipud(c / wave)
      inten = np.fliplr(inten)
      nainten = np.fliplr(nainten)
   
   return inten, wave, grav, temp, nainten, head

