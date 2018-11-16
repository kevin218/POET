# $Author: kevin $
# $Revision: 556 $
# $Date: 2011-09-13 11:32:07 -0400 (Tue, 13 Sep 2011) $
# $HeadURL: file:///home/esp01/svn/code/python/pipeline/trunk/p8tables.py $
# $Id: p8tables.py 556 2011-09-13 15:32:07Z kevin $

"""
 MODIFICATION HISTORY:
    Written by:	Kevin Stevenson, UCF  	2008-07-02
                kevin218@knights.ucf.edu
    Finished initial version:   kevin   2008-09-08
    Updated for multi events:   kevin   2009-11-01
    Added ip interpolation:     kevin   2010-06-28
    Added bjdutc and bjdtdb     kevin   2011-03-21
"""

from __future__ import print_function

def tables(event, fitnum, printout):

    import numpy as np
    import time
    import dec2sexa

    fit = event.fit[fitnum]

    print('MARK: ' + time.ctime() + ' : Print tables', file=printout)
    
    print('All interesting model parameters:', file=printout)
    print('%34s & %13s & %11s & %11s & %11s & %11s & %9s' % 
         ('Parameter', 'Value', 'LowError', 'HiError', 'TypErr', 'GaussErr', 'S/N'), file=printout)
    for i in range(len(fit.ifreepars)):
        print('%34s & %13.6f & %11.6f & %11.6f & %11.6f & %11.6f & %9.3f' % 
        (fit.parname[fit.ifreepars[i]], fit.bestp[fit.ifreepars[i]], fit.twosidederr[fit.ifreepars[i],0],
        fit.twosidederr[fit.ifreepars[i],1], fit.typerr[fit.ifreepars[i]], fit.medianp[fit.ifreepars[i],1],
        abs(fit.bestp[fit.ifreepars[i]] / fit.typerr[fit.ifreepars[i]])), file=printout)
    print('Original and rescaled reduced Chi^2 values:', file=printout)
    print(fit.oldredchisq, fit.redchisq, file=printout)
   
    # print wall-clock-format eclipse times
    if hasattr(fit.i, 'midpt'):
        print('Eclipse, h:mm:ss', file=printout)
        print('Deviation from ' + str(event.meanphase) + ' phase: ', \
            dec2sexa.dec2sexa1((fit.bestp[fit.i.midpt] - event.meanphase) * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.midpt,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.midpt,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.width] + ':         ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.width] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.width,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.width,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.t12] + ':             ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.t12] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.t12,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.t12,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.t34] + ':              ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.t34] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.t34,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.t34,1] * event.period * 24., 0), file=printout)
    if hasattr(fit.i, 'midpt2'):
        print('Eclipse, h:mm:ss', file=printout)
        print('Deviation from ' + str(event.meanphase) + ' phase: ', \
            dec2sexa.dec2sexa1((fit.bestp[fit.i.midpt2] - event.meanphase) * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.midpt2,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.midpt2,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.width2] + ':         ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.width2] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.width2,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.width2,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.t122] + ':             ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.t122] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.t122,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.t122,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.t342] + ':              ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.t342] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.t342,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.t342,1] * event.period * 24., 0), file=printout)
    if hasattr(fit.i, 'midpt3'):
        print('Eclipse, h:mm:ss', file=printout)
        print('Deviation from ' + str(event.meanphase) + ' phase: ', \
            dec2sexa.dec2sexa1((fit.bestp[fit.i.midpt3] - event.meanphase) * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.midpt3,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.midpt3,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.width3] + ':         ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.width3] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.width3,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.width3,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.t123] + ':             ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.t123] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.t123,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.t123,1] * event.period * 24., 0), file=printout)
        print(fit.parname[fit.i.t343] + ':              ', \
            dec2sexa.dec2sexa1(fit.bestp[fit.i.t343] * event.period * 24., 0), ' - ', \
            dec2sexa.dec2sexa1(abs(fit.twosidederr[fit.i.t343,0]) * event.period * 24., 0), ' + ', \
            dec2sexa.dec2sexa1(fit.twosidederr[fit.i.t343,1] * event.period * 24., 0), file=printout)

    # BJD of eclipse, ephemeris
    if hasattr(fit.i, 'midpt'):
        #julobsdate = event.j2kjd + fit.time / 86400.
        ctrphase = fit.bestp[fit.i.midpt]            # center phase
        ctrerr   = fit.twosidederr[fit.i.midpt,:]    # center phase error
        #CHOOSE eps SUCH THAT ctrind IS ONLY A FEW VALUES
        #eps      = (fit.time[1] - fit.time[0])/2e5
        foo      = np.where(fit.abscissauc >= ctrphase)[0][0]
        ctrind   = [foo-1, foo]
        print('Center phase and 2-sided error:', file=printout)
        print(ctrphase, ctrerr, file=printout)
        print('Indices of closest frames:', file=printout)
        print(ctrind, file=printout)
        print('Phases of closest frames:', file=printout)
        print(fit.abscissauc[ctrind], file=printout)
        try:
            print('BJD_UTC of closest frames:', file=printout)
            print(fit.bjdutcuc[ctrind], file=printout)
        except:
            print('Could not compute BJD_UTC of closest frames.', file=printout)
        try:
            print('BJD_TDB of closest frames:', file=printout)
            print(fit.bjdtdbuc[ctrind], file=printout)
        except:
            print('Could not compute BJD_TDB of closest frames.', file=printout)

    #BRIGHTNESS TEMPERATURES WITH ERROR
    if hasattr(fit, 'tbm'):
        print('Brightness Temperature = ' + str(round(fit.tbm,1)) + ' +/- ' + str(round(fit.tbsd,1)) + ' K', file=printout)

    #Recap of important parameters
    print('%9s %8s %8s %8s %8s %10s %10s %10s %10s' % 
         ('MODEL',  'SDNR', 'BIC', 'MIDPT', 'ERROR', 'DEPTH (%)', 'ERROR (%)', 'TEMP (K)', 'ERROR (K)'), file=printout)
    if hasattr(fit.i, 'midpt'):
        print('%9s %8.6f %8.1f %8.4f %8.5f %10.4f %10.4f %10.0f %10.0f' % 
             (fit.saveext, fit.sdnr, fit.bic, fit.bestp[fit.i.midpt], fit.typerr[fit.i.midpt],
              fit.bestp[fit.i.depth]*100, fit.typerr[fit.i.depth]*100, fit.tbm, fit.tbsd), file=printout)
    if hasattr(fit.i, 'midpt2'):
        print('%9s %8.6f %8.1f %8.4f %8.5f %10.4f %10.4f %10.0f %10.0f' % 
             (fit.saveext, fit.sdnr, fit.bic, fit.bestp[fit.i.midpt2], fit.typerr[fit.i.midpt2],
              fit.bestp[fit.i.depth2]*100, fit.typerr[fit.i.depth2]*100, fit.tbm2, fit.tbsd2), file=printout)
    if hasattr(fit.i, 'midpt3'):
        print('%9s %8.6f %8.1f %8.4f %8.5f %10.4f %10.4f %10.0f %10.0f' % 
             (fit.saveext, fit.sdnr, fit.bic, fit.bestp[fit.i.midpt3], fit.typerr[fit.i.midpt3],
              fit.bestp[fit.i.depth3]*100, fit.typerr[fit.i.depth3]*100, fit.tbm3, fit.tbsd2), file=printout)

    return

