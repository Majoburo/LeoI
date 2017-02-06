#!/usr/bin/env python
# srebin.py v0.2 by Majo
# What's new:
# Added loglogspl shifting.
#
# srebin.py v0.1 by Maximilian Fabricius 
# performs spectral rebinning

from astropy.io import fits as pyfits
from math import *
from scipy import *
from scipy import interpolate
import numpy as np
import sys


cL = 299792.


def loglogSpl(wl, fluxv, step=0., start=0., stop=0., zshift=0.):
    """log_rebin(wl, flux, start, stop=0, step)
    Performs rebinning an optional red/blue shifting
    wl an array of input wavelengths in log angstrom
    flux and array of fluxes in arbitrary untis
    start start wavelength in A
    stop (optional stop wavelength in A)
    step in (km/s)/c
    zshift (optional shift that the output will be shifted to in km/s)
    """
    warned = False

    # make sure the arrays are sorted in ascending order    
    ii = argsort(wl)
    wl = wl[ii]
    fluxv = fluxv[ii]

    if step == 0.:
        step = wl[1]-wl[0]
    if start == 0. or start < wl[0]:
        start = wl[0]
    if stop == 0. or stop > wl[-1]:
        stop = wl[-1]

    #print "data cover wl range: ", wl[0], wl[-1]
    #print "start: ", start
    #print "stop:  ", stop
    #print "step:  ", step
    tck = interpolate.splrep(wl,fluxv,s=0)

    #output
    out = []

    nmax = (stop + start + step/2.)/step
    beta = zshift/cL

    aalam  = (arange(nmax+1)*step + start)+beta
    aalam1 = aalam + (-step/2.)
    aalam2 = aalam + (+step/2.)
    n = 0
    while True:
        # calculate current window size    
        if aalam2[n] > stop:
            break
        if aalam1[n] < wl[0]:
            out.append(0)
            n+=1
            continue

        flx = interpolate.splint(aalam1[n],aalam2[n],tck)/(aalam2[n]-aalam1[n])
        out.append(flx)
        n+=1
    out_wls = arange(len(out))*step + start
    if len(fluxv)-len(out) < 0:
        out = out[:len(fluxv)-len(out)]
    else:
        out = np.pad(out,(0,len(fluxv)-len(out)),'edge')
    return wl, array(out)

def linlinSpl(wl, fluxv, step=0., start=0., stop=0., zshift=0.):
    """log_rebin(wl, flux, start, stop=0, step)
    Performs rebinning an optional red/blue shifting
    wl an array of input wavelengths in angstrom
    flux and array of fluxes in arbitrary untis
    start start wavelength in A
    stop (optional stop wavelength in A)
    step in (km/s)/c  
    zshift (optional shift that the output will be shifted to in km/s)
    """
    warned = False

    # make sure the arrays are sorted in ascending order    
    ii = argsort(wl)
    wl = wl[ii]
    fluxv = fluxv[ii]
    
    if step == 0.:
        step = wl[1]-wl[0]
    if start == 0. or start < wl[0]:
        start = wl[0]
    if stop == 0. or stop > wl[-1]:
        stop = wl[-1]
        
    tck = interpolate.splrep(wl,fluxv,s=0)
    
    #output
    out = []
    
    nmax = (stop + start + step/2.)/step
    #aalam  = (arange(nmax+1)*step + start)*(1. - zshift/cL)
    beta = zshift/cL

    aalam  = (arange(nmax+1)*step + start)*sqrt((1.-beta)/(1.+beta))
    #print start*(sqrt((1.-beta)/(1.+beta)) - 1)/step
    aalam1 = aalam + (-step/2.)
    aalam2 = aalam + (+step/2.)
    n = 0
    while True:
        # calculate current window size    
        if aalam2[n] > stop:
            break
        if aalam1[n] < wl[0]:
            out.append(0)
            n+=1
            continue
        
        flx = interpolate.splint(aalam1[n],aalam2[n],tck)/(aalam2[n]-aalam1[n])
        out.append(flx)
        n+=1
    out_wls = (arange(len(out))*step + start)    

    return out_wls, array(out)

def linlogSpl(wl, fluxv, start=0., step=0., stop=0.):
    """
    Performs logarithmic rebinning of arbitrarely sampled input spectra though spline integration.
    wl an array of input wavelengths in angstrom
    flux and array of fluxes in arbitrary untis
    start optional start wavelength in A, if not giving the starting wavelength wil be used
    stop (optional stop wavelength in A)
    step in (km/s)/c, if not given the (wl[1]-wl[0]) / wl[0] will be used
    """
    
    warned = False

    # make sure the array are sorted in ascending order    
    ii = argsort(wl)
    wl = wl[ii]
    fluxv = fluxv[ii]
    
    if stop == 0. or stop > wl[-1]:
        stop = wl[-1]
    if start == 0.:
        start = wl[0]
    if step == 0.:
        step = (wl[1]-wl[0])/wl[0]
        
    print "data cover wl range: ", wl[0], wl[-1]
    print "start: ", start
    print "stop:  ", stop
    print "step:  ", step
    tck = interpolate.splrep(wl,fluxv,s=0)

    #output
    out = []

    n = -1.
    while True:
        # calculate current window size    
        n += 1.
        alam  = start * exp(n*step)
        alam1 = alam  * exp(-step/2.)
        alam2 = alam  * exp(step/2.)
        #print "alam = %f + %f - %f " % (alam, alam1, alam2)
        #end while-loop if end of input range is reached
        if alam2 > stop:
            break

        if alam1 < wl[0]:
            out.append(0)
            continue
    
        flx = interpolate.splint(alam1,alam2,tck)
        out.append( flx/(alam2-alam1) )

    out_wls = arange(len(out))*step + log(start)    

    return out_wls, array(out)

def loglinSpl(wl, fluxv, start=0., step=0., stop=0.):
    """
    Performs logarithmic rebinning of arbitrarely sampled input spectra though spline integration.
    wl an array of input wavelengths in angstrom
    flux and array of fluxes in arbitrary untis
    start optional start wavelength in ln(lambda), if not given the starting wavelength will wl[0] be used
    stop (optional stop wavelength in ln(lambda) )
    step in lambda, if not given the (wl[1]-wl[0]) * wl[0] will be used
    """
    
    warned = False

    # make sure the array are sorted in ascending order    
    ii = argsort(wl)
    wl = wl[ii]
    fluxv = fluxv[ii]
    
    if stop == 0. or stop > wl[-1]:
        stop = wl[-1]
    if start == 0.:
        start = wl[0]
    if step == 0.:
        step = (wl[1]-wl[0])*exp(wl[0])
    #print "start: ", start
    #print "step: ", step
    #print "stop: ", stop
    tck = interpolate.splrep(wl,fluxv,s=0)
    
    #output
    out = []

    n = -1.
    while True:
        # calculate current window size    
        n += 1.
        alam  = log(exp(start) + n*step)
        alam1 = log(exp(alam)  + -step/2.)
        alam2 = log(exp(alam)  + step/2.)
        #print "alam = %f + %f - %f " % (alam, alam1, alam2)
        #end while-loop if end of input range is reached
        if alam2 > stop:
            break

        if alam1 < wl[0]:
            out.append(0)
            continue
        
        flx = interpolate.splint(alam1,alam2,tck)
        out.append(flx/(alam2-alam1))
    out_wls = arange(len(out))*step + exp(start)    

    return out_wls, array(out)

def linlog(wl, fluxv, start, step, stop=0.):
    """log_rebin(wl, flux, start, stop=0, step)
    Performs logarithmic rebinning of arbitrarely sampled input spectra.
    wl an array of input wavelengths in angstrom
    flux and array of fluxes in arbitrary untis
    start start wavelength in A
    stop (optional stop wavelength in A)
    step in km/s   
    function function for rebinning. Note for lin to log binning, "exp" has to be given, for 
    log to lin binning "log" has to be given!
    """
    
    warned = False

    # make sure the array are sorted in ascending order    
    ii = argsort(wl)
    wl = wl[ii]
    fluxv = fluxv[ii]
    
    if stop == 0. or stop > wl[-1]:
        stop = wl[-1]
    
    n1 = interpolate.interp1d( wl, arange(wl.shape[0]), 'linear')
    
    #output
    out = []

    n = -1.
    while True:
        # calculate current window size    
        n += 1.
        alam  = start * exp(n*step)
        alam1 = alam  * exp(-step/2.)
        alam2 = alam  * exp(step/2.)
        #print "alam = %f + %f - %f " % (alam, alam1, alam2)
        #end while-loop if end of input range is reached
        if alam2 > stop:
            break

        if alam1 < wl[0]:
            out.append(0)
            continue
        #tanslate wavelengths into pixel indices through linear interpolation
        x1  = n1(alam1)[0]
        x2  = n1(alam2)[0]
        mx1 = floor(x1)
        mx2 = floor(x2)
        if (x1 <= 0.) :
            mx1 = 0.
        if (x2 >= len(wl)):
            mx2 = len(wl)-2.

        # first integrate over all pixel which are fully contained in the
        #  current interval, integrate piecewise from the 
        #  center of on pixel to the next (therfore (pixel_i + pixel_(i+1))/2.)
        #
        flx=0.
        for ipix in range(mx1+1,mx2-1 + 1):
            flx=flx+( fluxv[ipix] + fluxv[ipix+1] )/2.
        
        # Now take care of the edges of the current interval using the
        # trapezium rule
        if (x1 < 0.):
            x1 = 0.

        if (x2 > len(wl)-1 ):
             x2 = len(wl)-1
        
        if mx1 < mx2:
            a1 = fluxv[mx1+1.]-fluxv[mx1]
            b1 = fluxv[mx1]
            f_prime1 = a1*(x1-mx1) + b1
            flx1 = (mx1+1.-x1) * (f_prime1 + fluxv[mx1+1])/2.
            
            a2 = fluxv[mx2+1]-fluxv[mx2]
            b2 = fluxv[mx2]
            f_prime2 = a2*(x2 - mx2) + b2
            flx2 = (x2-mx2) * ( f_prime2 + fluxv[mx2] )/2.
        elif mx1 == mx2:
            if not warned:
                print "Warning: Such a small step width leads to sub pixel bins."
                warned = True            
            a1 = fluxv[mx1+1.]-fluxv[mx1]
            b1 = fluxv[mx1]
            f_prime1 = a1*(x1-mx1) + b1
            f_prime2 = a1*(x2-mx1) + b1
            flx1 = (x2-x1) * (f_prime1+f_prime2)/2.
            flx2 = 0.
        else:
            print "Error mx1 > mx2, should never happen!"

        out.append(flx + flx1 + flx2)
    out_wls = arange(len(out))*step + log(start)    

    return out_wls, array(out)

def linlin(x, fluxv, start, step, stop=0.):
    """lin_rebin(x, flux, start, stop=0, step)
    Performs linear rebinning of arbitrarely sampled input spectra.
    x an array of input wavelengths in angstrom
    flux and array of fluxes in arbitrary untis
    start start wavelength in A
    stop (optional stop wavelength in A)
    step in km/s   
    function function for rebinning. Note for lin to log binning, "exp" has to be given, for 
    log to lin binning "log" has to be given!
    """
    
    warned = False

    # make sure the array are sorted in ascending order    
    ii = argsort(x)
    x = x[ii]
    fluxv = fluxv[ii]
    
    if stop == 0. or stop > x[-1]:
        stop = x[-1]

    # prepare pix(wavelength) interpolation class    
    n1 = interpolate.interp1d( x, arange(x.shape[0]), 'linear', bounds_error=False)
    
    #output
    out = []
    
    # do the interpolation here, just one, this speeds the whole process up
    nmax = (stop - -start - step/2.)/step
    aalam  = arange(nmax+1)*step + start
    aalam1 = aalam + (-step/2.)
    aalam2 = aalam + (+step/2.)
    
    xx1  = n1(aalam1)
    xx2  = n1(aalam2)

    n = -1.
    while True:
        # calculate current window size    
        n += 1.
        #print "alam = %f + %f - %f " % (alam, alam1, alam2)
        #end while-loop if end of input range is reached
        if aalam2[n] > stop:
            break

        if aalam1[n] < x[0]:
            out.append(0)
            continue
        #tanslate wavelengths into pixel indices through linear interpolation
        x1  = xx1[n]
        x2  = xx2[n]
        mx1 = floor(x1)
        mx2 = floor(x2)
        if (x1 <= 0.) :
            mx1 = 0.
        if (x2 >= len(x)):
            mx2 = len(x)-2.
        # first integrate over all pixel which are fully contained in the
        #  current interval, integrate piecewise from the 
        #  center of on pixel to the next (therefor (pixel_i + pixel_(i+1))/2.)
        #
        flx=0.
        for ipix in range(mx1+1,mx2-1 + 1):
            flx=flx+( fluxv[ipix] + fluxv[ipix+1] )/2.
        
        # Now take care of the edges of the current interval using the
        # trapezium rule
        if (x1 < 0.):
            x1 = 0.

        if (x2 > len(x)-1 ):
             x2 = len(x)-1
            
        if mx1 < mx2:
            a1 = fluxv[mx1+1.]-fluxv[mx1]
            b1 = fluxv[mx1]
            f_prime1 = a1*(x1-mx1) + b1
            flx1 = (mx1+1.-x1) * (f_prime1 + fluxv[mx1+1])/2.
            
            a2 = fluxv[mx2+1]-fluxv[mx2]
            b2 = fluxv[mx2]
            f_prime2 = a2*(x2 - mx2) + b2
            flx2 = (x2-mx2) * ( f_prime2 + fluxv[mx2] )/2.
        elif mx1 == mx2:
            if not warned:
                print "Warning: Such a small step width leads to sub pixel bins."
                warned = True            
            a1 = fluxv[mx1+1.]-fluxv[mx1]
            b1 = fluxv[mx1]
            f_prime1 = a1*(x1-mx1) + b1
            f_prime2 = a1*(x2-mx1) + b1
            flx1 = (x2-x1) * (f_prime1+f_prime2)/2.
            flx2 = 0.
        else:
            print "Error mx1 > mx2, should never happen!"

        out.append(flx + flx1 + flx2)
    out_xs = arange(len(out))*step + (start)    

    return out_xs, array(out)


