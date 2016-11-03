import numpy as np
from astropy.io import fits
import re
import math
from astropy import wcs
import matplotlib.pyplot as plt
from copy import copy as cp

"""
My list of tools:

updatefits
degtopix
getdata
kappasigma
maskemission_ln
regtodeg
regtopix

"""
def updatefits(hdu, newdata, newrange):
    """Updates fits files data and header.

    Parameters
    ----------
    hdu : fits file
    newdata: 2darray
    newrange: 1darray

    Returns
    -------
    h
        Updated fits file.
    """
    hdu[0].data = newdata
    h = hdu[0].header
    h['CRVAL1'] = newrange[0]
    h['CDELT1'] = newrange[1] - newrange[0]
    h['CD1_1'] = h['CDELT1']
    h['NAXIS1'] = len(newrange)
    hdu[0].header = h

    return  hdu

def degtopix(hdu,degrees):
    """Converts degrees to pixels

    Parameters
    ----------
    hdu : fits file
    degrees : float

    Returns
    -------
    float
        pixels corresponding to the given degrees
    """

    mat = np.round([[hdu[0].header['CD1_1'],hdu[0].header['CD1_2']],
                   [hdu[0].header['CD2_1'],hdu[0].header['CD2_2']]],10)
    CDELT1 = pow(abs(np.linalg.det(mat)),0.5) #Coordinate increment per pixel in DEGREES/PIXEL

    return degrees/CDELT1

def pixtodeg(hdu, pixels):
    """Converts pixels to degrees

    Parameters
    ----------
    hdu : fits file
    pixels : float

    Returns
    -------
    float
        degrees corresponding to the given pixels
    """

    mat = np.round([[hdu[0].header['CD1_1'],hdu[0].header['CD1_2']],
                   [hdu[0].header['CD2_1'],hdu[0].header['CD2_2']]],10)
    CDELT1 = pow(abs(np.linalg.det(mat)),0.5) #Coordinate increment per pixel in DEGREES/PIXEL

    return CDELT1*pixels


def getdata(filename, norm = False):

    """Gives wavelenght range and data for a given set of spectra

    Parameters
    ----------
    filename : string

    Options
    -------
    norm: bool
        Option for having the spectra normalized. Useful for cross - correlations.

    Returns
    -------
    data : 2darray
        fits file data

    datarange : 1darray
        wavelenght range
    """

    hdu = fits.open(filename)
    data = hdu[0].data
    if norm:
        if len(data.shape) == 1:
            data = data/np.mean(data, dtype=np.float64)
        else:
            for i, spectra in enumerate(data):
                data[i] = spectra/np.mean(spectra, dtype=np.float64)
    h = hdu[0].header
    datarange = h['CRVAL1'] + h['CDELT1']*np.arange(h['NAXIS1'])
    return data, datarange

def kappasigma(data, range, kappa = 2., iterations = 5):

    """Takes the continum of spectra. Does so by iteratively fitting a third order polynomial
    to the spectra, subtracting it, and then masking those values that deviate from the norm
    by kappa*sigma times.

    Parameters
    ----------
    data : 2darray
    range : 1darray
        Has to be a linear scale.
    kappa: float
        Value by .

    Returns
    -------
    data : ndarray
        fits file data

    datarange : 1darray
        wavelenght range
    """
    #Below is a piecewise polyfit. I fear this will just lower the correlations, so I don't use it.
     
    window = 1000. #roughly dividing the spectra by 4
    weight = np.ones(len(range))
    #weight[0]= 100000
    #weight[-1]= 100000
    #j = math.ceil(len(range)/window)
    pfitdata = np.poly1d(np.polyfit(range, data, 5, w = weight))
    data = data - pfitdata(range)
    #Taking top 10 %
    ordata = cp(data)
    #weight = np.ones(window)
    for j in xrange(int(math.ceil(len(range)/window))):
        for i in xrange(iterations):
            if window*(j+1) < len(range):
                window_range = range[window*j:window*(j+1)]
                window_data = data[window*j:window*(j+1)]
                if i == 0:
                    weight = np.zeros(len(window_data))
                weight[np.argsort(window_data)[:len(window_data)*8./10.]] = 1.
                sigma = np.std(window_data)
                weight *= (abs(window_data) < kappa*sigma)
                pfitdata = np.poly1d(np.polyfit(window_range, window_data, 3, w = weight))
                window_data = window_data - pfitdata(window_range)
                data[window*j:window*(j+1)]= window_data
                #weight = np.ones(len(window_range))
            else:
                window_range = range[window*j:]
                window_data = data[window*j:]
                if i == 0:
                    weight = np.zeros(len(window_data))
                weight[np.argsort(window_data)[:len(window_data)*8./10]] = 1
                sigma = np.std(window_data)
                weight *= (abs(window_data) < kappa*sigma)
                pfitdata = np.poly1d(np.polyfit(window_range, window_data, 3, w = weight))
                window_data = window_data - pfitdata(window_range)
                data[window*j:] = window_data
        #plt.plot(window_range,pfitdata(window_range))
        #plt.plot(window_range,kappa*sigma*np.ones(len(window_range)))
                #weight = np.ones(len(window_range))
    #plt.plot(range, data,"--")
    #plt.plot(range,ordata,'-')
    #plt.show()
    ''' 
    weight = np.ones(len(range))
    polyrange = np.arange(len(range))
    for i in xrange(iterations):

        pfitdata = np.poly1d(np.polyfit(polyrange, data, (i+1)*3, w = weight))
        data = data - pfitdata(polyrange)
        sigma = np.std(data)
        weight = (abs(data) < kappa*sigma)
    '''
    return data

def maskemission_ln(range, abslines):

    """Makes a rough cut of 20 angstorms for the given set of wavelenghts
    so they dont distort the polynomial


    Parameters
    ----------
    range : 1darray
    abslines: 1darray
        Array of absorption lines in Angstorms

    Returns
    -------
    data : ndarray
        fits file data

    datarange : 1darray
        wavelenght range
    """
    weight = np.ones(len(range))
    for line in abslines:
        weight *= (range > np.log(line + 10)) | (range < np.log(line - 10))
    return weight

def regtodeg(regionsfile = "vwinfo/fibers_numbered.reg"):

    """Tiny regular expressions code for extracting the
    fiber#, RA, DEC of vw fibers regions file.
    I used this to cross check that the contamination numbers of nearby stars were right.

    Parameters
    ----------
    regionsfile : .reg file
        VERY SPECIFIC REG FILE, WITH TXT FOR EVERY FIBER

    Returns
    -------
    data : 2darray
        [Fiber#, RA, DEC] per row

    """
    regions = []
    with open(regionsfile) as f:
        for line in f.readlines()[::2][4:]: #Read every two lines starting from the fourth 
            line = re.split("[,\(\)\{\}]",line)
            regions.append([int(line[4]),float(line[1]),float(line[2])])
    print regions
    return np.array(regions)

def regtopix(fitsfile, regionsfile = "vwinfo/fibers_numbered.reg"):

    """Get pixel positions of fibers for a given regions file on a fits image.

    Parameters
    ----------
    fitsfile: .fits file
        Image with RA DEC info in headers.

    regionsfile : .reg file
        VERY SPECIFIC REG FILE, WITH TXT FOR EVERY FIBER

    Returns
    -------
    data : 2darray
        [Fiber#, pixx, pixy] per row

    """

    fibers = regtodeg(regionsfile)
    with fits.open(fitsfile) as h:
        img_data = h[0].data
        w = wcs.WCS(h[0].header)
        pixcrd = w.wcs_world2pix(np.array(fibers[:,1:3],dtype='float64'),1)

    return np.column_stack([fibers[:,0],pixcrd])
