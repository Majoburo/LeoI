import numpy as np
from scipy.integrate import quad
from scipy.special import iv
from astropy.io import fits
#np.set_printoptions(precision=2)
import os
import sys
PATH = "/".join(os.getcwd().split("/")[:-1])
sys.path.insert(0,PATH)
import toolbox
import random
from scipy.signal import gaussian
import matplotlib.pyplot as plt
from astropy import wcs
import re

"""
DISCLAIMER:
    Do not trust anybody. Re-check everything. But remember, time is limited.
"""

#PARAMETERS
fluxnotdone = True
seeing = 1.5 #Degrees
fRdeg = 1.6 + 0.5 #Fiber radius in deg plus some error disk of radius 0.5
#Ellipse parameters
center = np.array([ 152.11436881,   12.304033  ])
centerpix = np.array([500.00,  286.00])
ab = 0.8 #a/b=0.08 but it was inclined 90 in KARL's image.
#I'll break the anulus from the center by which to classify the fibers into bins of 12 so its a multiple of 4 and 3.

##Volume under a 2d Gaussian fo Magnitude M and standard deviation sig
##intersected by a circle of radius fR at distance cD (center to center)
def integrand1(r,cD,fR,sig,M):
    return 2 * np.pi * r *  10 ** ((25.-M)/2.5) * iv(0, - cD * r/(sig*sig) ) * np.exp( - (cD * cD + r * r)/(2 * sig * sig))


def integrate(cD,fR,sig,M):
    return quad(integrand1, 0, fR, args = (cD, fR, sig, M))[0]


def binfibers(fitsfile,regionsfile):
    regions = []
    with open(regionsfile) as f:
        for line in f.readlines()[::2][4:]:
            line = re.split("[,\(\)\{\}]",line)
            regions.append([int(line[4]),float(line[1]),float(line[2])])
    regions = np.array(regions)
    regions = regions[np.argsort(regions, axis=0)[:,0]]

    with fits.open(fitsfile) as hdu:
        w = wcs.WCS(hdu[0].header)
        regpixcrd = w.wcs_world2pix(np.array(regions[:,1:3],dtype='float64'),1)

    furtherxy = (centerpix-regpixcrd)[np.argmax(np.sum((centerpix - regpixcrd)**2,axis=1))]
    b = np.sqrt(furtherxy[0]**2 + (furtherxy[1]/ab)**2)
    bins = 12
    fiberbin = np.zeros(len(regions))
    for i in xrange(bins):
        mask = ( (np.sqrt(   ((centerpix - regpixcrd)[:,0])**2
                        + ((centerpix - regpixcrd)[:,1]/ab)**2) <= (i+1)*b/bins)
                *(np.sqrt(   ((centerpix - regpixcrd)[:,0])**2
                        + ((centerpix - regpixcrd)[:,1]/ab)**2) > i*b/bins))
        fiberbin = fiberbin + mask*(i+1)
    np.savetxt("fiberbin.txt",fiberbin)
    return fiberbin


def getflux(fitsfile,catalogfile,regionsfile):
    catalog = np.loadtxt(catalogfile)
    regions = []
    with open(regionsfile) as f:
        for line in f.readlines()[::2][4:]:
            line = re.split("[,\(\)\{\}]",line)
            regions.append([int(line[4]),float(line[1]),float(line[2])])
    regions = np.array(regions)

    with fits.open(fitsfile) as hdu:
        w = wcs.WCS(hdu[1].header)
        catpixcrd = w.wcs_world2pix(np.array(catalog[:,2:4],dtype='float64'),1)
        regpixcrd = w.wcs_world2pix(np.array(regions[:,1:3],dtype='float64'),1)
        hdus = []
        hdus.append(hdu[1])
        pixsig = toolbox.degtopix(hdus, seeing/2.35/3600.) #Divided by 2.34 cause seeing is 2.35*sigma
        pixfR = toolbox.degtopix(hdus, fRdeg/3600.) #Plus some error disk of radius 0.5

    leophot = np.array(zip(catalog[:,7], catpixcrd[:,0], catpixcrd[:,1]))
    fiberpos = np.array(zip(regions[:,0], regpixcrd[:,0], regpixcrd[:,1]))
    fiberpos = fiberpos[np.argsort(fiberpos, axis=0)[:,0]]

    mask = ((leophot[:,2] < np.max(fiberpos[:,2]) + pixfR*2)*(leophot[:,2] > np.min(fiberpos[:,2]) - pixfR*2)*
            (leophot[:,1] < np.max(fiberpos[:,1]) + pixfR*2)*(leophot[:,1] > np.min(fiberpos[:,1]) - pixfR*2))

    fibers_flux = []
    for fiber in fiberpos:
        mask_nb_fib = (fiber[1] - leophot[mask][:,1])**2 + (fiber[2] - leophot[mask][:,2])**2 < pixfR**2
        vec_dist = (fiber[1:3] - leophot[:,1:3][mask][mask_nb_fib])
        cD = np.sqrt(vec_dist[:,0]**2 + vec_dist[:,1]**2)
        M = leophot[mask][mask_nb_fib][:,0]
        flux = []
        for j, dist in enumerate(cD):
            #print i,cD[j],integrate(cD=cD[j],fR = fR, sig = sig, M=M[j])/(2 * 10 ** ((25.-M[j])/2.5) * sig * sig * np.pi) * 100
            flux.append(integrate(dist,fR = pixfR, sig = pixsig, M=M[j]))
        fibers_flux.append([fiber[0],flux/np.sum(flux)*100.])

    return fibers_flux


fitsfile = "KARL/ds9.fits"
regionsfile = "KARL/LeoI_sky.reg"
#catalogfile = "KARL/LeoI_photometry.als"
fiberbin = binfibers(fitsfile,regionsfile)
#fibers_flux1 = getflux(fitsfile,catalogfile,regionsfile)
#np.savetxt("Karlmax.txt", np.max(fibers_flux1[:,1]))
#exit
#exit()
if fluxnotdone:

#FITS 1
    fitsfile = "hst_10520_02_acs_wfc_f435w/hst_10520_02_acs_wfc_f435w_drz.fits"
    catalogfile = "hst_10520_02_acs_wfc_f435w/hst_10520_02_acs_wfc_f435w_daophot_trm.cat"
    regionsfile = "hst_10520_02_acs_wfc_f435w/LeoI_sky.reg"
    fibers_flux1 = getflux(fitsfile,catalogfile,regionsfile)

#FITS 2

    fitsfile2 = "hst_12304_02_wfc3_uvis_f555w/hst_12304_02_wfc3_uvis_f555w_drz.fits"
    catalogfile2 = "hst_12304_02_wfc3_uvis_f555w/hst_12304_02_wfc3_uvis_f555w_daophot_trm.cat"
    regionsfile2 = "hst_12304_02_wfc3_uvis_f555w/LeoI_sky.reg"
    fibers_flux2 = getflux(fitsfile2,catalogfile2,regionsfile2)

#FITS 3

    fitsfile3 = "hst_12304_01_wfc3_uvis_f555w/hst_12304_01_wfc3_uvis_f555w_drz.fits"
    catalogfile3 = "hst_12304_01_wfc3_uvis_f555w/hst_12304_01_wfc3_uvis_f555w_daophot_trm.cat"
    regionsfile3 = "hst_12304_01_wfc3_uvis_f555w/LeoI_sky.reg"
    fibers_flux3 = getflux(fitsfile3,catalogfile3,regionsfile3)

    for fiberflux in fibers_flux1: #fibers_flux is formated: fibers_flux[0] = fiber#, fibers_flux[1] = fluxes
        fiber = int(fiberflux[0])
        if (len(fiberflux[1]) < len(fibers_flux2[fiber][1])): #comparing the amount of stars per fiber
            print "Fiber %d has too few stars, probably not fully covered by the catalog..."%fiber
            print "Using alternate catalog for this one..."
            fibers_flux1[fiber] = fibers_flux2[fiber]
        if (len(fiberflux[1]) < len(fibers_flux3[fiber][1])):
            print "Fiber %d has too few stars even in the second cat..."%fiber
            print "Using alternate catalog for this one..."
            fibers_flux1[fiber] = fibers_flux3[fiber]

    #np.save("fibersflux.txt",fibers_flux1)

else:
    fibers_flux1 = np.loadtxt("fibersflux.txt")


'''
    newgauss=np.zeros((len(flux),10000))
    maxgauss = []
    for f in xrange(1000):
        for i,w in enumerate(flux/np.sum(flux)):
            shift = int(random.gauss(0,100))
            newgauss[i] = np.take(w*gaussian(10000,270),np.arange(shift,10000+shift),mode = 'wrap')
        sumgauss = np.sum(newgauss,axis=0)
        maxgauss.append(np.argmax(sumgauss))
    fiberIMAGEstd.append(np.std(maxgauss))

    plt.plot(newgauss[0])
    plt.plot(newgauss[1])
    plt.plot(newgauss[2])
    plt.plot(sumgauss)
    plt.axvline(np.argmax(sumgauss), color='r')
    plt.show()
'''

corrfib = np.loadtxt("../compare/fibers.txt")

import time
start = time.time()

fiberIMAGEstd = []
maxflux = []
stdstd = []
iterations = int(sys.argv[2])
enhance = 100
guess_speed = float(sys.argv[1])*enhance
lenarray = 100*enhance
width = 27*(enhance/10)
#Doing normalized cross correlation (just as like with the data)
tmpmean = np.mean(gaussian(lenarray, width))
tmpstd = np.std(gaussian(lenarray, width))
tmpgauss = (gaussian(lenarray, width) - tmpmean)/tmpstd

for fibb in corrfib:
    flux = fibers_flux1[fibb[0].astype(int)][1]
    if len(flux)==0:
        print "Fiber %d has no cataloged stars inside or nearby..."%fibb[0]
        fiberIMAGEstd.append(0)
        maxflux.append(0)
    else:
        newgauss = np.zeros((len(flux), lenarray))
        maxgauss = []
        for f in xrange(iterations):
            for i,w in enumerate(flux/np.sum(flux)):
                shift = int(random.gauss(0, guess_speed))
                newgauss[i] = np.take(w*gaussian(lenarray, width), np.arange(shift, lenarray + shift),mode = 'wrap')
            sumgauss = np.sum(newgauss, axis = 0)
            sumgauss = (sumgauss - np.mean(sumgauss))/np.std(sumgauss)
            corr = np.correlate(sumgauss, tmpgauss, mode = 'same')[lenarray/4:lenarray*3/4]
            maxgauss.append(np.argmax(corr))
            #plt.plot(corr)
            #plt.show()
            #maxgauss.append(np.argmax(np.correlate(sumgauss, tmpgauss)[4000:6000]))


        fiberIMAGEstd.append(np.std(maxgauss))
        maxflux.append(np.max(flux))

fiberIMAGEstd = np.array(fiberIMAGEstd)
end = time.time()
print "This calculation took (in seconds):"
print(end - start)
fiberbin = fiberbin[(corrfib[:,0]).astype(int)]
np.savetxt("std%f_%d.txt"%(guess_speed/100,iterations),
            np.array(zip(corrfib[:,0],
                         corrfib[:,2],
                         fiberIMAGEstd/(10.),
                         maxflux,
                         fiberbin)), header = "fiber vel HSTstd HSTstdstd HSTmaxflux fiberbin")


'''
print corrfib
j=0
contamination = np.zeros((len(corrfib),3))
for k,fibb in enumerate(corrfib):
    if len(fibers_flux[fibb[0].astype(int)][1])==0:
        print "Fiber %d has no stars inside or nearby..."%fibb[0]
    else:
        contamination[k] = [np.max(fibers_flux[fibb[0].astype(int)][1]), -fibb[1],fibb[2], fiberIMAGEstd[fibb[0].astype(int)]]
    if np.any(np.array(fibers_flux[fibb[0].astype(int)][1] > 0)): #Only printing fibers that have a percentaje contrib by one star bigger than 70%
        print fibers_flux[fibb[0].astype(int)]
        j+=1
print "Fibers above 10%:"
print j #This is just a counter of fiber that met the criteria
contamination = contamination[np.argsort(contamination,axis=0)[:,0]]
#Now checking dispersion as a function of contamination...
div = 4
step = len(contamination)/div
table = []
print np.std(contamination[contamination[:,0]>0][:,2])
for i in xrange(div):
    if len(contamination)-(i+1)*step < step:
        table.append([len(contamination)-i*step,
                      np.std(contamination[:,1][i*step:len(contamination)]),
                      np.std(contamination[:,2][i*step:len(contamination)]),
                      np.min(contamination[:,0][i*step:len(contamination)]),
                      np.max(contamination[:,0][i*step:len(contamination)])])
    else:
        table.append([step,
                  np.std(contamination[:,1][i*step:(i+1)*step]),
                  np.std(contamination[:,2][i*step:(i+1)*step]),
                  np.min(contamination[:,0][i*step:(i+1)*step]),
                  np.max(contamination[:,0][i*step:(i+1)*step])])
np.savetxt("contamination.txt", np.array(table), header = "#kvel mvel min% max %")
np.savetxt("contamin.txt", np.array(contamination), header = "#kvel mvel min% max %")
'''
