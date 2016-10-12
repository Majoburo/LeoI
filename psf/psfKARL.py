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
"""
DISCLAIMER:
    This should be used with txt files that list pixel locations of fibers and stars on the same image.

"""
##Volume under a 2d Gaussian fo Magnitude M and standard deviation sig
##intersected by a circle of radius fR at distance cD (center to center)

def integrand1(r,cD,fR,sig,M):
    return 2 * np.pi * r *  10 ** ((25.-M)/2.5) * iv(0, - cD * r/(sig*sig) ) * np.exp( - (cD * cD + r * r)/(2 * sig * sig))


def integrate(cD,fR,sig,M):
        return quad(integrand1, 0, fR, args = (cD, fR, sig, M))[0]


#Pixel coords of the stars identified by daophot in ds9.fits
leophot = np.loadtxt("KARL/LeoI_photometry.als")
mstars = leophot[:,1:3]
mateostars = np.loadtxt("")
#Pixel coords of the fiber centers in ds9.fits
fiberpos1 = np.loadtxt("KARL/fibersposss.txt")
fiberpos = fiberpos1[:,3:5]
fibR = 35*2
#mask around the IFU to take only those stars into the selection
mask = (mstars[:,1] < np.max(fiberpos[:,1]) + fibR)*(mstars[:,1] > np.min(fiberpos[:,1]) - fibR)*(mstars[:,0] > np.min(fiberpos[:,0]) - fibR)*(mstars[:,0] < np.max(fiberpos[:,0]) + fibR)
hdu = fits.open("KARL/ds9.fits")

sig = toolbox.degtopix(hdu,1.5/2.35/3600.)
fR = toolbox.degtopix(hdu, (1.6 + 0.5)/3600.)
fibers_flux = []
cont_number = []
maxflux = []
totalflux =[]
for i,fiber in enumerate(fiberpos):
    mask_nb_fib = (fiber[0] - mstars[mask][:,0])**2+(fiber[1] - mstars[mask][:,1])**2 <= fR**2
    cont_number.append(sum(mask_nb_fib*1))
    vec_dist = (fiber - mstars[mask][mask_nb_fib])
    cD = np.sqrt(vec_dist[:,0]**2+vec_dist[:,1]**2)
    M = leophot[mask][mask_nb_fib][:,3]
    flux=[]
    for j, dist in enumerate(cD):
        #print i,cD[j],integrate(cD=cD[j],fR = fR, sig = sig, M=M[j])/(2 * 10 ** ((25.-M[j])/2.5) * sig * sig * np.pi) * 100
        flux.append(integrate(cD[j],fR = fR, sig = sig, M=M[j]))
    fibers_flux.append([i,flux/np.sum(flux)*100.])
    totalflux.append(np.sum(flux))
    if len(flux)!=0:
        maxflux.append(np.max(flux))
    else:
        maxflux.append(0)
maxflux = np.array(maxflux)
totalflux = np.array(totalflux) + 5.30523
maxflux = maxflux/totalflux*100.
print np.mean(cont_number), np.std(cont_number)
np.savetxt("mateostarsifu.txt",mstars[mask])
np.savetxt("maxfluxKARL.txt",np.array(maxflux))
corrfib = np.loadtxt("../compare/fibers.txt")
print corrfib
j=0
contamination = np.zeros((len(corrfib),3))
for k,fibb in enumerate(corrfib):
    if len(fibers_flux[fibb[0].astype(int)][1])==0:
        print "Fiber %d has no stars inside or nearby..."%fibb[0]
    else:
        contamination[k] = [np.max(fibers_flux[fibb[0].astype(int)][1]), -fibb[1],fibb[2]]
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
