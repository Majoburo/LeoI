import psf
import configobj
from astropy.table import Table
config = configobj.ConfigObj("config.ini")
import re
import numpy as np
from astropy.io import fits
from astropy import wcs
import toolbox
from scipy.integrate import quad
from scipy.special import iv

maxstars = 50
reg_order = np.loadtxt("reg_order.txt")

def integrand1(r,cD,fR,sig,M):
    return 2 * np.pi * r *  10 ** ((25.-M)/2.5)*iv(0, - cD * r/(sig*sig) )*np.exp( - (cD * cD + r * r)/(2 * sig * sig))

def integrate(cD,fR,sig,M):
    return quad(integrand1, 0, fR, args = (cD, fR, sig, M))[0]

def VW_IFU_pos(regionsfile):
    regions = []
    with open(regionsfile) as f:
        for line in f.readlines()[4:]:
            line = re.split("[,\(\)]",line)
            regions.append([float(line[1]),float(line[2])])
    regions = np.array(regions)
    return regions

def getflux(fitsfile, catalogfile, POS):
    catalog = np.loadtxt(catalogfile)
    catalog = catalog[catalog[:,-1]==0] #only taking detections without a flag warning
    catalog = catalog[catalog[:,5] < 26]
    #fiber = table[table['Bin'] > 0] #== float(sys.argv[3])]
    RA = POS[:,1]
    DEC = POS[:,2]
    #The median seeing for fibers in Hectochelle: 0.7.
    #The median seeing for fibers in VW: 1.5

    with fits.open(fitsfile) as hdu:
        w = wcs.WCS(hdu[1].header)
        catpixcrd = w.wcs_world2pix(np.array(catalog[:,2:4], dtype='float64'),1)
        #catpixcrd[:,0] = catpixcrd[:,0] - 2.1118999999998778
        #catpixcrd[:,1] = catpixcrd[:,1] - 48.02390000000014
        #catalog[:,2:4] = w.wcs_pix2world(catpixcrd, 1)
        #np.savetxt(catalogfile, catalog)
        starpixcrd = w.wcs_world2pix(zip(RA,DEC),1)
        hdus = []
        hdus.append(hdu[1])
        pixsig = toolbox.degtopix(hdus, 1.5/2.35/3600.) #Divided by 2.35 cause seeing is 2.35*sigma
        pixfR = toolbox.degtopix(hdus, (3.2/2.+0.5)/3600.) #Plus some error disk of radius 0.5
    mpixfR = np.max(pixfR)
    pixsig = np.ones(len(POS))*pixsig
    pixfR = np.ones(len(POS))*pixfR
    leophot = np.array(zip(catalog[:,7], catpixcrd[:,0], catpixcrd[:,1]))
    starpos = np.array(zip(POS[:,0], starpixcrd[:,0], starpixcrd[:,1], pixsig, pixfR))
    mask = ((leophot[:,2] < np.max(starpos[:,2]) + mpixfR*2)*(leophot[:,2] > np.min(starpos[:,2]) - mpixfR*2)*
            (leophot[:,1] < np.max(starpos[:,1]) + mpixfR*2)*(leophot[:,1] > np.min(starpos[:,1]) - mpixfR*2))

    stars_flux = np.zeros((len(starpos), maxstars))
    n_stars = np.zeros(len(starpos))
    #f=open("coord_deep_bin%s.txt"%sys.argv[3],'ab')
    for i, star in enumerate(starpos):
        mask_nb_fib = (star[1] - leophot[mask][:,1])**2 + (star[2] - leophot[mask][:,2])**2 < star[4]**2
        #f.write("%f /n"%leophot[:,0][mask][mask_nb_fib])
        #np.savetxt(f, leophot[:,0][mask][mask_nb_fib])
        #np.savetxt(f, np.arange(len(mask_nb_fib))[mask_nb_fib])
        vec_dist = (star[1:3] - leophot[:,1:3][mask][mask_nb_fib])
        cD = np.sqrt(vec_dist[:,0]**2 + vec_dist[:,1]**2)
        M = leophot[mask][mask_nb_fib][:,0]
        flux = np.zeros(maxstars)
        for j, dist in enumerate(cD):
            if j > maxstars:
                print "You have more than %d stars in a fiber! All hope is lost."%maxstars
                exit()
            flux[j] = integrate(dist, fR =  star[4], sig =  star[3], M=M[j])
            n_stars[i] = j+1
        if n_stars[i] > 0:
            stars_flux[i] = flux/np.sum(flux)*100.
    #f.close()
    #exit()
    return n_stars, stars_flux



POS = VW_IFU_pos("../regions/LEOI_field.reg")
#POS = POS[reg_order.astype(int)]

#t = Table.read("Leo_table.hdf5")
#print t[t['Diameter']>3]
POS = np.array(zip(reg_order, POS[:,0], POS[:,1]))
#print psf.VW_IFU_pos("../../LeoII/regions/LeoI_vwDec13b_field_true.reg")
#print np.read("../regions/LEOI_field.reg")
for field in config['fields']:
    n_stars, stars_flux = getflux(config['fields'][field]['fits'], config['fields'][field]['catalog'], POS)
    print n_stars, stars_flux


