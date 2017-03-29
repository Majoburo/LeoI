#!/usr/bin/env python
from astropy.table import Table
from astropy.coordinates import Angle
from astropy.io import fits
from astropy import wcs
from scipy.integrate import quad
from scipy.special import iv
from scipy.signal import gaussian
from scipy import interpolate
import matplotlib.pyplot as plt
import random
import os.path
import numpy as np
import configobj
import re
import toolbox
import srebin
import time
import sys
import h5py
from mpi4py import MPI

config = configobj.ConfigObj("config.ini")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

range_start = float(config["mc"]["min_speed"])
range_stop  = float(config["mc"]["max_speed"])
step        = float(config["mc"]["step_speed"])
steps_per_proc = np.abs((range_stop-range_start)/step)/size

my_start =  rank*steps_per_proc*step + range_start
my_stop = (rank+1)*steps_per_proc*step + range_start


maxstars = int(config["mc"]["max_stars"])

enhance = int(config["mc"]["enhance"])
lenarray = 1000*enhance

template, trange =  toolbox.getdata(config["vw"]["template"])
newrange = np.linspace(trange[0], trange[-1], num = lenarray)
newdata = interpolate.InterpolatedUnivariateSpline(trange, template)
template = newdata(newrange)
trange = newrange

#Flux integration functions
def integrand1(r,cD,fR,sig,M):

    return 2 * np.pi * r *  10 ** ((25.-M)/2.5)*iv(0, - cD * r/(sig*sig) )*np.exp( - (cD * cD + r * r)/(2 * sig * sig))

def integrate(cD,fR,sig,M):

    return quad(integrand1, 0, fR, args = (cD, fR, sig, M))[0]

#Reading functions
def VW_IFU_pos(regionsfile):
    regions = []
    with open(regionsfile) as f:
        for line in f.readlines()[::2][4:]:
            line = re.split("[,\(\)\{\}]",line)
            regions.append([int(line[4]),float(line[1]),float(line[2])])
    regions = np.array(regions)
    regions = regions[np.argsort(regions, axis=0)[:,0]]

    return regions

def mateo_data(mateotable):
    t = Table.read(mateotable, format="csv",
                       names=("Fiber_n","R_b","PAb","vel","svel","RA","DEC","M","sM","SRCc"),
                       comment='#')
    ra = Angle([a+" hours" for a in t['RA']]).degree
    dec = Angle([a+" degree" for a in t['DEC']]).degree
    vel = t['vel']
    mateo_data = np.column_stack((np.arange(len(ra)), ra, dec, vel))

    return mateo_data

def ellipse_major_axis(ab, angle, x, y):
    sin = np.sin(np.deg2rad(angle))
    cos = np.cos(np.deg2rad(angle))
    return np.sqrt(((x*cos+y*sin)/ab)**2 + (x*sin-y*cos)**2)


def binfibers(fitsfile, table):
#I use a fits file just for the convinience of calculating this in pixels.
    vwfibers = table[table['Diameter']==3.2]
    center = [np.array([config['Parameters']['ellipse_center_RA'], config['Parameters']['ellipse_center_DEC']],dtype = float)]
    ab = float(config['Parameters']['ab'])
    angle = float(config['Parameters']['angle'])
    bins = int(config['Parameters']['bins'])
    with fits.open(fitsfile) as hdu:
        w = wcs.WCS(hdu[0].header)
        vwpixcrd = w.wcs_world2pix(zip(vwfibers['RA'], vwfibers['DEC']),1)
        fiberpixcrd = w.wcs_world2pix(zip(table['RA'], table['DEC']), 1)
        centerpix = w.wcs_world2pix(center,1)
    furtherxy = vwpixcrd[np.argmax(np.sum((centerpix - vwpixcrd)**2,axis=1))]
    furtherxy = (furtherxy - centerpix)[0]
    b = np.sqrt((furtherxy[0]/ab)**2 + (furtherxy[1])**2)  #Taking an extra bins outside the ifu
    fiberbin = np.zeros(len(fibers))
    jj = 0
    for i in xrange(bins):
        x, y = (centerpix - fiberpixcrd)[:,1], (centerpix - fiberpixcrd)[:,0]
        ell_b = ellipse_major_axis(ab, angle, x, y)
        maskfibers = 0
        j = 0
        while np.sum(maskfibers) <= 27:
            j = j + 1
            maskfibers = (ell_b <= ((j+jj)/100.)*b/3.)*(ell_b > (jj/100.)*b/3.)

        jj = j+jj
        distance = ell_b[maskfibers]
        fiberbin = fiberbin + maskfibers * (i+1)
        distance = np.sum(distance)/len(distance)
        with fits.open(fitsfile) as hdu:
            print "Average distance in arcsec in bin %d: %f"%(i, toolbox.pixtodeg(hdu, distance)*3600)
    table['Bin'] = fiberbin
    return

def getflux(fitsfile, catalogfile, table):
    catalog = np.loadtxt(catalogfile)
    catalog = catalog[catalog[:,-1]==0] #only taking detections without a flag warning
    catalog = catalog[catalog[:,5] < 26]
    fiber = table[table['Bin'] > 0]
    RA = fiber['RA']
    DEC = fiber['DEC']
    #The median seeing for fibers in Hectochelle: 0.7.
    #The median seeing for fibers in VW: 1.5

    with fits.open(fitsfile) as hdu:
        w = wcs.WCS(hdu[1].header)
        catpixcrd = w.wcs_world2pix(np.array(catalog[:,2:4], dtype='float64'),1)
        starpixcrd = w.wcs_world2pix(zip(RA,DEC),1)
        hdus = []
        hdus.append(hdu[1])
        pixsig = toolbox.degtopix(hdus, fiber['Seeing']/2.35/3600.) #Divided by 2.35 cause seeing is 2.35*sigma
        pixfR = toolbox.degtopix(hdus, (fiber['Diameter']/2.+0.5)/3600.) #Plus some error disk of radius 0.5
    mpixfR = np.max(pixfR)
    leophot = np.array(zip(catalog[:,7], catpixcrd[:,0], catpixcrd[:,1]))
    starpos = np.array(zip(fiber['Fiber'], starpixcrd[:,0], starpixcrd[:,1], pixsig, pixfR))
    mask = ((leophot[:,2] < np.max(starpos[:,2]) + mpixfR*2)*(leophot[:,2] > np.min(starpos[:,2]) - mpixfR*2)*
            (leophot[:,1] < np.max(starpos[:,1]) + mpixfR*2)*(leophot[:,1] > np.min(starpos[:,1]) - mpixfR*2))

    stars_flux = np.zeros((len(starpos), maxstars))
    n_stars = np.zeros(len(starpos))
    for i, star in enumerate(starpos):
        mask_nb_fib = (star[1] - leophot[mask][:,1])**2 + (star[2] - leophot[mask][:,2])**2 < star[4]**2
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
    return n_stars, stars_flux

def calculatestd(t,guess_speed):

        start = time.time()
        fiberIMAGEstd = []
        iterations = int(sys.argv[1])
        guess_speed = float(guess_speed)*enhance
        tmpmean = np.mean(template)
        tmpstd = np.std(template)
        tmp = (template - tmpmean)/tmpstd
        for fiber in t[t['N_Stars'] > 0]:
            flux = fiber['Flux'][:fiber['N_Stars']]
            newtmp = np.zeros((len(flux), len(tmp)))
            maxgauss = []
            for f in xrange(iterations):
                for i,w in enumerate(flux/np.sum(flux)):
                    shift = int(random.gauss(0, guess_speed))
                    newtmp[i] = np.take(w*tmp, np.arange(shift, len(tmp) + shift),mode = 'wrap')
                    #newtmp[i] = srebin.loglogSpl(trange, tmp, zshift = shift)[1]
                    #plt.plot(newtmp[i])
                    #plt.show()
                sumtmp = np.sum(newtmp, axis = 0)
                sumtmp = (sumtmp - np.mean(sumtmp))/np.std(sumtmp)
                corr = np.correlate(sumtmp, tmp, mode = 'same')[len(tmp)/4:len(tmp)*3/4]
                #plt.plot(sumtmp[len(tmp)/4:len(tmp)*3/4])
                #plt.plot(tmp[len(tmp)/4:len(tmp)*3/4])
                #plt.plot(corr)
                #plt.show()
                maxgauss.append(np.argmax(corr))
            fiberIMAGEstd.append(np.std(maxgauss))
        fiberIMAGEstd = np.array(fiberIMAGEstd)
        astd = np.stack((np.ones(len(fiberIMAGEstd))*guess_speed/10., fiberIMAGEstd/10.),axis=-1)
        end = time.time()
        print "This calculation took (in seconds):"
        print(end - start)
        astd.shape = (len(fiberIMAGEstd), 1, 2)
        return astd

if __name__ == '__main__':
    if os.path.isfile("Leo_table.hdf5"):
        curr_rank=0
        t = Table.read("Leo_table.hdf5")
        astd = np.zeros((len(t), 1, 2))
        first_guess = True
        for guess in np.arange(my_start,my_stop,step):
            astd[t['N_Stars'] > 0] = calculatestd(t,guess)
            if not first_guess:
                ostd = np.array(t['Estimated_Std'])
                astd = np.concatenate((ostd, astd), axis=1)
                first_guess = False
            t.replace_column('Estimated_Std', astd)
        all_stds = comm.gather(t['Estimated_Std'],root=0)
        if rank == 0:
            astd = np.concatenate(all_stds,axis=1)
            t.replace_column('Estimated_Std', astd)
            t.write("Leo_table.hdf5", path='updated_data', overwrite=True)
    else:
        if rank==0:
            #Loading vw position an kinematic data:
            vwpos = VW_IFU_pos(config["vw"]["regionsfile"])
            vwkin = np.loadtxt(config['vw']["data"])
            # Selecting only star fibers in vw.
            vwstarfibers = np.in1d(vwpos[:,0], vwkin[:,0])
            vwpos = vwpos[vwstarfibers]
            # Appending fiber diameter (for only only distinction between vw and hectochelle fibers)
            vwfibers = np.column_stack((vwpos, vwkin[:,2], np.ones(len(vwpos))*3.2, np.ones(len(vwpos))*1.5))
            mateodata = mateo_data(config['mateo']["data"])
            mateofibers = np.column_stack((mateodata, np.ones(len(mateodata))*1.6, np.ones(len(mateodata))*0.7))

            fibers = np.concatenate((vwfibers, mateofibers), axis = 0)
            t = Table(np.zeros((len(fibers),11)),
                       names=('Fiber','RA','DEC','Velocity','Diameter','Seeing','Bin','Catalog','N_Stars',
                              'Flux', 'Estimated_Std'), dtype=('int16',float,float,float,float,float,'int16','int16','int16',float,float))
            t.replace_column('Flux', np.zeros((len(t),maxstars)))
            t.replace_column('Estimated_Std', np.zeros((len(t), 1, 2)))
            t['Fiber'] = fibers[:,0]
            t['RA'] = fibers[:,1]
            t['DEC'] = fibers[:,2]
            t['Velocity'] = fibers[:,3]
            t['Diameter'] = fibers[:,4]
            t['Seeing'] = fibers[:,5]
            binfibers(config["vw"]["fits"], t)
            n_starslist = []
            stars_fluxlist = []
            for field in config['fields']:
                n_stars, stars_flux = getflux(config['fields'][field]['fits'], config['fields'][field]['catalog'], t)
                n_starslist.append(n_stars)
                stars_fluxlist.append(stars_flux)
            n_starslist = np.array(n_starslist)
            stars_fluxlist = np.array(stars_fluxlist)

            usedcat = np.argmax(n_starslist,axis=0)

            stars_fluxlist = np.transpose(stars_fluxlist,(1, 0, 2))
            n_starslist = [row[usedcat[i]] for i,row in enumerate(n_starslist.T)]
            stars_fluxlist = [row[usedcat[i]] for i,row in enumerate(stars_fluxlist)]
            n_starslist = np.array(n_starslist)
            t['Catalog'][t['Bin'] > 0] = usedcat+1
            t['N_Stars'][t['Bin'] > 0] = n_starslist
            t['Flux'][t['Bin'] > 0] = stars_fluxlist
            t['Estimated_Std'][t['N_Stars'] > 0] = calculatestd(t, range_start)
            t.write("Leo_table.hdf5", path='updated_data')
            print t

