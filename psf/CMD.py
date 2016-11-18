import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.table import Table
import os
from scipy.stats import gaussian_kde
VERBOSE=True
if os.path.isfile("CMD.fits"):
    t = Table.read("CMD.fits")
else:
    shallowbluecat2 = "hst_12304_02_wfc3_uvis_f555w/hst_12304_02_wfc3_uvis_f555w_daophot_trm.cat"
    shallowbluecat1 = "hst_12304_01_wfc3_uvis_f555w/hst_12304_01_wfc3_uvis_f555w_daophot_trm.cat"
    shallowredcat2 = "hst_12304_02_wfc3_uvis_f814w/hst_12304_02_wfc3_uvis_f814w_daophot_trm.cat"
    shallowredcat1 = "hst_12304_01_wfc3_uvis_f814w/hst_12304_01_wfc3_uvis_f814w_daophot_trm.cat"


    deepredcat = "hst_10520_04_acs_wfc_f814w/hst_10520_04_acs_wfc_f814w_daophot_trm.cat"
    deepbluecat = "hst_10520_02_acs_wfc_f435w_WRONGAST/hst_10520_02_acs_wfc_f435w_daophot_trm.cat"

    deepred = np.loadtxt(deepredcat)
    deepblue = np.loadtxt(deepbluecat)
    shallowred = np.append(np.loadtxt(shallowredcat1), np.loadtxt(shallowredcat2),axis=0)
    shallowblue = np.append(np.loadtxt(shallowbluecat1), np.loadtxt(shallowbluecat1),axis=0)
    deepred = shallowred
    deepblue = shallowblue
    if VERBOSE:
        print "Using catalogs:"
        print deepredcat
        print deepbluecat
        print "Choosing the stars that carry no FLAG."

    deepred = deepred[deepred[:,-1]==0]
    deepblue = deepblue[deepblue[:,-1]==0]

    if VERBOSE:
        print "Getting closest star from blue catalog to red catalog and storing a Table."

    data=[]
    for star in deepred:
        ra_red, dec_red  = star[2], star[3]
        mag_red = star[5]
        
        ra_blue, dec_blue = deepblue[:,2], deepblue[:,3]
        dist = np.sqrt((ra_red-ra_blue)**2+(dec_red-dec_blue)**2)
        min_dist_arg = np.argmin(dist)
        min_dist = dist[min_dist_arg]
        
        ra_blue, dec_blue = deepblue[min_dist_arg,2], deepblue[min_dist_arg,3]
        mag_blue = deepblue[min_dist_arg, 5]
        data.append([min_dist_arg,
                     min_dist,
                     mag_red,
                     mag_blue,
                     ra_red,
                     dec_red,
                     ra_blue,
                     dec_blue])

    t = Table(np.array(data),names=('ARG_BLUE','MIN_DIST','M_RED','M_BLUE','RA_RED','DEC_RED','RA_BLUE','DEC_BLUE'))

    closest=[]
    if VERBOSE:
        print "Clearing for repetitions on the mapping of min distance RED -> BLUE. Choosing the smallest one."
#Red has more stars so I use it to go over the range of blue.
    for ARG_RED in np.arange(len(t)):
        #Flagging the ones that are not repeated matches for a given star in ARG_RED with 100
        #(so they won't be considered in np.argmin).
        dist_i = t['MIN_DIST'] * (ARG_RED==t['ARG_BLUE']) + 100 * (ARG_RED!=t['ARG_BLUE'])
        closest.append([ARG_RED,np.argmin(dist_i)])
    closest = np.array(closest)

    maskdist = (t['MIN_DIST'] < 0.00005)*(closest[:,1]!=0)
    t = t[maskdist]
    if VERBOSE:
        print t
    t.write("CMD.fits",format='fits')

leo1 = np.loadtxt("coord_shallow_bin1.txt")
leo2 = np.loadtxt("coord_shallow_bin2.txt")
leo3 = np.loadtxt("coord_shallow_bin3.txt")
leo4 = np.loadtxt("coord_shallow_bin4.txt")


t1 = t[np.in1d(t['ARG_BLUE'],leo1)]
t2 = t[np.in1d(t['ARG_BLUE'],leo2)]
t3 = t[np.in1d(t['ARG_BLUE'],leo3)]
t4 = t[np.in1d(t['ARG_BLUE'],leo4)]
#t=t[t['MIN_DIST'] < 0.00001]
plt.scatter(t['M_RED'],t['M_BLUE'],s=0.1)
plt.xlabel("$M_{814}$")
plt.ylabel("$M_{435}$")
plt.show()
x = t['M_BLUE']-t['M_RED']
y = t['M_RED']
xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
plt.gca().invert_yaxis()
#plt.scatter(x, y, c=z, s=0.1, edgecolor='')
plt.scatter(x, y, color = 'k', s=0.1, edgecolor='')
plt.xlabel("$M_{435}-M_{814}$")
plt.ylabel("$M_{814}$")
plt.show()

plt.scatter(t1['M_BLUE']-t1['M_RED'],t1['M_RED'],label='Bin1',alpha=.4)
plt.scatter(t2['M_BLUE']-t2['M_RED'],t2['M_RED'],label='Bin2',alpha=.4)
plt.scatter(t3['M_BLUE']-t3['M_RED'],t3['M_RED'],label='Bin3',alpha=.4)
plt.scatter(t4['M_BLUE']-t4['M_RED'],t4['M_RED'],label='Bin4',alpha=.4)
plt.xlabel("$M_{435}-M_{814}$")
plt.ylabel("$M_{814}$")
plt.gca().invert_yaxis()
plt.legend()
plt.show()
'''
	shallowblue2=np.loadtxt("hst_12304_02_wfc3_uvis_f555w/hst_12304_02_wfc3_uvis_f555w_daophot_trm.cat")
	shallowblue1=np.loadtxt("hst_12304_01_wfc3_uvis_f555w/hst_12304_01_wfc3_uvis_f555w_daophot_trm.cat")
	sdssk = np.loadtxt("../../LeoII/psf/KARL/LeoI_photometry.als")

pltdist.append(distmask)
    pltmagdiff.append(magdiff)
    pltmag.append(deepblue[mindistarg,5])
pltmagdiff,pltmag,pltdist = np.array(pltmagdiff),np.array(pltmag),np.array(pltdist)
#star[distmask]
#plt.hist(pltdist, bins = 1000)
#plt.hist(pltmag, bins=1000)
plt.scatter(pltmagdiff[pltdist],pltmag[pltdist],s=1)
plt.show()
exit()
hstsh = np.append(shallowblue2[:,5], shallowblue1[:,5])
weights = np.ones_like(hstsh)*6/float(len(hstsh))
hstdeep = deepblue[:,5]*float(sys.argv[1])
sdss = sdssk[:,3]

plt.scatter(deepred[:,5],deepred[:,5]-deepblue[:,5])
plt.show()
exit()
coord=np.loadtxt("coord.txt")
plt.hist(coord,bins=70,histtype='step', label="FIBERS",weights=np.ones(len(coord))*10)#,normed=True)
plt.hist(hstdeep, bins=100, histtype='step', label="Deep HST")#,normed=True)
plt.hist(hstsh,bins=170, histtype='step', label="Shallow HST")#, weights=weights)
#plt.hist(sdss, bins=100, histtype='step', label="SDSS",normed=True)
plt.legend()
plt.xlabel("Magnitude")
plt.ylabel("Counts")
plt.show()
'''
