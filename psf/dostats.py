import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import scipy.integrate
from numpy.polynomial.hermite import hermfit, hermval
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import re
import glob
from astropy.table import Table
#Putting all MC std into a Table


def add_estimated_std():
    stdest = []
    for a in glob.glob("std*txt"):
        data = np.loadtxt(a)
        stdest.append(data)

    stdest = np.transpose(np.array(stdest), (1, 0, 2))
    for i, fiber in enumerate(stdest):
        stdest[i] = fiber[np.argsort(fiber, axis = 0)[:, 0]]

    t = Table.read('Leo_table.fits')
    table = t[t['N_Stars'] > 0]
    table.replace_column('Estimated_Std', stdest)
    Table.write(table, 'Leo_tabule.fits', overwrite = True)
    return

add_estimated_std()
b_bins = 4
t = Table.read('Leo_tabule.fits')
loglikelihood = np.zeros((b_bins, t['Estimated_Std'].shape[1])) #first row is the xvalues
sigmarange = t['Estimated_Std'][1][:, 0]
vel1 = t['Velocity'][t['Diameter'] == 1.6]
vel2 = t['Velocity'][t['Diameter'] == 3.2]
for index in xrange(t['Estimated_Std'].shape[1]):
    std = t['Estimated_Std'][t['Diameter'] == 1.6][:, index, 1]
    mean1 = np.sum((1./std)**2 * vel1)/np.sum((1./std)**2)
    t['Velocity'][t['Diameter'] == 1.6] = abs(vel1 - mean1)
    std = t['Estimated_Std'][t['Diameter'] == 3.2][:, index, 1]
    mean2 = np.sum((1./std)**2 * vel2)/np.sum((1./std)**2)
    t['Velocity'][t['Diameter'] == 3.2] = abs(vel2 - mean1)
    binnedstdplt = []
    for i in xrange(b_bins):
        binnedstd = t['Estimated_Std'][(t['Bin'] == i+1)][:, index, 1]
        binnedvel = t['Velocity'][(t['Bin'] == i+1)]
        loglikelihood[i, index] = -np.sum(((binnedvel)/(binnedstd))**2/2) - np.log(np.prod(np.sqrt(2*np.pi)*binnedstd))

plt.plot(sigmarange, loglikelihood[3])
plt.show()
lowersigma = np.min(sigmarange)
uppersigma = np.max(sigmarange)
'''
range = np.array([float(re.split("std|_",a)[2]) for a in glob.glob("std*txt")])
range = np.sort(range)
print range
lowersigma = 0.5
uppersigma = 20.
#range = np.arange(lowersigma, uppersigma, 0.5)

b_bins = 4
loglikelihood = np.zeros((b_bins+1, len(range))) #first row is the xvalues
loglikelihood[0] = range
binnedstdplt=[]
'''
'''
for index, sigma in enumerate(range):
    file1 = np.loadtxt("m_std%.1f_%s.txt"%(sigma, sys.argv[1]))
    #data1 = file1[file1[:,2] > 0]
    std.append(file1[:,2]/10.)
plt.plot(range,std)
plt.show()
'''
'''
#exit()
table = Table.read("Leo_table.fits")
t = table[table['N_Stars']>0]
vel = t['Velocity']
bins = t['Bin']
for index, sigma in enumerate(range):
    realstd = float(sigma) # I fine-grained to a tenth of a km/s

    data1 = np.loadtxt("std_%.2f_%s.txt"%(sigma, sys.argv[1]))
    #data1 = file1[file1[:,2] > 0]
    #data1 = file1
    std = data1[:,1]
    #vel = data1[:,1]
    #std[std == 0.] = np.mean(std[std!=0.]) #Correcting by average std for the ones i have no data on.
    mean1 = np.sum((1./std)**2*vel)/np.sum((1./std)**2)
    data1[:,1] = abs(vel - mean1)
    #databinned = databinned[np.argsort(databinned, axis=0)[:,-1]]
    binnedstdplt = []
    for i in xrange(b_bins):
        #binmask = (databinned[:,-1]>i*b_bins)*(databinned[:,-1]<=(i+1)*b_bins)
        binnedstd = data1[bins == i][:,1]
        #databinned = databinned[databinned[:,2] > 0]
        binnedvel = vel[bins == i]
        #databinned = databinned[databinned[:,2] > 0]
        if sigma==10:
            print "Bin %d"%i
            print np.sqrt(np.sum(binnedvel**2)/(len(binnedvel)))
        #std[std == 0.] = np.mean(std[std!=0.]) #Correcting by average std for the ones i have no data on.
        loglikelihood[i+1, index] = -np.sum(((binnedvel)/(binnedstd))**2/2)-np.log(np.prod(np.sqrt(2*np.pi)*binnedstd))
        if i == 2:
            binnedstdplt.append([sigma, binnedstd])
print binnedstdplt
'''
fitrange = np.arange(lowersigma, uppersigma - 1, 1/100.)
print len(fitrange)/20
likelihoodfit = [fitrange]
for j, bin in enumerate(loglikelihood):
    #func = np.polyfit(sigmarange, bin, 20)
    f = interp1d(sigmarange, bin, kind = 'linear')
    window_size, poly_order = 1251, 4
    yy_sg = savgol_filter(f(fitrange), window_size, poly_order)
    #f = np.poly1d(func)
    plt.plot(sigmarange,bin)
    plt.plot(fitrange, yy_sg, 'k', label= "Smoothed curve")
    plt.plot(fitrange, f(fitrange))
    plt.show()
    test = 100. * np.exp(f(fitrange))/np.sum(np.exp(f(fitrange)))
    test = 100. * np.exp(yy_sg)/np.sum(np.exp(yy_sg))
    likelihoodfit.append(test)
    plt.plot(sigmarange, np.exp(bin)/np.sum(np.exp(bin)))
    plt.plot(fitrange, test, '-', label = 'BIN %d'%j)
    plt.show()
likelihoodfit =  np.array(likelihoodfit)
np.savetxt("likelihoodfit.txt", likelihoodfit.T)

x2 = np.arange(len(likelihoodfit.T[1:-1]))
h = 1
for j, bin in enumerate(likelihoodfit[1:b_bins + 1]):
    hessian2 = 10000*(bin[x2+h] - 2 * bin[x2] + bin[x2-h])/ h**2
    plt.plot(likelihoodfit[0][10:-1], bin[10:-1])
    plt.plot(likelihoodfit[0][10:-1], hessian2[9:])
    plt.plot(likelihoodfit[0][10:-1], bin[10:-1]*hessian2[9:])
    print likelihoodfit[0][np.argmax(bin)]
    print np.sqrt(1/(- np.sum(bin[1:-1] * hessian2, axis = 0)))

    plt.show()

