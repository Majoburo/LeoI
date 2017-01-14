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
    #for a in glob.glob(sys.argv[1]+"std*_1000.txt"):
    for a in glob.glob("std*_1000.txt"):
        data = np.loadtxt(a)
        stdest.append(data)

    stdest = np.transpose(np.array(stdest), (1, 0, 2))
    for i, fiber in enumerate(stdest):
        stdest[i] = fiber[np.argsort(fiber, axis = 0)[:, 0]]
    #t = Table.read(sys.argv[1]+'Leo_table.fits')
    t = Table.read('Leo_table.fits')
    est_std = np.zeros((len(t),140,2))
    est_std[t['N_Stars'] > 0] = stdest
    table = t
    table.replace_column('Estimated_Std', est_std)
    #table.add_row(table[0])
    #t2 = Table.read('std_1000_cat_1/Leo_table.fits')
    #print t2[60]
    print table
    Table.write(table, 'Leo_table.fits', overwrite = True)
    return

add_estimated_std()
exit()
def func_flux(x):
    return 27.9585 - 0.456327*x + 0.0101273*x**2



t = Table.read('std_10000_cat123/Leo_table.fits')
t = Table.read('Leo_table.fits')
t1 = Table.read('std_1000_cat_1/Leo_table.fits')
t23 = Table.read('std_1000_cat_23/Leo_table.fits')
catmask = (t1['N_Stars'] > 0)*(t23['N_Stars'] > 0)
'''
for i in range(140):
    y2 = t23['Estimated_Std'][catmask][:,i,1]
    y1 = t1['Estimated_Std'][catmask][:,i,1]
    p = np.poly1d(np.polyfit(y2, y1, 2))
    t['Estimated_Std'][:,i,1][t['Catalog'] > 1] = p(t['Estimated_Std'][:,i,1][t['Catalog'] > 1])
    #plt.plot(y2,y1,'ro')
    #plt.plot(y2,p(y2),'ro')
    #plt.show()
table = t[t['N_Stars'] > 0]
Table.write(table, 'std_10000_cat123/Leo_tabule.fits', overwrite = True)
'''
'''
for i, x in enumerate(t23["Estimated_Std"][catmask][:,:,1].T):
    y = np.array(x/np.max(t23["Flux"][catmask],axis=1)*func_flux(np.max(t23["Flux"],axis=1)[catmask])).T
    
    p = np.poly1d(np.polyfit(y, y2[i], 2))
    print i,np.sum(abs(y2[i]-p(y)))
    #range = np.linspace(5, 10, 100)
    #plt.scatter(y,y2[i])
    #plt.plot(range,p(range))
    #plt.show()
exit()
a = np.array([x/np.max(t23["Flux"][catmask],axis=1)*func_flux(np.max(t23["Flux"],axis=1)[catmask]) for x in t23["Estimated_Std"][catmask][:,:,1].T]).T
b = t23['Estimated_Std'][catmask][:,:,1]
print (b-a)[0]

exit()
'''
b_bins = 4
t = Table.read('Leo_table.fits')
loglikelihood = np.zeros((b_bins, t['Estimated_Std'].shape[1],100)) #first row is the xvalues
sigmarange = t['Estimated_Std'][1][:, 0]
vel1 = t['Velocity'][t['Diameter'] == 1.6]
vel2 = t['Velocity'][t['Diameter'] == 3.2]
original = np.zeros((t['Estimated_Std'].shape[1],b_bins,100))
for index in xrange(t['Estimated_Std'].shape[1]):
    std = t['Estimated_Std'][t['Diameter'] == 1.6][:, index, 1]
    mean1 = np.sum((1./std)**2 * vel1)/np.sum((1./std)**2)
    t['Velocity'][t['Diameter'] == 1.6] = abs(vel1 - mean1)
    std = t['Estimated_Std'][t['Diameter'] == 3.2][:, index, 1]
    mean2 = np.sum((1./std)**2 * vel2)/np.sum((1./std)**2)
    t['Velocity'][t['Diameter'] == 3.2] = abs(vel2 - mean1)
    binnedstdplt = []
    for j, minflux in enumerate(np.arange(0,100.,1.)):
        for i in xrange(b_bins):
            binnedstd = t['Estimated_Std'][(t['Bin'] == i+1)*(np.max(t['Flux'],axis=1) > minflux)][:, index, 1]
            binnedvel = t['Velocity'][(t['Bin'] == i+1)*(np.max(t['Flux'],axis=1) > minflux)]
            #binnedstd = t['Estimated_Std'][(t['Bin'] == i+1)*(np.max(t['Flux'],axis=1) > minflux)*(t["Diameter"]==1.6)][:, index, 1]
            #binnedvel = t['Velocity'][(t['Bin'] == i+1)*(np.max(t['Flux'],axis=1) > minflux)*(t["Diameter"]==1.6)]
            loglikelihood[i, index, j] = -np.sum(((binnedvel)/(binnedstd))**2/2) - np.log(np.prod(np.sqrt(2*np.pi)*binnedstd))
'''            original[index, i, j] = np.sqrt(np.sum(binnedvel*binnedvel)/len(binnedvel-1))


original = np.mean(original, axis=0)
plt.plot(original[0], label="Bin 1")
plt.plot(original[1], label="Bin 2")
plt.plot(original[2], label="Bin 3")
plt.plot(original[3], label="Bin 4")
plt.xlabel("Mininum Max Flux")
plt.ylabel("Mean Velocity Dispersion")
plt.title("PATCHED CATALOGS")
plt.legend()
plt.show()
'''
loglikelihood = np.mean(loglikelihood, axis=1)
plt.plot(loglikelihood[0], label="Bin 1")
plt.plot(loglikelihood[1], label="Bin 2")
plt.plot(loglikelihood[2], label="Bin 3")
plt.plot(loglikelihood[3], label="Bin 4")
plt.xlabel("Mininum Max Flux")
plt.ylabel("Mean Velocity Dispersion")
plt.title("PATCHED CATALOGS")
plt.legend()
plt.show()
#print loglikelihood
#exit()
'''
'''
lowersigma = np.min(sigmarange)
uppersigma = np.max(sigmarange)

fitrange = np.arange(lowersigma, uppersigma - 1, 1/100.)
maxfluxdependance = np.zeros((4,100))
for minflux in range(100):
    likelihoodfit = [fitrange]
    for j, bin in enumerate(loglikelihood[:, :, minflux]):
        #func = np.polyfit(sigmarange, bin, 20)
        #f = np.poly1d(func)
        f = interp1d(sigmarange, bin, kind = 'linear')
        window_size, poly_order = 1251, 4
        yy_sg = savgol_filter(f(fitrange), window_size, poly_order)
        plt.plot(sigmarange,bin)
        plt.plot(fitrange, yy_sg, 'k', label= "Smoothed curve")
        #plt.plot(fitrange, f(fitrange))
        plt.show()
        #test = 100. * np.exp(f(fitrange))/np.sum(np.exp(f(fitrange)))
        test = 100. * np.exp(yy_sg)/np.sum(np.exp(yy_sg))
        likelihoodfit.append(test)
        #plt.plot(sigmarange, np.exp(bin)/np.sum(np.exp(bin)))
        #plt.plot(fitrange, test, '-', label = 'BIN %d'%j)
        #plt.legend()
        #plt.show()
    likelihoodfit =  np.array(likelihoodfit)

    x2 = np.arange(len(likelihoodfit.T[1:-1]))
    h = 1
    for j, bin in enumerate(likelihoodfit[1:b_bins + 1]):
        hessian2 = 10000*(bin[x2+h] - 2 * bin[x2] + bin[x2-h])/ h**2
        #plt.plot(likelihoodfit[0][10:-1], bin[10:-1])
        #plt.plot(likelihoodfit[0][10:-1], hessian2[9:])
        #plt.plot(likelihoodfit[0][10:-1], bin[10:-1]*hessian2[9:])
        #print likelihoodfit[0][np.argmax(bin)]
        #print np.sqrt(1/(- np.sum(bin[1:-1] * hessian2, axis = 0)))
        #print likelihoodfit[0][np.argmax(bin)]
        maxfluxdependance[j, minflux] = likelihoodfit[0][np.argmax(bin)]
        #plt.show()
plt.plot(maxfluxdependance.T)
plt.plot(maxfluxdependance[0], label="Bin 1")
plt.plot(maxfluxdependance[1], label="Bin 2")
plt.plot(maxfluxdependance[2], label="Bin 3")
plt.plot(maxfluxdependance[3], label="Bin 4")
plt.xlabel("Minimum Max Flux")
plt.ylabel("Estimated Velocity Dispersion")
plt.title("PATCHED CATALOGS")
plt.legend()
plt.show()
