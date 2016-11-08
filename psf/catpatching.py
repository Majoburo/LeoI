from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
t1=Table.read("std_1000_cat_1/Leo_table.fits")
t23=Table.read("std_1000_cat_23/Leo_table.fits")

mask=(np.max(t1["Flux"],axis=1)*np.max(t23["Flux"],axis=1)>0)

def func_flux(x):
    return 27.9585 - 0.456327*x + 0.0101273*x**2

gs=gridspec.GridSpec(2, 1, height_ratios=[100, 70])
a0 = plt.subplot(gs[0])
a1 = plt.subplot(gs[1])
a0.plot(np.max(t23["Flux"],axis=1)[mask],np.max(t1["Flux"],axis=1)[mask],'ro')
a0.set_ylabel("Max Flux Deeper Cat")
a1.set_xlabel("Max Flux Shallow Cat")
a1.set_ylabel("Fit - Max Flux Deeper Cat")
a0.plot(np.arange(10,100,1), func_flux(np.arange(10,100,1)))
a1.plot(np.max(t23["Flux"],axis=1)[mask], -func_flux(np.max(t23["Flux"],axis=1)[mask])+np.max(t1["Flux"],axis=1)[mask],'ro')
a1.plot(np.arange(10,100),np.zeros(90))
a1.set_ylim([-50,20])
a0.set_ylim([0,100])
plt.show()
plt.scatter(np.mean(t23["Estimated_Std"][mask],axis=1)[:,1],np.mean(t1["Estimated_Std"][mask],axis=1)[:,1])

plt.ylabel(r'$\bar{{\sigma}}_{DEEP}$')
plt.xlabel(r'$\bar{{\sigma}}_{SHALLOW}$')
plt.show()
plt.plot(np.mean(t23["Estimated_Std"][mask],axis=1)[:,1]/np.max(t23["Flux"][mask],axis=1)*np.max(t1["Flux"][mask],axis=1),np.mean(t1["Estimated_Std"][mask],axis=1)[:,1],'ro',label=r'$MAXFLUX_{DEEP}$')
plt.plot(np.mean(t23["Estimated_Std"][mask],axis=1)[:,1]/np.max(t23["Flux"][mask],axis=1)*func_flux(np.max(t23["Flux"],axis=1)[mask]),np.mean(t1["Estimated_Std"][mask],axis=1)[:,1],'bo',label=r'$FIT(MAXFLUX_{DEEP})$')
plt.legend()
plt.xlabel(r'$\frac{\bar{{\sigma}}_{SHALLOW}*MAXFLUX_{DEEP}}{MAXFLUX_{SHALLOW}}$')
plt.ylabel(r'$\bar{{\sigma}}_{DEEP}$')
plt.show()
