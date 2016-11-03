import numpy as np
from astropy.table import Table
t=Table.read("Leo_table.fits")
np.savetxt("bin1.txt",np.array(zip(t["RA"][t["Bin"]==1],t["DEC"][t["Bin"]==1])))
np.savetxt("bin2.txt",np.array(zip(t["RA"][t["Bin"]==2],t["DEC"][t["Bin"]==2])))
np.savetxt("bin3.txt",np.array(zip(t["RA"][t["Bin"]==3],t["DEC"][t["Bin"]==3])))
np.savetxt("bin4.txt",np.array(zip(t["RA"][t["Bin"]==4],t["DEC"][t["Bin"]==4])))
