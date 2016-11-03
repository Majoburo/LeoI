import numpy as np
import sys

iterations = int(sys.argv[1])
sigmamin = float(sys.argv[2])
sigmamax = float(sys.argv[3])
sigmastep = float(sys.argv[4])
f = open("psfbatch",'w')
for i in np.arange(sigmamin, sigmamax, sigmastep):
    f.write('python psf.py %.2f %d\n'%(i,iterations))
f.close()
