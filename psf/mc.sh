#!/bin/sh
iterations=10000
sigmamin=5
sigmamax=40
sigmastep=0.25
#rm Leo_table.fits
#python psf.py 10 3
#rm std_10.00_3.txt
python writebatch.py $iterations $sigmamin $sigmamax $sigmastep
sbatch psfbatch.slurm
