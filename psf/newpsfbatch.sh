#!/bin/bash
#
#SBATCH --time=00:05:00
#SBATCH -A Hobby-Eberly-Telesco
#SBATCH -p vis
#SBATCH -N 7
#SBATCH -n 140
#SBATCH -o hello.out
export PYTHONPATH=/home/03237/majoburo/lib/python
PYCODE=$WORK/vw/LeoI/psf/psf.py
ibrun /home/03237/majoburo/anaconda2/bin/python $PYCODE 10 3 3
