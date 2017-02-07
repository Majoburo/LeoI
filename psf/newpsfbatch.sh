#!/bin/bash
#
#SBATCH --time=06:00:00
#SBATCH -A Hobby-Eberly-Telesco
#SBATCH -p gpu
#SBATCH -N 7
#SBATCH -n 140
#SBATCH -o hello.out
export PYTHONPATH=/home/03237/majoburo/lib/python
PYCODE=$WORK/vw/LeoI/psf/psf.py
ibrun /home/03237/majoburo/anaconda2/bin/python $PYCODE 10000
