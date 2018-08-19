#!/bin/bash
#BSUB -J dist_cross_kmeans
#BSUB -n 16
#BSUB -R span[ptile=16]
#BSUB -R rusage[mem=4]
#BSUB -W 48:00
#BSUB -o /home/rafal.wiewiora/job_outputs/%J.stdout
#BSUB -eo /home/rafal.wiewiora/job_outputs/%J.stderr

source /home/rafal.wiewiora/.bashrc
cd $LS_SUBCWD

export PYEMMA_NJOBS=16 OMP_NUM_THREADS=1

# Launch job.
python dist_cross_kmeans_newscheme.py kinetic 1 50 2
