#!/bin/bash
#BSUB -J dist_kmeans
#BSUB -n 32
#BSUB -R span[ptile=32]
#BSUB -R rusage[mem=7]
#BSUB -W 48:00
#BSUB -o /home/rafal.wiewiora/job_outputs/%J.stdout
#BSUB -eo /home/rafal.wiewiora/job_outputs/%J.stderr

source /home/rafal.wiewiora/.bashrc
cd $LS_SUBCWD

export PYEMMA_NJOBS=32 OMP_NUM_THREADS=1

# Launch job.
python dist_kmeans_newscheme.py 100 0
