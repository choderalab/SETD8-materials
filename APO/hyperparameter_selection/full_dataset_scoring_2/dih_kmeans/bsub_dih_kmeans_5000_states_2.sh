#!/bin/bash
#BSUB -J dih_kmeans
#BSUB -n 16
#BSUB -R span[ptile=16]
#BSUB -R rusage[mem=3]
#BSUB -W 48:00
#BSUB -o /home/rafal.wiewiora/job_outputs/%J.stdout
#BSUB -eo /home/rafal.wiewiora/job_outputs/%J.stderr

source /home/rafal.wiewiora/.bashrc
cd $LS_SUBCWD

export PYEMMA_NJOBS=16 OMP_NUM_THREADS=1

# Launch job.
python dih_kmeans_newscheme.py 5000 2
