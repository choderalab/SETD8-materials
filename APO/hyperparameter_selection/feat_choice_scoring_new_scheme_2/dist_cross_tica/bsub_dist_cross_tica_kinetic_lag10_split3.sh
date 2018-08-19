#!/bin/bash
#BSUB -J dist_cross_tica
#BSUB -n 1
#BSUB -R span[ptile=1]
#BSUB -R rusage[mem=50]
#BSUB -W 12:00
#BSUB -o /home/rafal.wiewiora/job_outputs/%J.stdout
#BSUB -eo /home/rafal.wiewiora/job_outputs/%J.stderr

source /home/rafal.wiewiora/.bashrc
cd $LS_SUBCWD
python dist_cross_tica.py kinetic 10 3
