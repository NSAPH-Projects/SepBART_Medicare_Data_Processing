#!/bin/bash
#
#SBATCH -p fasse # partition (queue)
#SBATCH -c 40 # number of cores
#SBATCH --mem 148000 # memory pool for all cores
#SBATCH -t 1-08:00 # time (D-HH:MM)

singularity exec --cleanenv --env R_LIBS_USER=$HOME/R/ifxrstudio/RELEASE_3_14 /n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_14.sif Rscript run_sepbart_slurm_2000.R