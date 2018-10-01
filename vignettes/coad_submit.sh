#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=coadclust
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4571mb


module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load R/3.4.4-X11-20180131

# call to app
R CMD BATCH vignettes/bafitcoad.R
