#!/bin/bash
#PBS -l nodes=8:ppn=28
#PBS -l walltime=48:00:00
#PBS -l pmem=4571mb

module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load R/3.4.4-X11-20180131

# call to app
MY_PARALLEL_OPTS="-N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel-${SLURM_JOBID}.log"
MY_EXEC="R CMD BATCH '--args {1}' ./lusc_mccv_clust.R"

parallel $MY_PARALLEL_OPTS srun $MY_EXEC ::: {1..25}
