#!/bin/bash
#PBS -N clust
#PBS -l walltime=24:00:00
#PBS -l nodes=5:ppn=25
#PBS -l pmem=4000mb

module load mro/3.5.1-foss-2017a

# call to app
MY_PARALLEL_OPTS="-N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel-${SLURM_JOBID}.log"
MY_SRUN_OPTS="-N 1 -n 1 --exclusive"
MY_EXEC="R CMD BATCH '--args {1}' lus/lusc_mccv_clust.R"

parallel $MY_PARALLEL_OPTS srun $MY_SRUN_OPTS $MY_EXEC ::: {1..25}
