#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=clustcoad
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH --ntasks=25
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4571mb
#SBATCH --cpus-per-task=8

module load mro/3.5.1-foss-2017a

# call to app
MY_PARALLEL_OPTS="-N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel-${SLURM_JOBID}.log"
MY_EXEC="R CMD BATCH '--args {1}' coad_mccv_clust.R"

parallel $MY_PARALLEL_OPTS srun $MY_EXEC ::: {1..25}
