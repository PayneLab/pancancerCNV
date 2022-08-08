#!/bin/bash

#SBATCH --time=07:00:00   # walltime
#SBATCH --ntasks=1000   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=8192M   # memory per CPU core
#SBATCH -J "perms1"   # job name
#SBATCH --mail-user=calebmlindgren@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

srun python 08_has_vs_not_has_tumor_normal_props_diffs_permutation_tests.py
