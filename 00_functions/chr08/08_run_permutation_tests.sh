#!/bin/bash

#SBATCH --time=14:00:00   # walltime
#SBATCH --ntasks=2000   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=4096M   # memory per CPU core
#SBATCH -J "perms1"   # job name
#SBATCH --mail-user=calebmlindgren@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

jupyter nbconvert --ExecutePreprocessor.timeout=None --to notebook --inplace --execute 08_has_vs_not_has_tumor_normal_props_diffs_permutation_tests.ipynb
