#!/bin/bash

#SBATCH --time=6:00:00   # walltime
#SBATCH --ntasks=11   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4000M   # memory per CPU core
#SBATCH --mail-user=calebmlindgren@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o raw_output/slrm_%a.out 
#SBATCH --array 0-9998:11 # 9999 factors to 3*3*11*101

for i in {0..10}; do # run 11 processes in parallel
     python 09_00_run_ith_test.py $((SLURM_ARRAY_TASK_ID + i)) &
done
wait
