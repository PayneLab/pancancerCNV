#!/bin/bash

#SBATCH --time=4:00:00   # walltime
#SBATCH --ntasks=11   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4000M   # memory per CPU core
#SBATCH --array=0-9998:11 # 9999 factors to 3*3*11*101
#SBATCH --mail-user=calebmlindgren@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output_raw/slrm_%a.out 

for i in {0..10}; do # run 11 processes in parallel
     python 09_00_run_ith_test.py $((SLURM_ARRAY_TASK_ID + i)) &
done
wait
