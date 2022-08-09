#!/bin/bash

#SBATCH --time=2:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH --mail-user=calebmlindgren@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output/slrm_%a.out 
#SBATCH --array=0-4999

python 09_00_run_ith_test.py $SLURM_ARRAY_TASK_ID &
python 09_00_run_ith_test.py $((SLURM_ARRAY_TASK_ID + 5000)) &
wait
