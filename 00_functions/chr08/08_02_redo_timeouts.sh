#!/bin/bash

#SBATCH --time=4:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2000M   # memory per CPU core
#SBATCH --mail-user=calebmlindgren@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output/slrm_%a.out 
#SBATCH --array=120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,4668,4669,4670,4671,4672,4673,4674,4675,4676,4677,4678,4679

python 08_00_run_ith_test.py $SLURM_ARRAY_TASK_ID &
python 08_00_run_ith_test.py $((SLURM_ARRAY_TASK_ID + 5000)) &
wait

