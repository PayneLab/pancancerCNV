#!/bin/bash

#SBATCH --time=56:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH --mail-user=calebmlindgren@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash 09_00_check_for_identical.sh all_run5/
