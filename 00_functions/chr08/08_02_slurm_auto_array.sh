#!/bin/bash

find indices -type f | \
    slurm-auto-array --ntasks 1 --mem 2000M --time 10:00:00 \
                     -o output/slrm_%a.out \
                     --mail-type BEGIN --mail-type END --mail-type FAIL --mail-user calebmlindgren@gmail.com \
                     --command ./08_01_run_ith_test.sh
