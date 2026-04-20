#!/bin/bash
#PBS -l mem=10gb,walltime=01:00:00,ncpus=1
#PBS -N Coloc_Obesity
#PBS -J 1-9
#PBS -o /working/lab_tracyo/kelsieB/coloc/OBESITY_EC/logs/
#PBS -e /working/lab_tracyo/kelsieB/coloc/OBESITY_EC/logs/

cd /working/lab_tracyo/kelsieB/coloc/OBESITY_EC/

VARIABLES=($(awk -v var=$PBS_ARRAY_INDEX 'NR==var' PBS_chr_pos.txt))
module load R
Rscript coloc.R ${VARIABLES[0]} ${VARIABLES[1]}
