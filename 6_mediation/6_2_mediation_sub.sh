#!/bin/bash
#PBS -N GenomicSEM_mediation
#PBS -l mem=4gb,walltime=03:00:00,ncpus=1
#PBS -l avx2=True
#PBS -J 1-696
#PBS -o /path/to/logs/
#PBS -e /path/to/logs/

module load R
cd $PBS_O_WORKDIR

VARIABLES=($(awk -v var=$PBS_ARRAY_INDEX 'NR==var' run_gsemGWAS_array.values))
Rscript pipeline.R ${VARIABLES[0]} ${VARIABLES[1]}
