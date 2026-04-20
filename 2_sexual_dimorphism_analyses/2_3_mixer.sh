#!/bin/bash
#PBS -N mixer_cross_trait
#PBS -l walltime=48:00:00
#PBS -l mem=64GB
#PBS -l ncpus=8
#PBS -l avx2=True
#PBS -o logs/mixer_cross_trait.out
#PBS -e logs/mixer_cross_trait.err

mkdir -p logs

module load mixer/20240923
export BGMG_SHARED_LIBRARY=/software/mixer/mixer-20240923/src/build/lib/libbgmg.so

WORK_DIR="/working/lab_tracyo/kelsieB/obesity_SD_test/mixer"
REF_DIR="/working/lab_tracyo/kelsieB/reference_files/LDREF"
OUT_DIR="${WORK_DIR}/mixer_results"
FEMALE_GWAS="${WORK_DIR}/Female_Obesity_GWAS_MiXeR.txt"
MALE_GWAS="${WORK_DIR}/Male_Obesity_Common_Factor_GWAS_MiXeR.txt"
EXTRACT="${REF_DIR}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep1.snps"
BIM="${REF_DIR}/1000G.EUR.@.bim"
LD="${REF_DIR}/mixer_ld/1000G.EUR.@.ld"

mkdir -p ${OUT_DIR} && cd ${OUT_DIR}

python $MIXER_ROOT/precimed/mixer.py fit1 \
    --trait1-file ${FEMALE_GWAS} --out female_fit \
    --extract ${EXTRACT} --bim-file ${BIM} --ld-file ${LD} --chr2use 1-22

python $MIXER_ROOT/precimed/mixer.py fit1 \
    --trait1-file ${MALE_GWAS} --out male_fit \
    --extract ${EXTRACT} --bim-file ${BIM} --ld-file ${LD} --chr2use 1-22

python $MIXER_ROOT/precimed/mixer.py fit2 \
    --trait1-file ${FEMALE_GWAS} --trait2-file ${MALE_GWAS} \
    --trait1-params-file female_fit.json --trait2-params-file male_fit.json \
    --out female_male_cross \
    --extract ${EXTRACT} --bim-file ${BIM} --ld-file ${LD} --chr2use 1-22

python $MIXER_ROOT/precimed/mixer.py test2 \
    --trait1-file ${FEMALE_GWAS} --trait2-file ${MALE_GWAS} \
    --load-params-file female_male_cross.json --out female_male_cross \
    --bim-file ${BIM} --ld-file ${LD} --chr2use 1-22
