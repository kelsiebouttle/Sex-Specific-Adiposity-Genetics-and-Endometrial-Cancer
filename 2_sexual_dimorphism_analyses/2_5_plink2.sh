#!/bin/bash
#PBS -N prs_calculation
#PBS -l ncpus=4,walltime=2:00:00,mem=16GB
#PBS -j oe

set -e

SEX=$1  # "male" or "female"
BFILE_DIR="/reference/data/UKBB_500k/versions/bfile_all_SNP"
SCORE_FILE="prs_scores/${SEX}_pgs_score.txt"

cd /working/lab_tracyo/kelsieB/PRS

for chr in $(seq 1 22); do
  plink2 \
    --bfile ${BFILE_DIR}/ukb_imp_chr${chr} \
    --score ${SCORE_FILE} 1 2 3 header \
    --out prs_scores/${SEX}_obesity_prs_chr${chr} \
    --threads 4
done

Rscript -e '
args <- commandArgs(trailingOnly = TRUE)
sex  <- args[1]
library(data.table)
scores    <- rbindlist(lapply(1:22, function(chr) fread(paste0("prs_scores/", sex, "_obesity_prs_chr", chr, ".sscore"))))
score_col <- grep("SCORE", names(scores), value = TRUE)[1]
combined  <- scores[, .(SCORE = sum(get(score_col))), by = .(FID = `#FID`, IID)]
fwrite(combined, paste0("prs_scores/", sex, "_obesity_prs.sscore"))
' --args $SEX
