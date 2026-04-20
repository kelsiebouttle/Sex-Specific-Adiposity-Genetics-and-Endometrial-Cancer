#!/bin/bash
#PBS -N susiex_direct
#PBS -l walltime=4:00:00,mem=16gb,ncpus=2
#PBS -o /working/lab_tracyo/kelsieB/publication/obesity_ec/genomicSEM/SUSIEX/OUT_ERR/
#PBS -e /working/lab_tracyo/kelsieB/publication/obesity_ec/genomicSEM/SUSIEX/OUT_ERR/

module load plink/1.90b7.4

lead_snp="$1"; chr="$2"; start="$3"; end="$4"; region_id="$5"

mkdir -p sugsig_finemapping_results_direct_v2
temp_ld_dir="sugsig_temp_ld_matrices_direct_v2/${region_id}_chr${chr}"
mkdir -p "${temp_ld_dir}"

ref_file="/reference/data/UKBB_500k/versions/lab_stuartma/LD_reference/MAF_0.001/LD_ref_chr${chr}_noDup"

plink --bfile "${ref_file}" --chr ${chr} --from-bp ${start} --to-bp ${end} \
      --maf 0.005 --r square --threads 2 \
      --out "${temp_ld_dir}/chr${chr}_region_LD"

plink --bfile "${ref_file}" --chr ${chr} --from-bp ${start} --to-bp ${end} \
      --maf 0.005 --make-bed --threads 2 \
      --out "${temp_ld_dir}/chr${chr}_region_ref"

/working/lab_tracyo/kelsieB/SuSiEx/bin/SuSiEx \
    --sst_file=Direct_Effects_Clean.txt \
    --n_gwas=38430 \
    --ref_file="${temp_ld_dir}/chr${chr}_region_ref" \
    --ld_file="${temp_ld_dir}/chr${chr}_region_LD" \
    --out_dir=./sugsig_finemapping_results_direct_v2/ \
    --out_name=direct_${region_id}_${lead_snp} \
    --chr=${chr} \
    --bp=${start},${end} \
    --chr_col=2 --snp_col=1 --bp_col=3 \
    --a1_col=5  --a2_col=6 --eff_col=7 \
    --se_col=8  --pval_col=9 \
    --max_iter=1000 \
    --plink=plink
