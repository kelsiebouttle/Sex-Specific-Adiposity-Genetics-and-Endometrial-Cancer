#!/bin/bash

PHESANT="/working/lab_tracyo/kelsieB/PRS/PHESANT"
WORKDIR="/working/lab_tracyo/kelsieB/PRS"

echo "Submitting 60 parallel PHESANT jobs for FEMALE PRS in combined population..."

for PART in {1..60}; do
cat > ${WORKDIR}/phesant_femalePRS_p${PART}.sh << INNEREOF
#!/bin/bash
#PBS -N femalePRS_p${PART}
#PBS -l ncpus=4,mem=16gb,walltime=20:00:00
#PBS -j oe

cd ${WORKDIR}
module load R

Rscript ${PHESANT}/WAS/phenomeScan.r \\
  --phenofile=${WORKDIR}/combined_ukb_phesant.tab \\
  --traitofinterestfile=${WORKDIR}/female_trait_phesant.csv \\
  --confounderfile=${WORKDIR}/combined_confounders_phesant.csv \\
  --variablelistfile=${PHESANT}/variable-info/outcome-info.tsv \\
  --datacodingfile=${PHESANT}/variable-info/data-coding-ordinal-info.txt \\
  --traitofinterest=obesity_prs \\
  --resDir=${WORKDIR}/phewas_femalePRS_combined/ \\
  --userId=userId \\
  --tab=TRUE \\
  --genetic=TRUE \\
  --partIdx=${PART} \\
  --numParts=60

echo "Female PRS combined part ${PART}/60 complete!"
INNEREOF

chmod +x ${WORKDIR}/phesant_femalePRS_p${PART}.sh
qsub ${WORKDIR}/phesant_femalePRS_p${PART}.sh

done

echo "Submitted 60 female PRS jobs (femalePRS_p1 through femalePRS_p60)"
echo "Results will be in: ${WORKDIR}/phewas_femalePRS_combined/"
echo "Should complete in ~4 hours (half the time of 30 parts!)"
