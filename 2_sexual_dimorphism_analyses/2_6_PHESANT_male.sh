#!/bin/bash

PHESANT="/working/lab_tracyo/kelsieB/PRS/PHESANT"
WORKDIR="/working/lab_tracyo/kelsieB/PRS"

echo "Submitting 60 parallel PHESANT jobs for MALE PRS in combined population..."

for PART in {1..60}; do
cat > ${WORKDIR}/phesant_malePRS_p${PART}.sh << INNEREOF
#!/bin/bash
#PBS -N malePRS_p${PART}
#PBS -l ncpus=4,mem=16gb,walltime=20:00:00
#PBS -j oe

cd ${WORKDIR}
module load R

Rscript ${PHESANT}/WAS/phenomeScan.r \\
  --phenofile=${WORKDIR}/combined_ukb_phesant.tab \\
  --traitofinterestfile=${WORKDIR}/male_trait_phesant.csv \\
  --confounderfile=${WORKDIR}/combined_confounders_phesant.csv \\
  --variablelistfile=${PHESANT}/variable-info/outcome-info.tsv \\
  --datacodingfile=${PHESANT}/variable-info/data-coding-ordinal-info.txt \\
  --traitofinterest=obesity_prs \\
  --resDir=${WORKDIR}/phewas_malePRS_combined/ \\
  --userId=userId \\
  --tab=TRUE \\
  --genetic=TRUE \\
  --partIdx=${PART} \\
  --numParts=60

echo "Male PRS combined part ${PART}/60 complete!"
INNEREOF

chmod +x ${WORKDIR}/phesant_malePRS_p${PART}.sh
qsub ${WORKDIR}/phesant_malePRS_p${PART}.sh

done


