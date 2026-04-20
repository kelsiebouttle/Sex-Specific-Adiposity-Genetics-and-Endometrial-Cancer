cd /working/lab_tracyo/kelsieB/PRS/PHESANT/resultsProcessing

module load R

Rscript mainCombineResults.r \
  --resDir=/working/lab_tracyo/kelsieB/PRS/phewas_malePRS_combined/ \
  --variablelistfile=../variable-info/outcome-info.tsv \
  --numParts=60
  
Rscript mainCombineResults.r \
  --resDir=/working/lab_tracyo/kelsieB/PRS/phewas_femalePRS_combined/ \
  --variablelistfile=../variable-info/outcome-info.tsv \
  --numParts=60
