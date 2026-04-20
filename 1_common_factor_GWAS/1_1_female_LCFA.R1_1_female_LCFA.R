library(data.table)
library(tidyverse)
library(GenomicSEM)

options(datatable.fread.datatable = FALSE)
setwd("/working/lab_tracyo/kelsieB/publication/obesity_ec/genomicSEM/commonfactor/female/")

trait.names <- c("BMI", "WHR", "BFP", "VAT", "ASAT", "GFAT")

# LDSC
LDSCoutput <- ldsc(
  traits = c("f_inputs/BMI.sumstats.gz",
             "f_inputs/WHR_munge.sumstats.gz",
             "f_inputs/BFP_munge.sumstats.gz",
             "f_inputs/VAT_munge.sumstats.gz",
             "f_inputs/ASAT_munge.sumstats.gz",
             "f_inputs/GFAT_munge.sumstats.gz"),
  sample.prev     = rep(NA, 6),
  population.prev = rep(NA, 6),
  ld  = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
  wld = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
  trait.names = trait.names
)
save(LDSCoutput, file = "outputs/Female_LDSCoutput.RData")

# Sumstats
p_sumstats <- sumstats(
  files = c("f_inputs/1A_LDSC_filtered_bmi.giant-ukbb.meta-analysis.females.23May2018.txt",
            "f_inputs/3A_LDSC_filtered_whr.giant-ukbb.meta-analysis.females.23May2018.txt",
            "f_inputs/BFP_females_cleaned_hg19.tsv",
            "f_inputs/6A_filtered_vat_female_stats_N",
            "f_inputs/7A_filtered_asat_female_stats_N",
            "f_inputs/8A_filtered_gfat_female_stats_N"),
  ref         = "/working/lab_tracyo/kelsieB/reference_files/reference.1000G.maf.0.005.txt",
  trait.names = trait.names,
  se.logit    = rep(FALSE, 6),
  OLS         = rep(TRUE, 6),
  info.filter = 0.6,
  maf.filter  = 0.01
)
save(p_sumstats, file = "outputs/ALLTraits_Female_Sumstats.RData")
