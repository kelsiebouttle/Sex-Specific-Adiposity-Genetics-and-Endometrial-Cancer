library(data.table)
library(tidyverse)
library(GenomicSEM)

options(datatable.fread.datatable = FALSE)
setwd("/working/lab_tracyo/kelsieB/publication/obesity_ec/genomicSEM/commonfactor/")

trait.names <- c("BMI", "WHR", "BFP", "VAT", "ASAT", "GFAT")

# LDSC
LDSCoutput <- ldsc(
  traits = c("male/m_inputs/BMI_munge.sumstats.gz",
             "male/m_inputs/WHR_munge.sumstats.gz",
             "male/m_inputs/BFP_munge.sumstats.gz",
             "male/m_inputs/VAT_munge.sumstats.gz",
             "male/m_inputs/ASAT_munge.sumstats.gz",
             "male/m_inputs/GFAT_munge.sumstats.gz"),
  sample.prev     = rep(NA, 6),
  population.prev = rep(NA, 6),
  ld  = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
  wld = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
  trait.names = trait.names
)
save(LDSCoutput, file = "male/outputs/Male_LDSCoutput.RData")

# Sumstats
p_sumstats <- sumstats(
  files = c("male/m_inputs/LDSC_filtered_bmi.giant-ukbb.meta-analysis.males.23May2018.txt",
            "male/m_inputs/LDSC_filtered_whr.giant-ukbb.meta-analysis.males.23May2018.txt",
            "male/m_inputs/BFP_males_cleaned_formatted.txt",
            "male/m_inputs/filtered_0321_vat_Male_formatted.txt",
            "male/m_inputs/filtered_0321_asat_Male_formatted.txt",
            "male/m_inputs/filtered_0321_gfat_Male_formatted.txt"),
  ref         = "/working/lab_tracyo/kelsieB/reference_files/reference.1000G.maf.0.005.txt",
  trait.names = trait.names,
  se.logit    = rep(FALSE, 6),
  OLS         = rep(TRUE, 6),
  info.filter = 0.6,
  maf.filter  = 0.01
)
save(p_sumstats, file = "male/outputs/ALLTraits_Male_Sumstats.RData")
