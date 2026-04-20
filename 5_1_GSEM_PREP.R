library(data.table)
library(tidyverse)
library(GenomicSEM)

options(datatable.fread.datatable = FALSE)
setwd("path")

trait.names <- c("EC", "BMI", "WHR", "BFP", "VAT", "ASAT", "GFAT")

# LDSC
LDSCoutput <- ldsc(
  traits = c("inputs/EC.sumstats.gz",
             "inputs/BMI.sumstats.gz",
             "inputs/WHR_munge.sumstats.gz",
             "inputs/BFP_munge.sumstats.gz",
             "inputs/VAT_munge.sumstats.gz",
             "inputs/ASAT_munge.sumstats.gz",
             "inputs/GFAT_munge.sumstats.gz"),
  sample.prev     = c(0.5, NA, NA, NA, NA, NA, NA),
  population.prev = c(0.034, NA, NA, NA, NA, NA, NA),
  ld  = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
  wld = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
  trait.names = trait.names
)
dir.create("output", showWarnings = FALSE)
save(LDSCoutput, file = "output/LDSCoutput.RData")

# Sumstats
p_sumstats <- sumstats(
  files = c("inputs/clean_2024_EC_sumstats_with_N.txt",
            "inputs/1A_LDSC_filtered_bmi.giant-ukbb.meta-analysis.females.23May2018.txt",
            "inputs/3A_LDSC_filtered_whr.giant-ukbb.meta-analysis.females.23May2018.txt",
            "inputs/BFP_females_cleaned_hg19.tsv",
            "inputs/6A_filtered_vat_female_stats_N",
            "inputs/7A_filtered_asat_female_stats_N",
            "inputs/8A_filtered_gfat_female_stats_N"),
  ref         = "/working/lab_tracyo/kelsieB/reference_files/reference.1000G.maf.0.005.txt",
  trait.names = trait.names,
  se.logit    = c(TRUE,  rep(FALSE, 6)),
  OLS         = c(FALSE, rep(TRUE,  6)),
  info.filter = 0.6,
  maf.filter  = 0.01
)
save(p_sumstats, file = "output/ALLTraits_EC_Sumstats.RData")
