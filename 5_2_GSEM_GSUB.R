#!/usr/bin/env Rscript
library(GenomicSEM)

args     <- commandArgs(trailingOnly = TRUE)
beginsub <- as.numeric(args[1])
endsub   <- as.numeric(args[2])
NameAM   <- as.character(args[3])
NameAI   <- as.character(args[4])

load("output/ALLTraits_EC_LDSCoutput.RData")
load("Final_ALLTraits_EC_Sumstats.RData")

model <- 'AM =~ 1*BMI + WHR + BFP + VAT + ASAT + GFAT + NA*EC
          AI =~ NA*EC
          AM ~ SNP
          AI ~ SNP
          AM ~~ 1*AM
          AI ~~ 1*AI
          AM ~~ 0*AI
          BMI ~~ BMI
          WHR ~~ WHR
          BFP ~~ BFP
          VAT ~~ VAT
          ASAT ~~ ASAT
          GFAT ~~ GFAT
          EC ~~ 0*EC
          WHR ~~ VAT
          WHR ~~ ASAT
          SNP ~~ SNP'

outputGWAS <- userGWAS(covstruc = LDSCoutput, SNPs = p_sumstats[beginsub:endsub, ],
                       estimation = "DWLS", model = model, sub = c("AM~SNP", "AI~SNP"),
                       fix_measurement = FALSE, toler = 1e-60, cores = 1, parallel = FALSE)

d.frame <- as.data.frame(do.call(rbind, outputGWAS))
snp_efx <- d.frame[d.frame$rhs == "SNP", ]
saveRDS(snp_efx[snp_efx$lhs == "AM", ], file = paste0(NameAM, ".RData"))
saveRDS(snp_efx[snp_efx$lhs == "AI", ], file = paste0(NameAI, ".RData"))
