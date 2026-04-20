library(GenomicSEM)

args <- commandArgs(trailingOnly = TRUE)
START_ROW <- if (length(args) >= 3) as.numeric(args[1]) else NULL
END_ROW   <- if (length(args) >= 3) as.numeric(args[2]) else NULL
BATCH_ID  <- if (length(args) >= 3) as.character(args[3]) else "full"

setwd("/working/lab_tracyo/kelsieB/publication/obesity_ec/genomicSEM/commonfactor/male/")

load("outputs/Male_LDSCoutput.RData")
load("outputs/ALLTraits_Male_Sumstats.RData")

obesity_sumstats <- p_sumstats

if (!is.null(START_ROW)) {
  END_ROW <- min(END_ROW, nrow(obesity_sumstats))
  obesity_sumstats <- obesity_sumstats[START_ROW:END_ROW, ]
}

model <- "F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
          WHR ~~ VAT
          WHR ~~ GFAT
          F1 ~ SNP"

obesity_GWAS <- userGWAS(covstruc = LDSCoutput,
                         SNPs     = obesity_sumstats,
                         model    = model,
                         sub      = "F1~SNP",
                         parallel = FALSE,
                         cores    = 1,
                         Q_SNP    = TRUE)

dir.create("male_results", showWarnings = FALSE)

out_stem <- if (BATCH_ID != "full") paste0("male_results/male_CF_GWAS_batch", BATCH_ID) else "male_results/male_CF_GWAS"

write.table(obesity_GWAS[[1]], file = paste0(out_stem, ".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
save(obesity_GWAS, file = paste0(out_stem, ".RData"))
