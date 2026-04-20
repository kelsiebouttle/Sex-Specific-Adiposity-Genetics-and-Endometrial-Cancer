library(coloc)

load("cleaned_obesity_ec_data.RData")

args <- commandArgs(trailingOnly = TRUE)
chr  <- as.numeric(args[1])
pos  <- as.numeric(args[2])
win  <- 100000

bmi_r <- bmi_final[bmi_final$CHR == chr & bmi_final$POS_numeric >= pos - win & bmi_final$POS_numeric <= pos + win, ]
ec_r  <- ec_final[ ec_final$chr_b37 == chr & ec_final$pos_numeric >= pos - win & ec_final$pos_numeric <= pos + win, ]

common <- intersect(bmi_r$SNP, ec_r$MarkerName)

bmi_r <- bmi_r[bmi_r$SNP %in% common, ]
ec_r  <- ec_r[ ec_r$MarkerName %in% common, ]

bmi_coloc <- list(beta = bmi_r$BETA, varbeta = bmi_r$SE^2,     snp = bmi_r$SNP,
                  position = bmi_r$POS_numeric, type = "quant", N = bmi_r$N, MAF = bmi_r$MAF)

ec_coloc  <- list(beta = ec_r$Effect, varbeta = ec_r$StdErr^2, snp = ec_r$MarkerName,
                  position = ec_r$pos_numeric,  type = "cc",    N = ec_r$TotalSampleSize,
                  s = 0.1, MAF = ec_r$MAF)

result <- coloc.abf(dataset1 = bmi_coloc, dataset2 = ec_coloc)

dir.create(paste0("CHR_", chr), showWarnings = FALSE)
save(result, file = paste0("CHR_", chr, "/chr", chr, "_pos", pos, "_coloc.RData"))
