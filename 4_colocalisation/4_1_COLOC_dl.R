library(coloc)

# ── Prepare data (run once) ───────────────────────────────────────────────────

adipose_data <- read.table("Female_Obesity_LDSC_ready.txt", header = TRUE, stringsAsFactors = FALSE)
ec_data      <- read.table("2024_EC_sumstats.txt",          header = TRUE, stringsAsFactors = FALSE)

adipose_data$POS_numeric <- as.numeric(adipose_data$BP)
ec_data$pos_numeric      <- as.numeric(ec_data$pos_b37)

adipose <- adipose_data[!is.na(adipose_data$POS_numeric) & !is.na(adipose_data$BETA) & !is.na(adipose_data$SE) & adipose_data$POS_numeric > 0, ]
ec      <- ec_data[     !is.na(ec_data$pos_numeric)      & !is.na(ec_data$Effect)    & !is.na(ec_data$StdErr) & ec_data$pos_numeric > 0, ]

adipose <- adipose[!duplicated(adipose$SNP), ]
ec      <- ec[     !duplicated(ec$MarkerName), ]

adipose$MAF <- pmin(adipose$MAF, 1 - adipose$MAF)
ec$MAF      <- pmin(ec$Freq1, 1 - ec$Freq1)

adipose_final <- adipose[!is.na(adipose$MAF) & adipose$MAF > 0 & adipose$MAF < 1, ]
ec_final      <- ec[     !is.na(ec$MAF)      & ec$MAF  > 0     & ec$MAF  < 1 & ec$Freq1 > 0 & ec$Freq1 < 1, ]

save(adipose_final, ec_final, file = "cleaned_obesity_ec_data.RData")

