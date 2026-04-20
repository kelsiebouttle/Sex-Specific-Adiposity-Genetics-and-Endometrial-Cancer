source("/working/lab_tracyo/kelsieB/cpassoc/CPASSOC/FunctionSet.R")

# ── EC vs Obesity ──────────────────────────────────────────────────────────────

Dt <- read.table("merged_zscores.txt", header = TRUE)

pruned  <- read.table("/working/lab_tracyo/kelsieB/reference_files/cpassoc/EUR_1000G_Merged.prune.in", header = FALSE)
Dtcorr  <- Dt[Dt$SNP %in% pruned$V1, ]
Dtcorr  <- Dtcorr[apply(abs(Dtcorr[, c("Z_EC", "Z_Obesity")]) <= 5, 1, all), ]

corMatrix      <- diag(2)
corMatrix[1,2] <- corMatrix[2,1] <- cor(Dtcorr$Z_EC, Dtcorr$Z_Obesity, use = "complete.obs")

Sumstat    <- subset(Dt, select = c("Z_EC", "Z_Obesity"))
Samplesize <- c(EC = 226607, Obesity = 305573)

folds    <- cut(seq_len(nrow(Sumstat)), breaks = 20, labels = FALSE)
para     <- EstimateGamma(N = 1E4, Samplesize, corMatrix)
res.shet <- list()

for (i in 1:20) {
  idx       <- which(folds == i)
  Test.shet <- SHet(Sumstat[idx, ], Samplesize, corMatrix, correct = 1, isAllpossible = TRUE)
  p.shet    <- pgamma(Test.shet - para[3], shape = para[1], scale = para[2], lower.tail = FALSE)
  res.shet[[i]] <- data.frame(SNP = Dt$SNP[idx], p.Shet = p.shet)
}

res.shet <- do.call(rbind, res.shet)
write.table(res.shet, "CPASSOC_SHet_EC_Obesity.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# ── Multi-trait ────────────────────────────────────────────────────────────────

Dt <- read.table("ALL_zscores_merged_FINAL.txt", header = TRUE)

Dtcorr  <- Dt[Dt$MarkerName %in% pruned$V1, ]
z_cols  <- c("Z_EC", "Z_BMI", "Z_WHR", "Z_BFP", "Z_VAT", "Z_ASAT", "Z_GFAT")
Dtcorr  <- Dtcorr[apply(abs(Dtcorr[, z_cols]) <= 1.96, 1, all), ]

corMatrix <- diag(7)
for (i in 1:6)
  for (j in (i+1):7)
    corMatrix[i,j] <- corMatrix[j,i] <- cor(Dtcorr[, i+1], Dtcorr[, j+1])

Sumstat    <- subset(Dt, select = z_cols)
Samplesize <- c(EC = 226607, BMI = 262817, WHRadjBMI = 262759, BFP = 154337,
                VAT = 20038, ASAT = 20038, GFAT = 20038)

folds    <- cut(seq_len(nrow(Sumstat)), breaks = 20, labels = FALSE)
para     <- EstimateGamma(N = 1E4, Samplesize, corMatrix)
res.shet <- list()

for (i in 1:20) {
  idx       <- which(folds == i)
  Test.shet <- SHet(Sumstat[idx, ], Samplesize, corMatrix, correct = 1, isAllpossible = TRUE)
  p.shet    <- pgamma(Test.shet - para[3], shape = para[1], scale = para[2], lower.tail = FALSE)
  res.shet[[i]] <- data.frame(MarkerName = Dt$MarkerName[idx], p.Shet = p.shet)
}

res.shet <- do.call(rbind, res.shet)
write.table(res.shet, "CPASSOC_SHet_Multitrait.txt", row.names = FALSE, quote = FALSE, sep = "\t")
