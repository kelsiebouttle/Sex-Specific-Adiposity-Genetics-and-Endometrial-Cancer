library(GenomicSEM)

trait.names <- c("BMI", "WHR", "BFP", "VAT", "ASAT", "GFAT")

residual_analysis <- function(sex) {
  load(paste0("outputs/", sex, "_LDSCoutput.RData"))
  load(paste0("outputs/", sex, "_ObesityFactor_Modeloutput.RData"))

  # Observed correlations
  obs_cor <- cov2cor(LDSCoutput$S)
  rownames(obs_cor) <- colnames(obs_cor) <- trait.names

  # Factor loadings
  fl <- output$results[output$results$op == "=~", ]
  loadings <- fl$STD_All[match(trait.names, fl$rhs)]

  # Implied correlations
  L <- matrix(loadings, ncol = 1)
  implied_cor <- L %*% t(L)
  diag(implied_cor) <- 1
  rownames(implied_cor) <- colnames(implied_cor) <- trait.names

  # Residuals
  ut <- upper.tri(obs_cor)
  idx <- which(ut, arr.ind = TRUE)
  residual_df <- data.frame(
    Trait1       = trait.names[idx[, 1]],
    Trait2       = trait.names[idx[, 2]],
    Observed     = round(obs_cor[ut], 3),
    Implied      = round(implied_cor[ut], 3),
    Residual     = round((obs_cor - implied_cor)[ut], 3),
    Abs_Residual = round(abs((obs_cor - implied_cor)[ut]), 3)
  )
  residual_df <- residual_df[order(residual_df$Abs_Residual, decreasing = TRUE), ]

  write.csv(residual_df, paste0("outputs/", sex, "_ResidualTable.csv"), row.names = FALSE)
  save(obs_cor, implied_cor, residual_df, file = paste0("outputs/", sex, "_ResidualAnalysis.RData"))

  residual_df
}

male_residuals   <- residual_analysis("Male")
female_residuals <- residual_analysis("Female")
