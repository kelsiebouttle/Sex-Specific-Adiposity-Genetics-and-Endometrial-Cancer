library(GenomicSEM)

load("outputs/Male_LDSCoutput.RData")

#adjust based on residual script ouput
models <- list(
  "Base"    = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT',
  "Model_1" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT',
  "Model_2" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT
               WHR ~~ GFAT',
  "Model_3" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT
               WHR ~~ GFAT
               VAT ~~ GFAT',
  "Model_4" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT
               WHR ~~ GFAT
               VAT ~~ GFAT
               BFP ~~ GFAT',
  "Model_5" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT
               WHR ~~ GFAT
               VAT ~~ GFAT
               BFP ~~ GFAT
               BFP ~~ VAT',
  "Model_6" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT
               WHR ~~ GFAT
               VAT ~~ GFAT
               BFP ~~ GFAT
               BFP ~~ VAT
               ASAT ~~ GFAT',
  "Model_7" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT
               WHR ~~ GFAT
               VAT ~~ GFAT
               BFP ~~ GFAT
               BFP ~~ VAT
               ASAT ~~ GFAT
               BMI ~~ GFAT',
  "Model_8" = 'F1 =~ BMI + WHR + BFP + VAT + ASAT + GFAT
               WHR ~~ VAT
               WHR ~~ GFAT
               VAT ~~ GFAT
               BFP ~~ GFAT
               BFP ~~ VAT
               ASAT ~~ GFAT
               BMI ~~ GFAT
               BMI ~~ VAT'
)

fit_results <- lapply(names(models), function(model_name) {
  out <- tryCatch(
    usermodel(LDSCoutput, estimation = "DWLS", model = models[[model_name]]),
    error = function(e) NULL
  )

  if (is.null(out)) return(data.frame(Model = model_name, CFI = NA, SRMR = NA, RMSEA = NA, df = NA))

  save(out, file = paste0("outputs/Male_", model_name, "_output.RData"))

  fit <- out$modelfit
  data.frame(Model = model_name, CFI = round(fit$CFI, 4),
             SRMR = round(fit$SRMR, 4), RMSEA = round(fit$RMSEA, 4), df = fit$df)
})

fit_results <- do.call(rbind, fit_results)
write.csv(fit_results, "outputs/Male_Progressive_Fit_Results.csv", row.names = FALSE)
save(fit_results, file = "outputs/Male_Progressive_Analysis.RData")
