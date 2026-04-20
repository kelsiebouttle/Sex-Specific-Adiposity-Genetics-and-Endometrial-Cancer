library(data.table)
library(GenomicSEM)

options(datatable.fread.datatable = FALSE)
setwd("/path/to/working/dir/")

trait.names <- c("EC", "BMI", "WHR", "BFP", "VAT", "ASAT", "GFAT")

# ── Step 1: LDSC ───────────────────────────────────────────────────────────────

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
  ld  = "/path/to/ld/",
  wld = "/path/to/ld/",
  trait.names = trait.names
)
dir.create("output", showWarnings = FALSE)
save(LDSCoutput, file = "output/LDSCoutput.RData")

# ── Step 2: Sumstats ───────────────────────────────────────────────────────────

p_sumstats <- sumstats(
  files = c("inputs/EC.txt",
            "inputs/BMI.txt",
            "inputs/WHR.txt",
            "inputs/BFP.txt",
            "inputs/VAT.txt",
            "inputs/ASAT.txt",
            "inputs/GFAT.txt"),
  ref         = "/path/to/reference.1000G.maf.0.005.txt",
  trait.names = trait.names,
  se.logit    = c(TRUE,  rep(FALSE, 6)),
  OLS         = c(FALSE, rep(TRUE,  6)),
  info.filter = 0.6,
  maf.filter  = 0.01
)
save(p_sumstats, file = "output/ALLTraits_EC_Sumstats.RData")

# ── Model definition ───────────────────────────────────────────────────────────

model <- 'F1 =~ 1*BMI + WHR + BFP + VAT + ASAT + GFAT
          WHR ~~ GFAT
          WHR ~~ VAT
          F1 ~ a*SNP
          EC ~ b*F1
          EC ~ c*SNP
          ab    := a*b
          total := c + (a*b)'

sub <- c("EC ~ SNP", "ab := a*b", "total := c + (a*b)")

# ── Step 3: Batched genome-wide GWAS ──────────────────────────────────────────

args    <- commandArgs(trailingOnly = TRUE)
n_start <- as.numeric(args[1])
n_stop  <- min(as.numeric(args[2]), nrow(p_sumstats))

results <- userGWAS(covstruc        = LDSCoutput,
                    SNPs            = p_sumstats[n_start:n_stop, ],
                    estimation      = "DWLS",
                    model           = model,
                    sub             = sub,
                    parallel        = FALSE,
                    fix_measurement = FALSE,
                    GC              = "none")

write.csv(results, file = paste0("output/GSEM_mediation_", n_start, "_", n_stop, ".csv"), row.names = FALSE)

# ── Step 4: Specific SNPs ──────────────────────────────────────────────────────

target_rsids      <- fread("EC_instruments_with_proxies.txt")$SNP
filtered_sumstats <- p_sumstats[p_sumstats$SNP %in% target_rsids, ]

results_specific <- userGWAS(covstruc        = LDSCoutput,
                             SNPs            = filtered_sumstats,
                             estimation      = "DWLS",
                             model           = model,
                             sub             = sub,
                             parallel        = FALSE,
                             fix_measurement = FALSE,
                             GC              = "none")

write.csv(results_specific, file = "output/GSEM_mediation_specific_snps.csv", row.names = FALSE)
