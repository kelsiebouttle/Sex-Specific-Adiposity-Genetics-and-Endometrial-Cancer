library(bigsnpr)
library(bigreadr)

ld_dir   <- "ld_reference/"
out_dir  <- "pgs_weights/"
NCORES   <- 10

dir.create(out_dir,    showWarnings = FALSE)
dir.create("tmp-data", showWarnings = FALSE)

map_ldref <- readRDS(paste0(ld_dir, "map.rds"))

run_ldpred2 <- function(gwas_file, sex) {
  sumstats  <- fread2(gwas_file)
  info_snp  <- snp_match(sumstats, map_ldref, join_by_pos = FALSE)

  is_bad <- info_snp$beta_se < 0.001 | info_snp$beta_se > 0.1 |
            info_snp$af_UKBB < 0.01  | info_snp$af_UKBB > 0.99
  df_beta <- info_snp[!is_bad, ]

  tmp <- tempfile(tmpdir = "tmp-data")
  for (chr in 1:22) {
    ind.chr  <- which(df_beta$chr == chr)
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
    corr_chr <- readRDS(paste0(ld_dir, "LD_chr", chr, ".rds"))[ind.chr3, ind.chr3]

    if (chr == 1) {
      ld   <- Matrix::colSums(corr_chr^2)
      corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
    } else {
      ld   <- c(ld, Matrix::colSums(corr_chr^2))
      corr$add_columns(corr_chr, nrow(corr))
    }
  }

  ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                 chi2 = (beta / beta_se)^2,
                                 sample_size = n_eff, ncores = NCORES))

  multi_auto <- snp_ldpred2_auto(corr, df_beta,
                                 h2_init       = ldsc[["h2"]],
                                 vec_p_init    = seq_log(1e-4, 0.9, length.out = 30),
                                 burn_in       = 500,
                                 num_iter      = 500,
                                 report_step   = 20,
                                 allow_jump_sign = FALSE,
                                 shrink_corr   = 0.95,
                                 ncores        = NCORES)

  range_auto <- sapply(multi_auto, function(a) diff(range(a$corr_est)))
  keep       <- which(range_auto > (0.95 * quantile(range_auto, 0.95, na.rm = TRUE)))
  beta_auto  <- rowMeans(sapply(multi_auto[keep], function(a) a$beta_est))

  out <- data.frame(chr  = df_beta$chr,  pos  = df_beta$pos,
                    rsid = df_beta$rsid, a0   = df_beta$a0,
                    a1   = df_beta$a1,   beta = beta_auto)

  fwrite2(out, paste0(out_dir, sex, "_obesity_pgs_weights.txt"))
}

run_ldpred2("male_obesity_ldpred2_input.txt",   "male")
run_ldpred2("female_obesity_ldpred2_input.txt", "female")
