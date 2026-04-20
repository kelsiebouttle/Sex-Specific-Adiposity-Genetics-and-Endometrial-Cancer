library(GenomicSEM)

simple_gencor <- function(file1, file2, sample_prev = NA, pop_prev = NA,
                          trait1_name = "Trait1", trait2_name = "Trait2",
                          file1_is_binary = FALSE, file1_sample_prev = NA, file1_pop_prev = NA) {

  sample.prev     <- c(if (file1_is_binary) file1_sample_prev else NA, sample_prev)
  population.prev <- c(if (file1_is_binary) file1_pop_prev    else NA, pop_prev)

  output_text <- capture.output({
    LDSCoutput <- ldsc(traits = c(file1, file2), sample.prev = sample.prev,
                       population.prev = population.prev,
                       ld  = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
                       wld = "/working/lab_tracyo/kelsieB/LDSC/eur_w_ld_chr/",
                       trait.names = c(trait1_name, trait2_name))
  })

  line <- output_text[grepl("Genetic Correlation between", output_text)][1]
  nums <- regmatches(line, gregexpr("-?[0-9]+\\.?[0-9]*", line))[[1]]

  if (length(nums) >= 2) {
    list(rG = as.numeric(nums[length(nums) - 1]), SE = as.numeric(nums[length(nums)]))
  } else {
    list(rG = NA, SE = NA)
  }
}

trait_info        <- read.csv("trait_list.csv")
corrected_mapping <- readRDS("corrected_traits_mapping.rds")

all_results <- do.call(rbind, lapply(1:nrow(trait_info), function(i) {
  trait_name  <- trait_info$Trait_name[i]
  trait_file  <- corrected_mapping[[trait_name]]
  sample_prev <- trait_info$Sample_prevalence[i]
  pop_prev    <- trait_info$Population_prevalence[i]

  if (is.null(trait_file)) return(NULL)

  ec <- simple_gencor("EC.sumstats.gz",                  trait_file, sample_prev, pop_prev, "EC",                 trait_name, file1_is_binary = TRUE, file1_sample_prev = 0.5, file1_pop_prev = 0.034)
  oi <- simple_gencor("Obesity_Independent.sumstats.gz", trait_file, sample_prev, pop_prev, "ObesityIndependent", trait_name)
  od <- simple_gencor("Obesity_Dependent.sumstats.gz",   trait_file, sample_prev, pop_prev, "ObesityDependent",   trait_name)

  data.frame(
    Trait    = trait_name,
    Analysis = c("Direct_EC", "Obesity_Independent", "Obesity_Dependent"),
    rG       = c(as.numeric(ec$rG), as.numeric(oi$rG), as.numeric(od$rG)),
    SE       = c(as.numeric(ec$SE), as.numeric(oi$SE), as.numeric(od$SE))
  )
}))

write.csv(all_results, "EC_obesity_genetic_correlations.csv", row.names = FALSE)
