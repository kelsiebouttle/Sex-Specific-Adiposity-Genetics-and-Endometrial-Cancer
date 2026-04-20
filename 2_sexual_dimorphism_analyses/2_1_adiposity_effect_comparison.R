library(data.table)
library(dplyr)

female_data <- fread("Female_Obesity_GWAS_Short.txt")
male_data   <- fread("Male_Obesity_Common_Factor_GWAS_Short.txt")

# Genome-wide correlation
common_snps <- inner_join(
  female_data %>% select(SNP, est_f = est),
  male_data   %>% select(SNP, est_m = est),
  by = "SNP"
)
r_global <- cor(common_snps$est_f, common_snps$est_m, method = "spearman")

# SNP categories
female_snps <- fread("Female_GWS_SNPs.txt")$SNP
male_snps   <- fread("Male_GWS_SNPs.txt")$SNP
shared_snps <- intersect(female_snps, male_snps)
n_total     <- length(union(female_snps, male_snps))

# Sex dimorphism function
calculate_dimorphism <- function(snp_list, category) {
  inner_join(
    female_data %>% filter(SNP %in% snp_list) %>% select(SNP, CHR, BP, beta_f = est, se_f = SE, p_f = P),
    male_data   %>% filter(SNP %in% snp_list) %>% select(SNP, beta_m = est, se_m = SE, p_m = P),
    by = "SNP"
  ) %>%
    mutate(
      se_diff              = sqrt(se_f^2 + se_m^2 - 2 * r_global * se_f * se_m),
      t_stat               = (beta_f - beta_m) / se_diff,
      p_diff               = 2 * pnorm(-abs(t_stat)),
      significant_dimorphism = p_diff < (0.05 / n_total),
      same_direction       = sign(beta_f) == sign(beta_m),
      stronger_in          = case_when(
        same_direction & abs(beta_f) > abs(beta_m) ~ "Female",
        same_direction & abs(beta_m) > abs(beta_f) ~ "Male",
        !same_direction & p_f < p_m               ~ "Female",
        !same_direction & p_m < p_f               ~ "Male",
        TRUE                                       ~ "Unclear"
      ),
      category = category
    )
}

all_results <- bind_rows(
  calculate_dimorphism(setdiff(female_snps, shared_snps), "Female-specific"),
  calculate_dimorphism(setdiff(male_snps,   shared_snps), "Male-specific"),
  calculate_dimorphism(shared_snps,                       "Shared")
)

write.csv(all_results,                                  "sex_dimorphism_all_leadsnps.csv",      row.names = FALSE)
write.csv(filter(all_results, significant_dimorphism),  "significant_dimorphic_effects.csv",    row.names = FALSE)
