=female_h2 <- 0.2124; female_se <- 0.0071
male_h2   <- 0.2211; male_se   <- 0.0073

z      <- (female_h2 - male_h2) / sqrt(female_se^2 + male_se^2)
p_diff <- 2 * pnorm(-abs(z))

cat("Z:", round(z, 3), "\n")
cat("P:", round(p_diff, 3), "\n")
