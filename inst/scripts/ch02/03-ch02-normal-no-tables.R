# Chapter 2 Script 03: Normal probabilities without tables using pnorm/qnorm

# Standard normal probabilities
p_le_1_96 <- pnorm(1.96)
p_gt_1_96 <- 1 - pnorm(1.96)

p_abs_le_1_96 <- pnorm(1.96) - pnorm(-1.96)

z_0_975 <- qnorm(0.975)

cat("P(Z <= 1.96) =", p_le_1_96, "\n")
cat("P(Z >  1.96) =", p_gt_1_96, "\n")
cat("P(|Z| <= 1.96) =", p_abs_le_1_96, "\n")
cat("z such that P(Z <= z) = 0.975:", z_0_975, "\n")
