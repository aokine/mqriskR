# Chapter 15 R Check: Example 15.7
library(mqriskR)

# Forward rates from Table 15.11:
# f_{0,5}, f_{1,4}, f_{2,3}, f_{3,2}, f_{4,1}
f <- c(0.04, 0.05, 0.06, 0.07, 0.08)

# Mortality rates q60,...,q64
qx <- c(0.02, 0.03, 0.04, 0.05, 0.06)

# -----------------------------
# Retrospective calculation
# -----------------------------
E_5_60 <- prod(1 - qx[1:5]) / (1 + f[1])^5
E_4_61 <- prod(1 - qx[2:5]) / (1 + f[2])^4
E_3_62 <- prod(1 - qx[3:5]) / (1 + f[3])^3
E_2_63 <- prod(1 - qx[4:5]) / (1 + f[4])^2
E_1_64 <- prod(1 - qx[5])   / (1 + f[5])^1

sddot_factor_ret <- sum(1 / c(E_5_60, E_4_61, E_3_62, E_2_63, E_1_64))
P_ret <- 10000 / sddot_factor_ret

# -----------------------------
# Prospective calculation
# Recover spot rates z1,...,z5 using Equation (15.3)
# -----------------------------
z <- numeric(5)
z[5] <- f[1]

for (n in 1:4) {
  z[n] <- ((1 + z[5])^5 / (1 + f[n + 1])^(5 - n))^(1 / n) - 1
}

# APV of benefit
benefit_apv <- 10000 * prod(1 - qx) / (1 + z[5])^5

# Premium annuity-due factor
surv_start <- c(1, cumprod(1 - qx))[1:5]
premium_disc <- c(
  1,
  (1 + z[1])^(-1),
  (1 + z[2])^(-2),
  (1 + z[3])^(-3),
  (1 + z[4])^(-4)
)
annuity_due_factor <- sum(premium_disc * surv_start)

P_pros <- benefit_apv / annuity_due_factor

cat("Retrospective premium =", round(P_ret, 2), "\n")
cat("Prospective premium   =", round(P_pros, 2), "\n")
cat("Recovered spot rates  =", round(100 * z, 4), "\n")

# Plot recovered spot rates
plot(
  1:5, 100 * z, type = "b", pch = 19,
  xlab = "Maturity n",
  ylab = "Spot rate z_n (%)",
  main = "Example 15.7: recovered spot rates"
)





# Exercise 15-2 R Check
library(mqriskR)

qx <- c(.03, .04, .05, .06, .07)

scenarios <- list(
  Scenario1 = c(.06, .07, .08, .09, .10),
  Scenario2 = rep(.06, 5),
  Scenario3 = c(.06, .05, .04, .03, .02)
)

results <- sapply(scenarios, function(i) Axn1_var(qx = qx, i = i))
print(results)

# Contribution breakdown for Scenario 3
i <- scenarios$Scenario3
vt <- vt_var(i)
surv <- c(1, cumprod(1 - qx))[1:5]
contrib <- vt * surv * qx

barplot(
  contrib,
  names.arg = 1:5,
  xlab = "Year of death",
  ylab = "Contribution to APV",
  main = "Exercise 15-2: Term insurance contributions (Scenario 3)"
)



# Exercise 15-3 R Check
library(mqriskR)

qx <- rep(0.02, 5)

i_level <- rep(0.06, 5)
i_rise  <- c(0.06, 0.07, 0.08, 0.09, 0.10)
i_fall  <- c(0.06, 0.05, 0.04, 0.03, 0.03)

annuity <- c(
  level = axn_var(qx, i_level),
  rise  = axn_var(qx, i_rise),
  fall  = axn_var(qx, i_fall)
)

pure_endowment <- c(
  level = nEx_var(qx, i_level),
  rise  = nEx_var(qx, i_rise),
  fall  = nEx_var(qx, i_fall)
)

print(annuity)
print(pure_endowment)

# Regulatory test
ratio_annuity <- annuity["fall"] / annuity["level"]
ratio_pe <- pure_endowment["fall"] / pure_endowment["level"]

cat("Annuity ratio:", ratio_annuity, "\n")
cat("Pure endowment ratio:", ratio_pe, "\n")

# Plot comparison
matplot(
  1:3,
  cbind(annuity, pure_endowment),
  type = "b", pch = 19,
  xaxt = "n",
  col = 1:2,
  xlab = "Scenario",
  ylab = "APV",
  main = "Exercise 15-3: Scenario comparison"
)
axis(1, at = 1:3, labels = c("Level", "Rising", "Falling"))
legend("right",col = 1:2,
       legend = c("Annuity", "Pure Endowment"),
       pch = 19, lty = 1)



# Exercise 15-4 R Check
library(mqriskR)

qx <- c(0.10, 0.15, 0.20, 0.25, 0.30)

i_inc <- c(0.07, 0.08, 0.09, 0.10, 0.11)
i_dec <- c(0.07, 0.06, 0.05, 0.04, 0.03)

# Annuity-certain
ann_certain_inc <- sum(vt_var(i_inc))
ann_certain_dec <- sum(vt_var(i_dec))

# Life annuity
ann_life_inc <- axn_var(qx, i_inc)
ann_life_dec <- axn_var(qx, i_dec)

# Term insurance
term_inc <- Axn1_var(qx, i_inc)
term_dec <- Axn1_var(qx, i_dec)

results <- data.frame(
  Product = c("Annuity Certain", "Life Annuity", "Term Insurance"),
  Increasing = c(ann_certain_inc, ann_life_inc, term_inc),
  Decreasing = c(ann_certain_dec, ann_life_dec, term_dec)
)

results$Ratio <- results$Decreasing / results$Increasing
print(results)

# Plot sensitivity
barplot(
  t(as.matrix(results[,2:3])),
  beside = TRUE,
  names.arg = results$Product,
  legend.text = c("Increasing", "Decreasing"),
  main = "Exercise 15-4: Interest sensitivity comparison",
  ylab = "APV"
)



# Exercise 15-14 R Check
library(mqriskR)

fn1 <- c(0.04, 0.05, 0.06, 0.07, 0.08)

z <- z_from_fn1(fn1)
print(z)

# Plot term structure
plot(
  1:length(z), 100 * z,
  type = "b", pch = 19,
  xlab = "Maturity (years)",
  ylab = "Spot rate (%)",
  main = "Exercise 15-14: Term structure from forward rates"
)


