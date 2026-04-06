

# Example 17.7 inputs
library(mqriskR)
G <- 19250
expense <- 240
i <- 0.06

q <- c(0.015, 0.017, 0.019, 0.021, 0.024)
p <- 1 - q

# Given reserves (Table 17.1)
V <- c(0, 2500, 4000, 5000, 4000, 0)

# Pre-contract expense
Pr0 <- -5000

# Profit vector
Pr <- numeric(6)
Pr[1] <- Pr0

for (t in 0:4) {
  inflow <- V[t+1] + G - expense
  accumulation <- inflow * (1 + i)

  outflow <- 1000000 * q[t+1] + V[t+2] * p[t+1]

  Pr[t+2] <- accumulation - outflow
}

round(Pr, 2)

# Or use Package function
Pr_pkg<- Pr_vector_disc(
  V = V,
  G = 19250,
  i = 0.06,
  r = 0,
  e = 240,
  q1 = q,
  b1 = 1000000,
  pre_contract_expense = 5000)

round(Pr_pkg, 2)

# Total expected profit (undiscounted)
total_profit <- sum(Pr)
round(total_profit, 2)

# Plot profit vector
plot(
  0:5, Pr,
  type = "b", pch = 19,
  main = "Example 17.7: Profit Vector",
  xlab = "Policy year",
  ylab = "Expected profit"
)
abline(h = 0, lty = 2)





# Example 17.8
library(mqriskR)
Pr <- c(-5000.00, 2688.10, 1868.60, 485.60, 534.60, 390.60)
p_tau <- c(0.98500, 0.98300, 0.98100, 0.97900, 0.97600)

Pi <- Pi_signature(Pr, p_tau = p_tau)
npv <- NPV_profit(Pi, r = 0.10)
apv_gp <- APV_gross_premiums(
  G = rep(19250, 5),
  r = 0.10,
  p_tau = p_tau
)
pm <- profit_margin(npv, apv_gp)

round(Pi, 2)
round(npv, 2)
round(apv_gp, 2)
round(pm, 5)

plot(
  0:5, Pi,
  type = "b", pch = 19,
  main = "Example 17.8: Profit Signature",
  xlab = "Policy year",
  ylab = "Profit per policy issued"
)
abline(h = 0, lty = 2)







# Example 17.9
library(mqriskR)
V <- c(0, 0, 0, 0, 0, 0)
q1 <- c(0.015, 0.017, 0.019, 0.021, 0.024)
p_tau <- 1 - q1

Pr <- Pr_vector_disc(
  V = V,
  G = 19279,
  i = 0.06,
  r = 0,
  e = 240,
  q1 = q1,
  b1 = 1000000,
  pre_contract_expense = 5000
)

Pi <- Pi_signature(Pr, p_tau = p_tau)
npv <- NPV_profit(Pi, r = 0.10)

round(Pr, 2)
round(Pi, 2)
round(npv, 2)

plot(
  0:5, Pr,
  type = "b", pch = 19,
  main = "Example 17.9: Profit Vector with No Reserves",
  xlab = "Policy year",
  ylab = "Expected profit"
)
abline(h = 0, lty = 2)

plot(
  0:5, Pi,
  type = "b", pch = 19,
  main = "Example 17.9: Profit Signature with No Reserves",
  xlab = "Policy year",
  ylab = "Profit per policy issued"
)
abline(h = 0, lty = 2)





# Example 17.11
library(mqriskR)
Vt <- 6559.89
Vt1 <- 7460.99
G <- 2264.12
B <- 100000

i <- 0.04
i_actual <- 0.06
q <- 0.00812
q_actual <- 0.00700
r <- 0.06
r_actual <- 0.08

# Expected profit
Pr_exp <- (Vt + G * (1 - r)) * (1 + i) -
  B * q - Vt1 * (1 - q)

# Actual profit
Pr_act <- (Vt + G * (1 - r_actual)) * (1 + i_actual) -
  B * q_actual - Vt1 * (1 - q_actual)

GT <- Pr_act - Pr_exp

# Gain by source, calculating interest first
GI <- (Vt + G * (1 - r)) * (i_actual - i)
GM <- (B - Vt1) * (q - q_actual)
GE <- GT - GI - GM

round(c(
  expected_profit = Pr_exp,
  actual_profit = Pr_act,
  total_gain = GT,
  gain_interest = GI,
  gain_mortality = GM,
  gain_expense = GE
), 2)

barplot(
  c(GI, GM, GE, GT),
  names.arg = c("Interest", "Mortality", "Expense", "Total"),
  main = "Example 17.11: Gain by Source",
  ylab = "Amount"
)
abline(h = 0, lty = 2)






# Example 17.12

library(mqriskR)

B_basic <- 50000
PA20 <- 10000

G <- 280
expense <- 50
i <- 0.07
q50 <- 0.00356
p50 <- 1 - q50

A30 <- 0.07570
A50 <- 0.19908
A51 <- 0.20820

V20_basic <- B_basic * (A50 - A30) / (1 - A30)
V21_basic <- B_basic * (A51 - A30) / (1 - A30)

V20_pua <- PA20 * A50
V21_pua <- PA20 * A51

# This reproduces the printed solution
Pr21_book <- (V20_basic + V20_pua + G - expense) * (1 + i) -
  q50 * B_basic -
  p50 * (V21_basic + V21_pua)

Div21_book <- 0.85 * Pr21_book
X21_book <- Div21_book / A51
p_compound_book <- X21_book / (B_basic + PA20)


round(c(
  V20_basic = V20_basic,
  V21_basic = V21_basic,
  V20_pua = V20_pua,
  V21_pua = V21_pua,
  profit_21 = Pr21_book,
  dividend_21 = Div21_book,
  paid_up_addition_21 = X21_book,
  compound_bonus_rate = p_compound_book
), 5)









# Exercise 17-1
library(mqriskR)
Pi <- c(-300, 100, 90, 80, 70)
r <- 0.10

npv <- NPV_profit(Pi, r = r)
irr <- IRR_profit(Pi, interval = c(0.001, 0.20))
payback <- discounted_payback_period(Pi, r = r)

# survival probabilities for APV of gross premiums
q <- c(0.05, 0.06, 0.07, 0.08)
p_tau <- 1 - q

apv_gp <- APV_gross_premiums(
  G = rep(1000, 4),
  r = r,
  p_tau = p_tau
)

pm <- profit_margin(npv, apv_gp)

round(c(
  NPV = npv,
  IRR = irr,
  APV_GP = apv_gp,
  profit_margin = pm
), 5)

payback

# Plot partial NPV path
partial <- NPV_partial(Pi, r = r)

plot(
  0:4, partial,
  type = "b", pch = 19,
  main = "Exercise 17-1: Partial Net Present Values",
  xlab = "t",
  ylab = "NPV(t)"
)
abline(h = 0, lty = 2)




library(mqriskR)

# Exercise 17-2
V <- c(0, 800, 800, 0)
q1 <- c(0.015, 0.017, 0.019)
p_tau <- 1 - q1

Pr <- Pr_vector_disc(
  V = V,
  G = 17500,
  i = 0.05,
  r = 0,
  e = 100,
  q1 = q1,
  b1 = 1000000,
  pre_contract_expense = 3000
)

Pi <- Pi_signature(Pr, p_tau = p_tau)
npv <- NPV_profit(Pi, r = 0.10)
irr <- IRR_profit(Pi, interval = c(0.001, 0.50))
payback <- discounted_payback_period(Pi, r = 0.10)

apv_gp <- APV_gross_premiums(
  G = rep(17500, 3),
  r = 0.10,
  p_tau = p_tau
)

pm <- profit_margin(npv, apv_gp)

round(Pr, 5)
round(Pi, 5)
round(c(
  NPV = npv,
  IRR = irr,
  APV_GP = apv_gp,
  profit_margin = pm
), 5)

payback

# Plot profit vector and profit signature
plot(
  0:3, Pr,
  type = "b", pch = 19,
  main = "Exercise 17-2: Profit Vector",
  xlab = "t",
  ylab = "Pr_t"
)
abline(h = 0, lty = 2)

plot(
  0:3, Pi,
  type = "b", pch = 19,
  main = "Exercise 17-2: Profit Signature",
  xlab = "t",
  ylab = "Pi_t"
)
abline(h = 0, lty = 2)






# Exercise 17-3

library(mqriskR)

Pr <- c(-300, 75, 250)
p40 <- 1 - 0.00278
target_pm <- -0.10

Pi <- c(
  Pr[1],
  Pr[2],
  Pr[3] * p40
)

f <- function(r) {
  npv <- NPV_profit(Pi, r = r)
  apv_gp <- APV_gross_premiums(
    G = c(100, 100),
    r = r,
    p_tau = c(p40)
  )
  profit_margin(npv, apv_gp) - target_pm
}

root <- uniroot(f, interval = c(0.001, 0.30))$root
round(root, 5)

# Plot objective function
r_grid <- seq(0.001, 0.30, length.out = 400)
y_grid <- sapply(r_grid, f)

plot(
  r_grid, y_grid,
  type = "l",
  main = "Exercise 17-3: Solving for the Risk Discount Rate",
  xlab = "r",
  ylab = "profit margin(r) - target"
)
abline(h = 0, lty = 2)
abline(v = root, lty = 3)








# Exercise 17-7
library(mqriskR)
q <- 0.013 + 0.001 * (0:4)
p <- 1 - q
i <- 0.05

# Given no-reserve profit vector beyond Pr0
profit_no_reserve <- c(450, 230, 100, -50, -220)

# Backward zeroization using the textbook logic
Vz <- numeric(6)
Vz[6] <- 0  # terminal reserve

# t = 4
Vz[5] <- 220 / 1.05

# t = 3
Vz[4] <- (50 + Vz[5] * p[4]) / 1.05

# t = 2
Vz[3] <- (-100 + Vz[4] * p[3]) / 1.05

# t = 1
Vz[2] <- (-230 + Vz[3] * p[2]) / 1.05
Vz[2] <- max(Vz[2], 0)

# t = 0
Vz[1] <- (-450 + Vz[2] * p[1]) / 1.05
Vz[1] <- max(Vz[1], 0)

names(Vz) <- paste0("V", 0:5)
round(Vz, 2)

# Revised profit vector
Pr_new <- numeric(6)
Pr_new[1] <- -300

# Pr1 uses V0 = V1 = 0
Pr_new[2] <- 450

# Pr2 uses reset reserve V1 = 0 and computed V2
Pr_new[3] <- 230 - Vz[3] * p[2]

# By construction, later profits are zero
Pr_new[4:6] <- 0

names(Pr_new) <- paste0("Pr", 0:5)
round(Pr_new, 2)

# Plot zeroized reserves
plot(
  0:5, Vz,
  type = "b", pch = 19,
  main = "Exercise 17-7: Zeroized Reserves",
  xlab = "t",
  ylab = expression(t * V^Z)
)
abline(h = 0, lty = 2)

# Plot revised profit vector
plot(
  0:5, Pr_new,
  type = "b", pch = 19,
  main = "Exercise 17-7: Revised Profit Vector",
  xlab = "t",
  ylab = expression(Pr[t])
)
abline(h = 0, lty = 2)









# Exercise 17-14
library(mqriskR)
q_d <- c(.00390, .00420, .00454, .00491, .00532,
         .00583, .00632, .00687, .00747, .00812)

q_w <- c(.04990, .04989, .04988, .04987, .04986,
         .02991, .02990, .02990, .02989, .02988)

p_tau <- c(.94620, .94591, .94558, .94522, .94482,
           .96426, .96378, .96323, .96264, .96200)

tp_tau <- c(.94620, .89502, .84631, .79995, .75581,
            .72880, .70240, .67657, .65130, .62655)

V <- c(0.00, 0.00, 0.00, 0.00, 1415.00, 3035.00,
       4738.00, 6476.00, 8251.00, 10062.00, 11908.00)

CV <- c(0.00, 0.00, 0.00, 0.00, 5708.00, 7310.40,
        8912.80, 10515.20, 12117.60, 13720.00)

G <- rep(2264.12, 10)
r_exp <- c(0.06, rep(0.03, 9))

Pr <- Pr_vector_disc(
  V = V,
  G = G,
  i = 0.04,
  r = r_exp,
  e = 0,
  q1 = q_d,
  q2 = q_w,
  b1 = 100000,
  b2 = CV,
  p_tau = p_tau,
  pre_contract_expense = 5000
)

Pi <- c(Pr[1], Pr[-1] * c(1, tp_tau[-length(tp_tau)]))

npv <- NPV_profit(Pi, r = 0.10)
apv_gp <- APV_gross_premiums(G = G, r = 0.10, p_tau = p_tau)
pm <- profit_margin(npv, apv_gp)

round(Pr, 5)
round(Pi, 5)
round(c(
  NPV = npv,
  APV_GP = apv_gp,
  profit_margin = pm
), 5)

# Plot profit signature
plot(
  0:10, Pi,
  type = "b", pch = 19,
  main = "Exercise 17-14: Profit Signature",
  xlab = "Policy year",
  ylab = "Pi_t"
)
abline(h = 0, lty = 2)




