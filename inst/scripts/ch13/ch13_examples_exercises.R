library(mqriskR)

# R Check for Example 13.3
# Rebuild the table from Example 13.1
x <- 45:50
qmat <- cbind(
  q1 = c(.011, .012, .013, .014, .015, .016),
  q2 = c(.100, .100, .100, .100, .100, .100)
)

tbl <- md_table(x = x, qxj = qmat, radix = 1000)

# Person age 46
x0 <- 46
l0 <- tbl$ltau[tbl$x == x0]

# Joint probabilities Pr(K_x = k, J_x = j)
k_vals <- 0:(max(tbl$x) - x0)
joint_q1 <- sapply(k_vals, function(k) {
  row <- tbl[tbl$x == x0 + k, ]
  row$d1 / l0
})
joint_q2 <- sapply(k_vals, function(k) {
  row <- tbl[tbl$x == x0 + k, ]
  row$d2 / l0
})

joint_df <- data.frame(
  k = k_vals,
  cause1 = joint_q1,
  cause2 = joint_q2,
  total = joint_q1 + joint_q2
)

print(joint_df)

# Requested probabilities from Example 13.3
pr_k2_j1 <- joint_df$cause1[joint_df$k == 2]
pr_k3 <- joint_df$total[joint_df$k == 3]

c(
  `Pr(K=2, J=1)` = pr_k2_j1,
  `Pr(K=3)` = pr_k3
)

# Plot the joint and marginal probabilities
matplot(
  joint_df$k,
  cbind(joint_df$cause1, joint_df$cause2, joint_df$total),
  type = "b", pch = 1,col = 1:3,lty = 1:3,
  xlab = "k",
  ylab = "Probability",
  main = "Joint and marginal probabilities for Example 13.3"
)
legend(
  "topright",
  legend = c("Pr(K=k, J=1)", "Pr(K=k, J=2)", "Pr(K=k)"),
  lty = 1:3, pch = 1, bty = "n", col = 1:3
)



# R Check for Example 13.7
# Decrement 1: move away during the year
# Decrement 2: pass to next grade at year-end only

q1_dep <- 0.02
q2_dep <- 0.96

# Out of 100 starting students:
# 2 leave by Decrement 1 during the year, so 98 remain exposed to year-end passing.
# Of those 98, 96 pass.
q2_prime <- 96 / 98

c(
  q1_dep = q1_dep,
  q2_dep = q2_dep,
  q2_prime = q2_prime
)

# Simple comparison plot
vals <- c(q1_dep, q2_dep, q2_prime)
names(vals) <- c("q^(1)", "q^(2)", "q'^(2)")
barplot(
  vals,
  ylab = "Probability",
  main = "Dependent vs independent probabilities in Example 13.7"
)

# A small sensitivity check:
# Suppose Decrement 1 varied from 0 to 0.10 while observed q^(2)=0.96 stayed fixed.
q1_grid <- seq(0, 0.10, by = 0.005)
q2prime_grid <- 0.96 / (1 - q1_grid)

plot(
  q1_grid, q2prime_grid, type = "l",
  xlab = "Dependent probability of Decrement 1",
  ylab = "Implied single-decrement probability q'^(2)",
  main = "Sensitivity of q'^(2) to competing Decrement 1"
)



# R Check for Example 13.10

library(mqriskR)

# Example-style numerical comparison for the summary table in Example 13.10

q1_prime <- 0.20
q2_prime <- 0.10

# SUDD: independent -> dependent
dep_sudd <- qx_dep_sudd(q1prime = q1_prime, q2prime = q2_prime)

# Constant force: independent -> dependent
dep_cf <- qx_dep_cf(c(q1_prime, q2_prime))

# MUDD / CF: dependent -> independent
ind_from_mudd_on_sudd <- qxprime_mudd(unname(dep_sudd))
ind_from_mudd_on_cf <- qxprime_mudd(unname(dep_cf))

# SUDD inverse: dependent -> independent
ind_from_sudd <- qxprime_sudd(q1 = dep_sudd["q1"], q2 = dep_sudd["q2"])

out <- rbind(
  `Independent input` = c(q1_prime, q2_prime),
  `Dependent under SUDD` = unname(dep_sudd),
  `Dependent under CF` = unname(dep_cf),
  `Recovered by MUDD from SUDD-dependent` = unname(ind_from_mudd_on_sudd),
  `Recovered by MUDD from CF-dependent` = unname(ind_from_mudd_on_cf),
  `Recovered by SUDD from SUDD-dependent` = unname(ind_from_sudd)
)

colnames(out) <- c("Cause 1", "Cause 2")
print(round(out, 6))




# Exercise 13-1
# Reconstruct the school table numerically

radix <- 1000

# Year 1
q_fail_1 <- 0.40
q_withdraw_1 <- 0.20
p_complete_1 <- 1 - q_fail_1 - q_withdraw_1

start_yr2 <- radix * p_complete_1

# Year 2
# Let F = failures in year 2, C = complete year 2
# F + C = 0.70 * start_yr2
# F = 0.40 * C
C2 <- (0.70 * start_yr2) / 1.40
F2 <- 0.40 * C2
W2 <- start_yr2 - F2 - C2

# Year 3
start_yr3 <- C2
p_complete_3 <- 0.60
C3 <- start_yr3 * p_complete_3

# ten times as many complete Year 2 as fail during Year 3
F3 <- C2 / 10
W3 <- start_yr3 - C3 - F3

total_withdraw <- radix * q_withdraw_1 + W2 + W3
prob_withdraw <- total_withdraw / radix

out <- c(
  start_yr2 = start_yr2,
  C2 = C2,
  F2 = F2,
  W2 = W2,
  start_yr3 = start_yr3,
  C3 = C3,
  F3 = F3,
  W3 = W3,
  total_withdraw = total_withdraw,
  prob_withdraw = prob_withdraw
)

print(round(out, 4))

barplot(
  c(Year1 = radix * q_withdraw_1, Year2 = W2, Year3 = W3),
  ylab = "Expected number of withdrawals",
  main = "Exercise 13-1: Withdrawals by year"
)



# Exercise 13-3

t_grid <- seq(0, 12, by = 0.1)

# Cause 1: mortality under uniform model with omega = 100, current age 50
tpx1 <- function(t) {
  out <- 1 - t / 50
  out[t > 50] <- 0
  out
}

# Cause 2: leaving academic employment with constant force 0.05
tpx2 <- function(t) exp(-0.05 * t)

# Total survival in multiple-decrement model
tpx_tau <- function(t) tpx1(t) * tpx2(t)

p5 <- tpx_tau(5)
p10 <- tpx_tau(10)
ans <- p5 - p10

c(
  p5 = p5,
  p10 = p10,
  probability = ans
)

plot(
  t_grid, tpx_tau(t_grid), type = "l",
  xlab = "t",
  ylab = expression(t*p[50]^{(tau)}),
  main = "Exercise 13-3: Survival in academic employment"
)
abline(v = c(5, 10), lty = 2)



# Exercise 13-5

t_grid <- seq(0, 49.9, by = 0.1)

tpx1_prime <- function(t) (75 - t) / 75
tpx2_prime <- function(t) ((50 - t) / 50)^2
tpx_tau <- function(t) tpx1_prime(t) * tpx2_prime(t)

mu1 <- function(t) 1 / (75 - t)
mu2 <- function(t) 2 / (50 - t)

fTJ1 <- function(t) tpx_tau(t) * mu1(t)
fTJ2 <- function(t) tpx_tau(t) * mu2(t)

# One-year decrement probabilities
qx1 <- integrate(fTJ1, lower = 0, upper = 1)$value
qx2 <- integrate(fTJ2, lower = 0, upper = 1)$value
qxtau <- qx1 + qx2

# Single-decrement probabilities
qx1_prime <- 1 - tpx1_prime(1)
qx2_prime <- 1 - tpx2_prime(1)

# Ultimate probabilities by cause
fJ1 <- integrate(fTJ1, lower = 0, upper = 50)$value
fJ2 <- integrate(fTJ2, lower = 0, upper = 50)$value

# Conditional cause probabilities at t = 1
fJ_given_T_1 <- mu1(1) / (mu1(1) + mu2(1))
fJ_given_T_2 <- mu2(1) / (mu1(1) + mu2(1))

out <- c(
  qx1 = qx1,
  qx2 = qx2,
  qxtau = qxtau,
  qx1_prime = qx1_prime,
  qx2_prime = qx2_prime,
  fJ1 = fJ1,
  fJ2 = fJ2,
  fJ_given_T_1 = fJ_given_T_1,
  fJ_given_T_2 = fJ_given_T_2
)

print(round(out, 6))

matplot(
  t_grid,
  cbind(fTJ1(t_grid), fTJ2(t_grid)),
  type = "l", lty = 1,col=1:2,
  xlab = "t",
  ylab = "Density",
  main = "Exercise 13-5: Cause-specific joint densities"
)
legend(
  "topright",
  legend = c(expression(f[T,J](t,1)), expression(f[T,J](t,2))),
  lty = 1, bty = "n",col=1:2,
)



# Exercise 13-6

radix <- 1000
total_leaving_first_year <- 48
survive_first_year <- 1 - total_leaving_first_year / radix

# total annual force
mu_tau <- -log(survive_first_year)

# given mu^(2) = 0.04
mu1 <- mu_tau - 0.04

# annual dependent probability for Cause 1
qx1 <- (mu1 / mu_tau) * (1 - exp(-mu_tau))

# expected survivors entering fourth year
l3 <- radix * survive_first_year^3

expected_failures_fourth_year <- l3 * qx1

c(
  survive_first_year = survive_first_year,
  mu_tau = mu_tau,
  mu1 = mu1,
  qx1 = qx1,
  l3 = l3,
  expected_failures_fourth_year = expected_failures_fourth_year
)

years <- 0:4
survivors <- radix * survive_first_year^years

plot(
  years, survivors, type = "b",
  xlab = "Completed years",
  ylab = "Expected survivors",
  main = "Exercise 13-6: Expected survivors by year"
)





# Exercise 13-13

library(mqriskR)

q1_prime <- c(0.10, 0.20, 0.20)
q2_prime <- c(0.25, 0.20, 0.10)

dep <- t(mapply(
  function(a, b) qx_dep_sudd(q1prime = a, q2prime = b),
  q1_prime, q2_prime
))

out <- data.frame(
  year = 1:3,
  q1_prime = q1_prime,
  q2_prime = q2_prime,
  q1 = dep[, "q1"],
  q2 = dep[, "q2"]
)

print(round(out, 4))

matplot(
  out$year,
  cbind(out$q1_prime, out$q2_prime, out$q1, out$q2),
  type = "b", pch = 1, col=1:4,
  xlab = "Year",
  ylab = "Probability",
  main = "Exercise 13-13: Independent vs dependent probabilities"
)
legend(
  "topright",
  legend = c("q1 prime", "q2 prime", "q1", "q2"),
  lty = 1, pch = 1, bty = "n", col=1:4
)



# Exercise 13-17

radix <- 1000
q1_prime <- 0.100
q2_prime <- 0.125
t_event <- 0.70

# Decrement 1 acts alone up to time 0.70 under UDD
q1_upto_event <- t_event * q1_prime
survivors_at_event <- radix * (1 - q1_upto_event)

# Decrement 2 occurs at time 0.70
d2 <- survivors_at_event * q2_prime
q2_dep <- d2 / radix

c(
  q1_upto_event = q1_upto_event,
  survivors_at_event = survivors_at_event,
  d2 = d2,
  q2_dep = q2_dep
)

barplot(
  c(
    start = radix,
    survivors_at_event = survivors_at_event,
    decrement2 = d2
  ),
  main = "Exercise 13-17: Cohort tracking to the point decrement",
  ylab = "Count"
)
