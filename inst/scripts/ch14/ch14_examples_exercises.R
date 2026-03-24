# Chapter 14 R Check: Example 14.2
# Whole life insurance with natural-cause and accidental-cause death benefits

library(mqriskR)

delta <- 0.05
mu_nc <- 0.01
mu_ac <- 0.002
mu_tau <- mu_nc + mu_ac

t <- seq(0, 200, by = 0.01)
ptau <- exp(-mu_tau * t)

Abar_nc <- Abarxj_md(
  t = t,
  ptau = ptau,
  muj = rep(mu_nc, length(t)),
  delta = delta,
  benefit = 1000
)

Abar_ac <- Abarxj_md(
  t = t,
  ptau = ptau,
  muj = rep(mu_ac, length(t)),
  delta = delta,
  benefit = 2000
)

APV_total <- Abar_nc + Abar_ac

out <- c(
  Abar_nc = Abar_nc,
  Abar_ac = Abar_ac,
  APV_total = APV_total
)
print(out)

# Plot discounted contribution densities
density_nc <- 1000 * exp(-delta * t) * ptau * mu_nc
density_ac <- 2000 * exp(-delta * t) * ptau * mu_ac
density_total <- density_nc + density_ac

plot(
  t, density_total, type = "l", lwd = 2,
  xlab = "t", ylab = "discounted contribution density",
  main = "Example 14.2: contribution densities"
)
lines(t, density_nc, lwd = 2, lty = 2)
lines(t, density_ac, lwd = 2, lty = 3)
legend(
  "topright",
  legend = c("total", "natural cause", "accidental cause"),
  lwd = 2, lty = c(1, 2, 3), bty = "n"
)



# Chapter 14 R Check: Example 14.10
# APVs and premium using Euler-approximated state probabilities

library(mqriskR)

mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

out <- tp00_tp01_euler(
  h = 0.10,
  n = 2.00,
  mu01 = mu01,
  mu02 = mu02,
  mu10 = mu10,
  mu12 = mu12
)

j <- 10
i_nom <- 0.05
i_step <- i_nom / j
v <- 1 / (1 + i_step)

q_death_step <- 1 - exp(-0.20 * 0.10)

# Part (a): death benefit APV
death_apv <- 10000 * sum(v^(1:20) * out$tp00[1:20] * q_death_step)

# Part (b): disability benefit APV
disab_apv <- 1000 * sum(v^(1:20) * out$tp01[2:21])

# Part (c): premium APV factor and premium
prem_factor <- 1 + sum(v^(1:19) * out$tp00[2:20])
P <- (death_apv + disab_apv) / prem_factor

ans <- c(
  death_apv = death_apv,
  disab_apv = disab_apv,
  premium_factor = prem_factor,
  P = P
)
print(ans)

# Plot state probabilities
plot(
  out$t, out$tp00, type = "l", lwd = 2,
  ylim = c(0, 1),
  xlab = "t", ylab = "probability",
  main = "Example 14.10: Euler state probabilities"
)
lines(out$t, out$tp01, lwd = 2, lty = 2)
lines(out$t, out$tp02, lwd = 2, lty = 3)
legend(
  "right",
  legend = c(expression(tp[x]^{"00"}),
             expression(tp[x]^{"01"}),
             expression(tp[x]^{"02"})),
  lwd = 2, lty = c(1, 2, 3), bty = "n"
)

# Plot discounted contribution streams
k <- 1:20
death_contrib <- 10000 * v^k * out$tp00[1:20] * q_death_step
disab_contrib <- 1000 * v^k * out$tp01[2:21]

plot(
  k / 10, death_contrib, type = "h", lwd = 2,
  xlab = "time", ylab = "discounted contribution",
  main = "Example 14.10: APV contribution streams"
)
lines(k / 10, disab_contrib, type = "h", lwd = 2, lty = 2)
legend(
  "topright",
  legend = c("death benefit", "disability benefit"),
  lwd = 2, lty = c(1, 2), bty = "n"
)


# Chapter 14 R Check: Example 14.13
# Constant-force permanent disability model

delta <- 0.05
mu01 <- 0.005
mu02 <- 0.02
mu12 <- 0.10

t <- seq(0, 200, by = 0.01)

# State-0 survival
tp00 <- exp(-(mu01 + mu02) * t)

# (a) Probability of ever becoming disabled
prob_disab <- sum(diff(t) * (head(tp00 * mu01, -1) + tail(tp00 * mu01, -1)) / 2)

# (b) Expected time spent disabled
# = P(enter state 1) * E(length in state 1 | entered)
exp_time_disabled <- prob_disab * (1 / mu12)

# (c) APV
integrand_death <- 100000 * exp(-delta * t) * tp00 * mu02
integrand_disab <- 500000 * exp(-delta * t) * tp00 * mu01

apv_death <- sum(diff(t) * (head(integrand_death, -1) + tail(integrand_death, -1)) / 2)
apv_disab <- sum(diff(t) * (head(integrand_disab, -1) + tail(integrand_disab, -1)) / 2)
apv_total <- apv_death + apv_disab

out <- c(
  prob_disab = prob_disab,
  exp_time_disabled = exp_time_disabled,
  apv_death = apv_death,
  apv_disab = apv_disab,
  apv_total = apv_total
)
print(out)

# Plot APV integrands
plot(
  t, integrand_death, type = "l", lwd = 2,
  xlab = "t", ylab = "discounted contribution density",
  main = "Example 14.13: APV integrands"
)
lines(t, integrand_disab, lwd = 2, lty = 2)
legend(
  "topright",
  legend = c("death while healthy", "disability payment"),
  lwd = 2, lty = c(1, 2), bty = "n"
)




# =========================================================
# Chapter 14 R Check: Exercise 14-8
# Cumulative withdrawal probability from Example 14.3 data
# =========================================================

library(mqriskR)

q1 <- c(.02, .03, .04, .05, .06)
q2 <- c(.30, .20, .20, .10, .00)
p_tau <- c(.68, .77, .76, .85, .94)

# cumulative withdrawal probability through duration n
cum_withdrawal <- cumsum(c(
  q2[1],
  p_tau[1] * q2[2],
  p_tau[1] * p_tau[2] * q2[3],
  p_tau[1] * p_tau[2] * p_tau[3] * q2[4],
  p_tau[1] * p_tau[2] * p_tau[3] * p_tau[4] * q2[5]
))

out <- data.frame(
  n = 1:5,
  cum_withdrawal = cum_withdrawal
)

print(out)

# value requested in part (b)
q3_withdrawal <- cum_withdrawal[3]
print(q3_withdrawal)

plot(
  out$n, out$cum_withdrawal, type = "b", pch = 19,
  xlab = "Duration n",
  ylab = "Cumulative withdrawal probability",
  main = "Exercise 14-8: cumulative withdrawal probability"
)
abline(h = q3_withdrawal, lty = 2)



# =========================================================
# Chapter 14 R Check: Exercise 14-14
# Euler approximation for employment / unemployment model
# =========================================================

library(mqriskR)

mu01 <- function(t) 0.20 + 0.0002 * t^2   # employed -> unemployed
mu02 <- function(t) 0.05                  # employed -> deceased
mu10 <- function(t) 0.80 - 0.04 * t       # unemployed -> employed
mu12 <- function(t) 0.05                  # unemployed -> deceased

out_h1 <- tp00_tp01_euler(
  h = 1, n = 20,
  mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
)

out_h05 <- tp00_tp01_euler(
  h = 0.5, n = 10,
  mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
)

out_h01 <- tp00_tp01_euler(
  h = 0.1, n = 10,
  mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
)

# comparison at t = 10
comp <- data.frame(
  step = c(1, 0.5, 0.1),
  tp00_t10 = c(
    subset(out_h1, abs(t - 10) < 1e-12)$tp00,
    subset(out_h05, abs(t - 10) < 1e-12)$tp00,
    subset(out_h01, abs(t - 10) < 1e-12)$tp00
  ),
  tp01_t10 = c(
    subset(out_h1, abs(t - 10) < 1e-12)$tp01,
    subset(out_h05, abs(t - 10) < 1e-12)$tp01,
    subset(out_h01, abs(t - 10) < 1e-12)$tp01
  )
)

print(comp)

# plot for h = 0.1
plot(
  out_h01$t, out_h01$tp00, type = "l", lwd = 2,
  ylim = c(0, 1),
  xlab = "t", ylab = "Probability",
  main = "Exercise 14-14: Euler approximations (h = 0.1)"
)
lines(out_h01$t, out_h01$tp01, lwd = 2, lty = 2)
lines(out_h01$t, out_h01$tp02, lwd = 2, lty = 3)
legend(
  "topright",
  legend = c(expression(tp[x]^{"00"}),
             expression(tp[x]^{"01"}),
             expression(tp[x]^{"02"})),
  lwd = 2, lty = c(1, 2, 3), bty = "n"
)


# =========================================================
# Chapter 14 R Check: Exercise 14-15
# APV of monthly unemployment benefit
# =========================================================

library(mqriskR)

mu01 <- function(t) 0.20 + 0.0002 * t^2
mu02 <- function(t) 0.05
mu10 <- function(t) 0.80 - 0.04 * t
mu12 <- function(t) 0.05

h <- 1 / 12
n <- 20

out <- tp00_tp01_euler(
  h = h, n = n,
  mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
)

j <- 0.04 / 12
k <- 1:240

apv_unemp <- 1400 * sum((1 + j)^(-k) * out$tp01[2:(240 + 1)])
print(apv_unemp)

contrib <- 1400 * (1 + j)^(-k) * out$tp01[2:(240 + 1)]

plot(
  k / 12, contrib, type = "l", lwd = 2,
  xlab = "Time (years)",
  ylab = "Discounted monthly contribution",
  main = "Exercise 14-15: APV contribution pattern"
)
abline(h = 0, lty = 3)



# =========================================================
# Chapter 14 R Check: Exercises 14-19 and 14-20
# Discrete-time Markov chain probabilities for the CCRC model
# =========================================================

library(mqriskR)

P <- matrix(
  c(0.94, 0.03, 0.02, 0.01,
    0.50, 0.30, 0.18, 0.02,
    0.00, 0.00, 0.93, 0.07,
    0.00, 0.00, 0.00, 1.00),
  nrow = 4, byrow = TRUE
)

# Exercise 14-19: P(X3 = 0 | X0 = 0)
p00_3 <- markov_nstep_prob(P, n = 3, i = 1, j = 1)

# Exercise 14-20: P(X3 = 2 | X0 = 0)
p02_3 <- markov_nstep_prob(P, n = 3, i = 1, j = 3)

print(c(p00_3 = p00_3, p02_3 = p02_3))

# plot selected n-step probabilities from state 0
n_vals <- 0:12
p00 <- sapply(n_vals, function(n) markov_nstep_prob(P, n = n, i = 1, j = 1))
p01 <- sapply(n_vals, function(n) markov_nstep_prob(P, n = n, i = 1, j = 2))
p02 <- sapply(n_vals, function(n) markov_nstep_prob(P, n = n, i = 1, j = 3))
p03 <- sapply(n_vals, function(n) markov_nstep_prob(P, n = n, i = 1, j = 4))

plot(
  n_vals, p00, type = "l", lwd = 2,
  ylim = c(0, 1),
  xlab = "Number of months",
  ylab = "Probability",
  main = "Exercises 14-19 and 14-20: n-step state probabilities"
)
lines(n_vals, p01, lwd = 2, lty = 2)
lines(n_vals, p02, lwd = 2, lty = 3)
lines(n_vals, p03, lwd = 2, lty = 4)
legend(
  "right",
  legend = c("State 0", "State 1", "State 2", "State 3"),
  lwd = 2, lty = c(1, 2, 3, 4), bty = "n"
)


# =========================================================
# Chapter 14 R Check: Exercise 14-23
# Risk-class transition matrix, projections, and average claim costs
# =========================================================

library(mqriskR)

P <- matrix(
  c(.826, .168, .006,
    .345, .617, .038,
    .111, .500, .389),
  nrow = 3, byrow = TRUE
)

pi1 <- c(.675, .303, .022)
claim_costs_Z <- c(500, 10000, 50000)

pi2 <- as.numeric(pi1 %*% P)
pi3 <- as.numeric(pi2 %*% P)

avg_Z  <- sum(c(.695, .287, .018) * claim_costs_Z)
avg_Z3 <- (1.05)^3 * sum(pi3 * claim_costs_Z)

print(P)
print(rbind(pi1 = pi1, pi2 = pi2, pi3 = pi3))
print(c(avg_Z = avg_Z, avg_Z3 = avg_Z3))

# plot projected class proportions
year <- 1:3
matplot(
  year,
  rbind(pi1, pi2, pi3),
  type = "l", lwd = 2, lty = 1:3,
  xlab = "Year index (1 = Z+1)",
  ylab = "Proportion",
  main = "Exercise 14-23: projected risk-class proportions"
)
legend(
  "topright",
  legend = c("Low", "Moderate", "High"),
  lwd = 2, lty = 1:3, bty = "n"
)

# plot average claim cost path
avg_Z1 <- sum(pi1 * claim_costs_Z) * 1.05
avg_Z2 <- sum(pi2 * claim_costs_Z) * (1.05)^2
avg_Z3_path <- sum(pi3 * claim_costs_Z) * (1.05)^3

avg_path <- c(avg_Z, avg_Z1, avg_Z2, avg_Z3_path)
plot(
  0:3, avg_path, type = "b", pch = 19,
  xlab = "Years after Z",
  ylab = "Average claim cost",
  main = "Exercise 14-23: average claim cost path"
)



# =========================================================
# Chapter 14 R Check: Exercise 14-25
# Mortality gain and withdrawal gain
# =========================================================

library(mqriskR)

# assumed values
V9 <- 115.00
V10 <- 128.83
G <- 16
e <- 3
i <- 0.06
DB <- 1000
WB <- 110

# assumed and actual decrement rates
q_d_assumed <- 0.01
q_w_assumed <- 0.10
q_d_actual <- 15 / 1000
q_w_actual <- 100 / 985

# gain from mortality: actual mortality, assumed withdrawal
GM <- gain_loss_md(
  Vt = V9,
  G = G,
  r = 0,
  e = e,
  i = i,
  b1 = DB,
  b2 = WB,
  s1 = 0,
  s2 = 0,
  q1 = q_d_assumed,
  q2 = q_w_assumed,
  Vt1 = V10,
  year_end_cause2 = TRUE,
  q1prime = q_d_actual,
  q2prime = q_w_assumed
)

# gain from mortality + withdrawal: actual mortality, actual withdrawal
GMW <- gain_loss_md(
  Vt = V9,
  G = G,
  r = 0,
  e = e,
  i = i,
  b1 = DB,
  b2 = WB,
  s1 = 0,
  s2 = 0,
  q1 = q_d_assumed,
  q2 = q_w_assumed,
  Vt1 = V10,
  year_end_cause2 = TRUE,
  q1prime = q_d_actual,
  q2prime = q_w_actual
)

GW <- GMW - GM

out <- c(
  gain_mortality = GM,
  gain_mortality_and_withdrawal = GMW,
  gain_withdrawal = GW
)
print(out)

# simple comparison plot
assumed_actual <- rbind(
  assumed = c(q_d_assumed, q_w_assumed),
  actual  = c(q_d_actual, q_w_actual)
)
colnames(assumed_actual) <- c("Mortality", "Withdrawal")

barplot(
  t(assumed_actual), beside = TRUE,
  main = "Exercise 14-25: assumed vs actual decrement rates",
  ylab = "Rate",
  legend.text = TRUE,
  args.legend = list(bty = "n", x = "topright")
)

