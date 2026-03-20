# R Check for Example 11.1
# Verify: beta = P + (beta - alpha)/adotxn

library(mqriskR)

x <- 40
k <- 20
i <- 0.05

P <- Pxn(x = x, n = k, i = i, model = "uniform", omega = 100)
adot_temp <- adotxn(x = x, n = k, i = i, model = "uniform", omega = 100)

gap <- 0.02
beta_from_formula <- P + gap / adot_temp
alpha_from_gap <- beta_from_formula - gap

print(round(c(
  P = P,
  adot_xn = adot_temp,
  gap = gap,
  beta = beta_from_formula,
  alpha = alpha_from_gap
), 6))

# Sensitivity plot
gap_grid <- seq(0, 0.10, by = 0.002)
beta_grid <- P + gap_grid / adot_temp

plot(
  gap_grid, beta_grid, type = "l", lwd = 2,
  xlab = expression(beta - alpha),
  ylab = expression(beta),
  main = "Example 11.1: Renewal Modified Premium vs (beta - alpha)"
)
abline(h = P, lty = 2)




# R Check for Example 11.4
# Fully continuous whole life with increasing benefit b_t = 1000 e^(0.04 t)

lambda <- 0.02
delta <- 0.04
t <- 2

# Net premium rate
abar <- 1 / (lambda + delta)
Pbar <- 1000 / abar

# APV of future benefit at duration t
apvb_t <- integrate(
  function(s) 1000 * lambda * exp(0.04 * t) * exp(-lambda * s),
  lower = 0, upper = Inf
)$value

apvb_t

# APV of future premiums at duration t
apvp_t <- Pbar * abar

Vt <- apvb_t - apvp_t

print(round(c(
  abar = abar,
  Pbar = Pbar,
  apvb_t = apvb_t,
  apvp_t = apvp_t,
  reserve_t = Vt
), 6))

# Reserve as a function of duration
t_grid <- seq(0, 10, by = 0.25)
reserve_grid <- 1000 * exp(0.04 * t_grid) - Pbar * abar

plot(
  t_grid, reserve_grid, type = "l", lwd = 2,
  xlab = "Duration t",
  ylab = expression(t * bar(V)),
  main = "Example 11.4: Continuous Reserve by Duration"
)
abline(v = t, lty = 2)
abline(h = Vt, lty = 2)


# R Check for Example 11.7
# Gross gain on a single policy

V10G <- 3950.73
V11G <- 4607.07
G <- 685.00

i_actual <- 0.065
q_actual <- 0.005
r_actual <- 0.06
s_actual <- 100
b <- 50000

GT <- (V10G + G * (1 - r_actual)) * (1 + i_actual) -
  ((b + s_actual) * q_actual + (1 - q_actual) * V11G)

print(round(c(
  total_gain = GT
), 6))

# Sensitivity to actual interest
i_grid <- seq(0.04, 0.08, by = 0.001)
GT_i <- sapply(i_grid, function(ii) {
  (V10G + G * (1 - r_actual)) * (1 + ii) -
    ((b + s_actual) * q_actual + (1 - q_actual) * V11G)
})

plot(
  i_grid, GT_i, type = "l", lwd = 2,
  xlab = "Actual interest rate",
  ylab = "Total gain",
  main = "Example 11.7: Total Gain vs Actual Interest Rate"
)
abline(v = i_actual, lty = 2)

# Sensitivity to actual mortality
q_grid <- seq(0.003, 0.008, by = 0.0001)
GT_q <- sapply(q_grid, function(qq) {
  (V10G + G * (1 - r_actual)) * (1 + i_actual) -
    ((b + s_actual) * qq + (1 - qq) * V11G)
})

plot(
  q_grid, GT_q, type = "l", lwd = 2,
  xlab = "Actual mortality rate",
  ylab = "Total gain",
  main = "Example 11.7: Total Gain vs Actual Mortality Rate"
)
abline(v = q_actual, lty = 2)



# R Check for Example 11.9
# Backward Euler approximation for gross reserve

dVdt <- function(t, V) {
  19.60 + 0.04 * V - 0.01 * t^2 * (110 - V)
}

# Backward recursion using derivative at left endpoint
backward_left <- function(times, V_terminal) {
  out <- numeric(length(times))
  out[length(times)] <- V_terminal

  for (k in seq(length(times) - 1, 1, by = -1)) {
    h <- times[k + 1] - times[k]
    t_left <- times[k]
    # solve: V_k + h*dVdt(t_left, V_k) = V_{k+1}
    a <- 1 + h * (0.04 + 0.01 * t_left^2)
    rhs <- out[k + 1] - h * 19.60 + h * 110 * 0.01 * t_left^2
    out[k] <- rhs / a
  }
  out
}

# Backward recursion using derivative at right endpoint
backward_right <- function(times, V_terminal) {
  out <- numeric(length(times))
  out[length(times)] <- V_terminal

  for (k in seq(length(times) - 1, 1, by = -1)) {
    h <- times[k + 1] - times[k]
    t_right <- times[k + 1]
    out[k] <- out[k + 1] - h * dVdt(t_right, out[k + 1])
  }
  out
}

times1 <- seq(4, 5, by = 0.5)
left_path <- backward_left(times1, 100)
right_path <- backward_right(times1, 100)

print(round(data.frame(
  time = times1,
  backward_left = left_path,
  backward_right = right_path
), 6))

# Smaller step size
times2 <- seq(4, 5, by = 0.25)
left_path_small <- backward_left(times2, 100)

plot(
  times1, left_path, type = "b", lwd = 2, pch = 16,
  ylim = range(c(left_path, right_path, left_path_small)),
  xlab = "Time",
  ylab = expression(t * bar(V)^G),
  main = "Example 11.9: Backward Euler Reserve Approximations"
)
lines(times1, right_path, type = "b", lwd = 2, pch = 1, lty = 2)
lines(times2, left_path_small, type = "b", lwd = 2, pch = 3, lty = 3)
legend(
  "topleft",
  legend = c("h = 0.50, left derivative", "h = 0.50, right derivative", "h = 0.25, left derivative"),
  lty = c(1, 2, 3),
  pch = c(16, 1, 3),
  bty = "n"
)




# R Check for Exercise 11-1
# Verify:
# (a) V_NLP - V_M = (beta - P) * adot_{x+t:n-t}
# (b) V_M = 1 - (beta + d) * adot_{x+t:n-t}

library(mqriskR)

x <- 40
n <- 20
t <- 10
i <- 0.05
beta <- 0.08

P <- Pxn(x = x, n = n, i = i, model = "uniform", omega = 100)
d <- i / (1 + i)
adot_future <- adotxn(x = x + t, n = n - t, i = i, model = "uniform", omega = 100)
A_future <- Axn(x = x + t, n = n - t, i = i, model = "uniform", omega = 100)

V_nlp <- tVxn(x = x, n = n, t = t, i = i, model = "uniform", omega = 100)
V_mod <- A_future - beta * adot_future

lhs_a <- V_nlp - V_mod
rhs_a <- (beta - P) * adot_future

lhs_b <- V_mod
rhs_b <- 1 - (beta + d) * adot_future

round(c(
  P = P,
  V_nlp = V_nlp,
  V_mod = V_mod,
  lhs_a = lhs_a,
  rhs_a = rhs_a,
  lhs_b = lhs_b,
  rhs_b = rhs_b
), 6)

# Sensitivity plot
beta_grid <- seq(P - 0.01, P + 0.06, by = 0.001)
V_mod_grid <- A_future - beta_grid * adot_future

plot(
  beta_grid, V_mod_grid, type = "l", lwd = 2,
  xlab = expression(beta),
  ylab = "Modified reserve",
  main = "Exercise 11-1: Modified Reserve vs Renewal Premium"
)
abline(v = P, lty = 2)
abline(h = V_nlp, lty = 3)




# R Check for Exercise 11-5
# Verify r = (s|1-s q_x+t) / q_x+t = 1 - s under UDD

s <- seq(0, 1, by = 0.01)

# Under UDD:
r_formula <- 1 - s

# Numerical illustration using a sample q value
q <- 0.12
deferred_piece <- (1 - s) * q
r_numeric <- deferred_piece / q

round(head(data.frame(
  s = s,
  r_numeric = r_numeric,
  r_formula = r_formula
), 10), 6)

plot(
  s, r_numeric, type = "l", lwd = 2,
  xlab = "s",
  ylab = "r",
  main = "Exercise 11-5: Interpolation Weight under UDD"
)
lines(s, r_formula, lty = 2, lwd = 2)
legend(
  "topright",
  legend = c("Numerical ratio", "1 - s"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)




# R Check for Exercise 11-6
# Verify:
# V_{t-1/2} = 1/2 * [ (1 + v p_{x+t-1}) V_t + v q_{x+t-1} ]

library(mqriskR)

x <- 40
t <- 10
i <- 0.05
v <- 1 / (1 + i)

Vt <- tVx(x = x, t = t, i = i, model = "uniform", omega = 100)
px <- tpx(x = x + t - 1, t = 1, model = "uniform", omega = 100)
qx <- 1 - px

# meanVx(x, t, ...) returns {}_{t+1/2}V, so use t - 1 to match {}_{t-1/2}V
lhs <- meanVx(x = x, t = t - 1, i = i, model = "uniform", omega = 100)
rhs <- 0.5 * ((1 + v * px) * Vt + v * qx)

round(c(
  Vt = Vt,
  px = px,
  qx = qx,
  lhs = lhs,
  rhs = rhs
), 6)

# Check across a range of durations
t_grid <- 1:20
lhs_grid <- sapply(t_grid, function(tt) {
  meanVx(x = x, t = tt - 1, i = i, model = "uniform", omega = 100)
})
rhs_grid <- sapply(t_grid, function(tt) {
  Vtt <- tVx(x = x, t = tt, i = i, model = "uniform", omega = 100)
  ptt <- tpx(x = x + tt - 1, t = 1, model = "uniform", omega = 100)
  qtt <- 1 - ptt
  0.5 * ((1 + v * ptt) * Vtt + v * qtt)
})

plot(
  t_grid, lhs_grid, type = "l", lwd = 2,
  xlab = "t",
  ylab = "Mean reserve",
  main = "Exercise 11-6: Mean Reserve Identity"
)
lines(t_grid, rhs_grid, lty = 2, lwd = 2)
legend(
  "topleft",
  legend = c("Package value for {}_{t-1/2}V", "Formula value"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)




# R Check for Exercise 11-8
# 3-year term insurance with decreasing benefits

i <- 0.06
v <- 1 / (1 + i)

b1 <- 200
b2 <- 150
b3 <- 100

qx <- 0.03
qx1 <- 0.06
qx2 <- 0.09

px <- 1 - qx
px1 <- 1 - qx1

P <- (
  b1 * v * qx +
    b2 * v^2 * px * qx1 +
    b3 * v^3 * px * px1 * qx2
) / (
  1 + v * px + v^2 * px * px1
)

V1 <- ((P) * (1 + i) - b1 * qx) / px
initial_second_year <- V1 + P

round(c(
  premium = P,
  terminal_reserve_year1 = V1,
  initial_reserve_year2 = initial_second_year
), 6)

# Sensitivity plot for b2
b2_grid <- seq(100, 250, by = 2)
init2_grid <- sapply(b2_grid, function(bb2) {
  Ptmp <- (
    b1 * v * qx +
      bb2 * v^2 * px * qx1 +
      b3 * v^3 * px * px1 * qx2
  ) / (
    1 + v * px + v^2 * px * px1
  )
  V1tmp <- (Ptmp * (1 + i) - b1 * qx) / px
  V1tmp + Ptmp
})

plot(
  b2_grid, init2_grid, type = "l", lwd = 2,
  xlab = expression(b[2]),
  ylab = "Second-year initial reserve",
  main = "Exercise 11-8: Initial Reserve vs Middle-Year Benefit"
)
abline(v = b2, lty = 2)




# R Check for Exercise 11-9
# 2-year endowment with death benefit = 1000 + reserve at end of year of death

i <- 0.10
q0 <- 0.10
q1 <- 0.11

# From the solution:
# V1 = P(1+i) - 100
# V2 = (V1 + P)(1+i) - 110
# Set V2 = 1000 and solve for P

P <- (1000 + 220) / (((1 + i) + 1) * (1 + i))

V1 <- P * (1 + i) - 100
V2 <- (V1 + P) * (1 + i) - 110

round(c(
  premium = P,
  reserve_1 = V1,
  reserve_2 = V2
), 6)

# Sensitivity to first-year mortality
q0_grid <- seq(0.02, 0.20, by = 0.002)
P_grid <- sapply(q0_grid, function(qq0) {
  # Generalized version:
  # V1 = P(1+i) - 1000 q0
  # V2 = (V1 + P)(1+i) - 1000 q1
  # Set V2 = 1000
  (1000 + 1000 * (qq0 * (1 + i) + q1)) / (((1 + i) + 1) * (1 + i))
})

plot(
  q0_grid, P_grid, type = "l", lwd = 2,
  xlab = expression(q[x]),
  ylab = "Net level premium",
  main = "Exercise 11-9: Premium vs First-Year Mortality"
)
abline(v = q0, lty = 2)




# R Check for Exercise 11-10
# Solve for b_11 from the recursion:
# (V_10 + P)(1+i) = b_11 q_30' + p_30' V_11
#
# To avoid inconsistency, compute P, V10, and V11 from the same model basis.

library(mqriskR)

i <- 0.05

# Question data
q30 <- 0.008427
q30_prime <- q30 + 0.01
p30_prime <- 1 - q30_prime

# Since premiums and reserves are stated to match those of issue age 20,
# compute all model-based quantities from the same age-20 basis.
P <- Px(x = 20, i = i, model = "uniform", omega = 100)
V10 <- tVx(x = 20, t = 10, i = i, model = "uniform", omega = 100)
V11 <- tVx(x = 20, t = 11, i = i, model = "uniform", omega = 100)

# Implied benefit in year 11
b11 <- ((V10 + P) * (1 + i) - p30_prime * V11) / q30_prime

round(c(
  P = P,
  V10 = V10,
  V11 = V11,
  q30 = q30,
  q30_prime = q30_prime,
  b11 = b11
), 6)

# Check the recursion numerically
lhs <- (V10 + P) * (1 + i)
rhs <- b11 * q30_prime + p30_prime * V11

round(c(
  lhs = lhs,
  rhs = rhs,
  difference = lhs - rhs
), 10)

# Sensitivity to excess mortality adjustment
excess_grid <- seq(0, 0.03, by = 0.0005)
b_grid <- sapply(excess_grid, function(extra_q) {
  q_prime <- q30 + extra_q
  p_prime <- 1 - q_prime
  ((V10 + P) * (1 + i) - p_prime * V11) / q_prime
})

plot(
  excess_grid, b_grid, type = "l", lwd = 2,
  xlab = "Excess mortality adjustment",
  ylab = expression(b[11]),
  main = "Exercise 11-10: Implied Benefit vs Excess Mortality"
)
abline(v = 0.01, lty = 2)
abline(h = b11, lty = 2)




# R Check for Exercise 11-19
# Ordered gain decomposition for 990 policies:
# (a) interest, (b) mortality, (c) expense

npol <- 990
G <- 90
V3G <- 97.12
V4G <- 164.13
benefit <- 10000

# Assumed
i_assumed <- 0.05
q_assumed <- 0.003
r_assumed <- 0.03

# Actual
i_actual <- 0.04
q_actual <- 0.002
r_actual <- 0.025

# Per-policy total gain under actual experience
GT_per <- (V3G + G * (1 - r_actual)) * (1 + i_actual) -
  (benefit * q_actual + (1 - q_actual) * V4G)

# Interest gain first
GI_per <- (V3G + G * (1 - r_assumed)) * (1 + i_actual) -
  (benefit * q_assumed + (1 - q_assumed) * V4G)

# Mortality gain second
GIM_per <- (V3G + G * (1 - r_assumed)) * (1 + i_actual) -
  (benefit * q_actual + (1 - q_actual) * V4G)
GM_per <- GIM_per - GI_per

# Expense gain last
GE_per <- GT_per - GI_per - GM_per

out <- c(
  total_gain = npol * GT_per,
  interest_gain = npol * GI_per,
  mortality_gain = npol * GM_per,
  expense_gain = npol * GE_per,
  check = npol * (GI_per + GM_per + GE_per)
)

round(out, 6)


## Use package function
out_per_policy <- decompGg_disc(
  VtG = 97.12,
  Vt1G = 164.13,
  G = 90,
  i_assumed = 0.05,
  q_assumed = 0.003,
  r_assumed = 0.03,
  e_assumed = 0,
  s_assumed = 0,
  i_actual = 0.04,
  q_actual = 0.002,
  r_actual = 0.025,
  e_actual = 0,
  s_actual = 0,
  b = 10000,
  order = c("interest", "mortality", "expense")
)

out_block <- npol * out_per_policy

round(out_block, 6)

barplot(
  height = c(out["interest_gain"], out["mortality_gain"], out["expense_gain"]),
  names.arg = c("Interest", "Mortality", "Expense"),
  main = "Exercise 11-19: Aggregate Gain by Source",
  ylab = "Gain for 990 policies"
)
abline(h = 0)





# R Check for Exercise 11-20
# Payout-phase deferred annuity with expense on survival benefit

library(mqriskR)

x <- 70
i_assumed <- 0.06
i_actual <- 0.055
q_assumed <- 0.02041
q_actual <- 0.025
p_assumed <- 1 - q_assumed
p_actual <- 1 - q_actual

benefit <- 10000
expense_rate <- 0.05
net_survival_outgo <- benefit * (1 + expense_rate)   # annuity + expense if alive

# Reserve at age 70 and 71 from annual life annuity-immediate values
# payout phase: reserve is APV of future survival-contingent payments
V70 <- net_survival_outgo * ax(x = x, i = i_assumed, model = "uniform", omega = 100)
V71 <- net_survival_outgo * ax(x = x + 1, i = i_assumed, model = "uniform", omega = 100)

# Ordered gains: mortality first, then interest
GM <- (V70 - net_survival_outgo) * (1 + i_assumed) - p_actual * V71
GI_total <- (V70 - net_survival_outgo) * (1 + i_actual) - p_actual * V71
GI <- GI_total - GM

round(c(
  V70 = V70,
  V71 = V71,
  mortality_gain = GM,
  interest_gain = GI,
  total_gain = GM + GI
), 6)

# Sensitivity to actual mortality
q_grid <- seq(0.01, 0.04, by = 0.0005)
GM_grid <- sapply(q_grid, function(qa) {
  pa <- 1 - qa
  (V70 - net_survival_outgo) * (1 + i_assumed) - pa * V71
})

plot(
  q_grid, GM_grid, type = "l", lwd = 2,
  xlab = "Actual mortality rate",
  ylab = "Mortality gain",
  main = "Exercise 11-20: Mortality Gain in Payout Phase"
)
abline(v = q_actual, lty = 2)

