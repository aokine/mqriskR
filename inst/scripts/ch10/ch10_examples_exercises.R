# R Check for Example 10.1
# Verify {}_{n-1}V_{x:\angl{n}} both ways:
# (1) prospective reserve formula
# (2) shortcut v - P_{x:\angl{n}}

A_xn <- 0.20
d <- 0.08
i <- d / (1 - d)
v <- 1 / (1 + i)

# Premium for the n-year endowment
adot_xn <- (1 - A_xn) / d
P_xn <- A_xn / adot_xn

# Reserve from the prospective formula at t = n - 1
# For one year remaining in an endowment:
# A_{x+n-1:\angl{1}} = v and ä_{x+n-1:\angl{1}} = 1
reserve_prospective <- v - P_xn * 1

# Reserve from the shortcut expression
reserve_shortcut <- v - P_xn

out <- c(
  i = i,
  v = v,
  premium = P_xn,
  reserve_prospective = reserve_prospective,
  reserve_shortcut = reserve_shortcut,
  difference = reserve_prospective - reserve_shortcut
)

print(round(out, 6))


# R Check for Example 10.2
# Solve for P_{x:\angl{n}}^1 from
# {}_nV_x = P_x / P_{x:\angl{n}}^{1(endowment)} - P_{x:\angl{n}}^1 / P_{x:\angl{n}}^{1(endowment)}

Px_whole <- 0.090
Vn <- 0.563
PnEx <- 0.00864   # this is P_{x:\angl{n}}^{ 1 } for the pure endowment part

Pxn1 <- Px_whole - Vn * PnEx

print(round(c(
  Px_whole = Px_whole,
  Vn = Vn,
  PnEx = PnEx,
  Pxn1 = Pxn1
), 6))

# Sensitivity plot: implied term premium as reserve changes
V_grid <- seq(0, 1.2, by = 0.01)
Pxn1_grid <- Px_whole - V_grid * PnEx

plot(
  V_grid, Pxn1_grid, type = "l", lwd = 2,
  xlab = "Reserve value",
  ylab = "Implied term premium",
  main = "Example 10.2: Implied Term Premium vs Reserve"
)
abline(v = Vn, lty = 2)
abline(h = Pxn1, lty = 2)


# R Check for Example 10.5
# Illustrate that {}_tV = P * s̈_t for t <= n
# and P = ä_{x+n} / s̈_n

library(mqriskR)

x <- 55
n <- 10
i <- 0.05
omega <- 100

d <- i / (1 + i)

adot_at_defer_end <- adotx(
  x = x + n,
  i = i,
  model = "uniform",
  omega = omega
)

sdot_n <- ((1 + i)^n - 1) / d
P <- adot_at_defer_end / sdot_n

t_vals <- 0:n
sdot_t <- ((1 + i)^t_vals - 1) / d
reserve_vals <- P * sdot_t

out <- data.frame(
  t = t_vals,
  sdot_t = sdot_t,
  reserve = reserve_vals
)

print(round(out, 6))
print(round(c(
  premium = P,
  reserve_at_n = tail(reserve_vals, 1),
  target_annuity_value = adot_at_defer_end
), 6))

plot(
  t_vals, reserve_vals, type = "b", pch = 19,
  xlab = "Duration t",
  ylab = "Reserve",
  main = "Example 10.5: Reserve Path During Deferred Period"
)
abline(h = adot_at_defer_end, lty = 2)



# R Check for Example 10.6
# q' - q = (.01)(1.03) / (1 - {}_{20}V_x)

V20 <- 0.427
premium_extra <- 0.01
i <- 0.03

q_excess <- premium_extra * (1 + i) / (1 - V20)

print(round(c(
  V20 = V20,
  premium_extra = premium_extra,
  i = i,
  q_excess = q_excess
), 6))

# Sensitivity to reserve level
V_grid <- seq(0.05, 0.90, by = 0.01)
q_excess_grid <- premium_extra * (1 + i) / (1 - V_grid)

plot(
  V_grid, q_excess_grid, type = "l", lwd = 2,
  xlab = "Reserve at duration 20",
  ylab = expression(q[x+19] * "'" - q[x+19]),
  main = "Example 10.6: Implied Excess Mortality vs Reserve"
)
abline(v = V20, lty = 2)
abline(h = q_excess, lty = 2)




# R Check for Example 10.8
# Verify {}_tV_x^{(m)} = [1 + beta(m) P_x^{(m)}] {}_tV_x under UDD

library(mqriskR)

x <- 40
t <- 10
i <- 0.05
m_vals <- c(2, 4, 12)
omega <- 100

annual_reserve <- tVx(
  x = x,
  t = t,
  i = i,
  model = "uniform",
  omega = omega
)

check_df <- lapply(m_vals, function(m) {
  ic <- interest_convert(i = i, m = m)
  im <- ic$im
  dm <- (im / m) / (1 + im / m)
  beta <- (i - im) / (im * dm)

  prem_m <- Px_m(
    x = x,
    m = m,
    i = i,
    model = "uniform",
    omega = omega
  )

  lhs <- tVx_m(
    x = x,
    t = t,
    m = m,
    i = i,
    model = "uniform",
    omega = omega
  )

  rhs <- (1 + beta * prem_m) * annual_reserve

  data.frame(
    m = m,
    annual_reserve = annual_reserve,
    premium_m = prem_m,
    lhs = lhs,
    rhs = rhs,
    difference = lhs - rhs
  )
})

check_df <- do.call(rbind, check_df)
print(round(check_df, 10))

plot(
  check_df$m, check_df$difference, type = "b", pch = 19,
  xlab = "m",
  ylab = "LHS - RHS",
  main = "Example 10.8: UDD Reserve Identity Check"
)
abline(h = 0, lty = 2)




# R Check for Exercise 10-3
# Compute {}_{20}V_{45} from the given premium relationships.

P45 <- 0.014
P45_20 <- 0.030
P45_20_term <- 0.022

reserve_20 <- P45 / P45_20_term - (P45_20 - P45_20_term) / P45_20_term

print(round(c(
  P45 = P45,
  P45_20 = P45_20,
  P45_20_term = P45_20_term,
  reserve_20 = reserve_20
), 6))

# Sensitivity plot: reserve as P45 varies slightly
P45_grid <- seq(0.010, 0.018, by = 0.0005)
reserve_grid <- P45_grid / P45_20_term - (P45_20 - P45_20_term) / P45_20_term

plot(
  P45_grid, reserve_grid, type = "l", lwd = 2,
  xlab = "Whole life premium P_45",
  ylab = "Reserve at duration 20",
  main = "Exercise 10-3: Retrospective Reserve Sensitivity"
)
abline(v = P45, lty = 2)
abline(h = reserve_20, lty = 2)





# R Check for Exercise 10-6
# Compute 1000({}_2V_{x:\angl{3}} - {}_1V_{x:\angl{3}})

P_x3 <- 0.33251
i <- 0.06
v <- 1 / (1 + i)

lx <- 100
lx1 <- 90
qx <- 1 - lx1 / lx

# Reserve at t = 1, retrospective
V1 <- P_x3 / ((1 / (1 + i)) * (1 - qx)) - (v * qx) / ((1 / (1 + i)) * (1 - qx))

# Reserve at t = 2, prospective
V2 <- v - P_x3

diff_val <- 1000 * (V2 - V1)

print(round(c(
  qx = qx,
  V1 = V1,
  V2 = V2,
  difference_1000 = diff_val
), 6))

# Bar plot
barplot(
  c(1000 * V1, 1000 * V2),
  names.arg = c("1000 V1", "1000 V2"),
  main = "Exercise 10-6: Reserves at Durations 1 and 2"
)
abline(h = c(1000 * V1, 1000 * V2), lty = 2)



# R Check for Exercise 10-9
# Find Pr(L < 190) for a 2-year term insurance of amount 400.

P <- 0.185825
V1 <- 0.04145
i <- 0.10
v <- 1 / (1 + i)

# From {}_1V_{x:\angl{2}}^1 = v q_{x+1} - P
q_x1 <- (V1 + P) / v

# Solve for p_x from premium equation
f <- function(px) {
  ((1 - px) / 1.10 + (px * q_x1) / (1.10^2)) / (1 + px / 1.10) - P
}

px <- uniroot(f, interval = c(0, 1))$root

# Losses by outcome
L1 <- 400 / 1.10 - 400 * P
L2 <- 400 / (1.10^2) - 400 * P - (400 * P) / 1.10

print(round(c(
  q_x1 = q_x1,
  px = px,
  loss_first_year = L1,
  loss_second_year = L2,
  probability_loss_lt_190 = px
), 6))

# Plot the loss outcomes
barplot(
  c(L1, L2, 0),
  names.arg = c("Death in year 1", "Death in year 2", "Survive 2 years"),
  main = "Exercise 10-9: Loss at Issue by Outcome",
  ylab = "Loss"
)
abline(h = 190, lty = 2)





# R Check for Exercise 10-20
# Given Var(L) = 0.20, A2bar_x = 0.30, and Abar_{x+20} = 0.70,
# solve for Abar_x and compute {}_{20}V(barA_x).

var_L <- 0.20
A2bar_x <- 0.30
Abar_x20 <- 0.70

# Solve:
# (A2bar_x - Abar_x^2) / (1 - Abar_x)^2 = var_L
f <- function(Abar_x) {
  (A2bar_x - Abar_x^2) / (1 - Abar_x)^2 - var_L
}

root <- uniroot(f, interval = c(1e-8, 0.999999))$root

reserve_20 <- 1 - (1 - Abar_x20) / (1 - root)

print(round(c(
  Abar_x = root,
  reserve_20 = reserve_20
), 6))

# Plot objective function
grid <- seq(0.01, 0.95, by = 0.001)
obj <- sapply(grid, f)

plot(
  grid, obj, type = "l", lwd = 2,
  xlab = "Abar_x",
  ylab = "Variance equation residual",
  main = "Exercise 10-20: Solving for Abar_x"
)
abline(h = 0, lty = 2)
abline(v = root, lty = 2)




# R Check for Exercise 10-22
# Compute {}_{10}V(barA_40) under a uniform survival model with omega = 100 and i = 0.05.

library(mqriskR)

x <- 40
t <- 10
i <- 0.05
omega <- 100

reserve_pkg <- tVbarAbarx(
  x = x,
  t = t,
  i = i,
  model = "uniform",
  omega = omega
)

reserve_ratio <- 1 - abarx(
  x = x + t,
  i = i,
  model = "uniform",
  omega = omega
) / abarx(
  x = x,
  i = i,
  model = "uniform",
  omega = omega
)

print(round(c(
  reserve_pkg = reserve_pkg,
  reserve_ratio = reserve_ratio
), 6))

# Plot reserve as a function of t
t_grid <- seq(0, 20, by = 0.25)
reserve_grid <- tVbarAbarx(
  x = x,
  t = t_grid,
  i = i,
  model = "uniform",
  omega = omega
)

plot(
  t_grid, reserve_grid, type = "l", lwd = 2,
  xlab = "Duration t",
  ylab = "Fully continuous reserve",
  main = "Exercise 10-22: Continuous Reserve Path"
)
abline(v = t, lty = 2)
abline(h = reserve_pkg, lty = 2)



# R Check for Exercise 10-25
# Compute 1000({}_5V_35^(4) - {}_5V_35) using Example 10.8 relationship.

A35 <- 0.17092
V5 <- 0.04471
i <- 0.05

# Interest conversion helper from package
library(mqriskR)

calc_diff <- function(m) {
  ic <- interest_convert(i = i, m = m)
  im <- ic$im
  dm <- (im / m) / (1 + im / m)

  beta_m <- (i - im) / (im * dm)
  alpha_m <- (i * (i / (1 + i))) / (im * dm)

  d <- i / (1 + i)
  adot35 <- (1 - A35) / d
  adot35_m <- alpha_m * adot35 - beta_m
  P35_m <- A35 / adot35_m

  1000 * beta_m * P35_m * V5
}

m_vals <- c(2, 4, 12)
diff_vals <- sapply(m_vals, calc_diff)

print(round(c(
  diff_m2 = diff_vals[1],
  diff_m4 = diff_vals[2],
  diff_m12 = diff_vals[3]
), 6))

barplot(
  diff_vals,
  names.arg = paste0("m=", m_vals),
  main = "Exercise 10-25: Reserve Difference by Payment Frequency",
  ylab = "1000(V^(m) - V)"
)



# R Check for Exercise 10-30
# Verify G^M = (q_assumed - q_actual)(1 - V_{t+1})

Vt <- 0.085697
Vt1 <- 0.096115
P <- 0.008013
i_assumed <- 0.06
q_assumed <- 0.00356
q_actual <- 0.00300

GM_formula <- (Vt + P) * (1 + i_assumed) - (q_actual + (1 - q_actual) * Vt1)
GM_short <- (q_assumed - q_actual) * (1 - Vt1)

print(round(c(
  GM_formula = GM_formula,
  GM_short = GM_short,
  difference = GM_formula - GM_short
), 10))

# Sensitivity in actual mortality
q_actual_grid <- seq(0.002, 0.005, by = 0.00005)
GM_grid <- (q_assumed - q_actual_grid) * (1 - Vt1)

plot(
  q_actual_grid, GM_grid, type = "l", lwd = 2,
  xlab = "Actual mortality q'_x+t",
  ylab = "Mortality gain",
  main = "Exercise 10-30: Mortality Gain Sensitivity"
)
abline(v = q_actual, lty = 2)
abline(h = GM_short, lty = 2)



# R Check for Exercise 10-31
# Verify G^I = (V_t + P)(i_actual - i_assumed)

Vt <- 0.085697
P <- 0.008013
Vt1 <- 0.096115
i_assumed <- 0.06
i_actual <- 0.065
q_assumed <- 0.00356

GI_formula <- (Vt + P) * (1 + i_actual) - (q_assumed + (1 - q_assumed) * Vt1)
GI_short <- (Vt + P) * (i_actual - i_assumed)

print(round(c(
  GI_formula = GI_formula,
  GI_short = GI_short,
  difference = GI_formula - GI_short
), 10))

# Sensitivity in actual interest
i_actual_grid <- seq(0.04, 0.08, by = 0.0005)
GI_grid <- (Vt + P) * (i_actual_grid - i_assumed)

plot(
  i_actual_grid, GI_grid, type = "l", lwd = 2,
  xlab = "Actual interest i'",
  ylab = "Interest gain",
  main = "Exercise 10-31: Interest Gain Sensitivity"
)
abline(v = i_actual, lty = 2)
abline(h = GI_short, lty = 2)
