# R Check for Example 12.2
# Verify {}_{2|}q_{xy} for independent lives

qx  <- 0.08
qy  <- 0.10
qx1 <- 0.09
qy1 <- 0.15
qx2 <- 0.10
qy2 <- 0.20

px  <- 1 - qx
py  <- 1 - qy
px1 <- 1 - qx1
py1 <- 1 - qy1
px2 <- 1 - qx2
py2 <- 1 - qy2

p2x <- px * px1
p2y <- py * py1

q_xy_age2 <- 1 - px2 * py2
value <- p2x * p2y * q_xy_age2

round(c(
  p2x = p2x,
  p2y = p2y,
  q_xy_age2 = q_xy_age2,
  `2|q_xy` = value
), 6)

# Small sensitivity plot: vary q_{x+2}
qx2_grid <- seq(0.02, 0.30, by = 0.01)
value_grid <- p2x * p2y * (1 - (1 - qx2_grid) * py2)

plot(
  qx2_grid, value_grid, type = "l", lwd = 2,
  xlab = expression(q[x+2]),
  ylab = expression({2*"|"}*q[xy]),
  main = "Example 12.2: Third-Year Joint-Life Failure Probability"
)
abline(v = qx2, lty = 2)
abline(h = value, lty = 2)


# R Check for Example 12.3
# Verify e̊_xy = 1 / (lambda_x + lambda_y)

lambda_x <- 0.03
lambda_y <- 0.05

closed_form <- 1 / (lambda_x + lambda_y)

numeric_val <- integrate(
  function(t) exp(-(lambda_x + lambda_y) * t),
  lower = 0, upper = Inf
)$value

round(c(
  closed_form = closed_form,
  numeric = numeric_val
), 6)

# Sensitivity plot
lambda_sum <- seq(0.02, 0.20, by = 0.002)
e_grid <- 1 / lambda_sum

plot(
  lambda_sum, e_grid, type = "l", lwd = 2,
  xlab = expression(lambda[x] + lambda[y]),
  ylab = expression(ring(e)[xy]),
  main = "Example 12.3: Joint-Life Complete Expectation"
)
abline(v = lambda_x + lambda_y, lty = 2)
abline(h = closed_form, lty = 2)




# R Check for Example 12.5
# Verify e̊_{overline{50:60}} under uniform model with omega = 100

library(mqriskR)

x <- 50
y <- 60
omega <- 100

ex <- ex_complete(x = x, model = "uniform", omega = omega)
ey <- ex_complete(x = y, model = "uniform", omega = omega)

exy <- integrate(
  function(t) {
    p1 <- ifelse(t <= omega - x, (omega - x - t) / (omega - x), 0)
    p2 <- ifelse(t <= omega - y, (omega - y - t) / (omega - y), 0)
    p1 * p2
  },
  lower = 0, upper = omega - y
)$value

e_last <- ex + ey - exy

round(c(
  ex = ex,
  ey = ey,
  exy = exy,
  e_last = e_last
), 6)

# Plot joint-life and last-survivor survival functions
t_grid <- seq(0, 40, by = 0.25)
p_joint <- ifelse(t_grid <= 40,
                  ((50 - t_grid) / 50) * ((40 - t_grid) / 40), 0)
p_last <- ifelse(t_grid <= 40,
                 (50 - t_grid) / 50 + (40 - t_grid) / 40 - p_joint,
                 ifelse(t_grid <= 50, (50 - t_grid) / 50, 0))

plot(
  t_grid, p_joint, type = "l", lwd = 2,
  xlab = "t",
  ylab = "Survival probability",
  main = "Example 12.5: Joint-Life and Last-Survivor Survival"
)
lines(t_grid, p_last, lwd = 2, lty = 2)
legend("topright", legend = c("Joint-life", "Last-survivor"), lty = c(1, 2), lwd = 2, bty = "n")


# R Check for Example 12.13
# Joint density f(x,y) = (x+y)/27 on (0,3)x(0,3)

fxy <- function(tx, ty) {
  ifelse(tx > 0 & tx < 3 & ty > 0 & ty < 3, (tx + ty) / 27, 0)
}

fx <- function(tx) {
  sapply(tx, function(z) integrate(function(ty) fxy(z, ty), 0, 3)$value)
}

fy <- function(ty) {
  sapply(ty, function(z) integrate(function(tx) fxy(tx, z), 0, 3)$value)
}

EX <- integrate(function(tx) tx * fx(tx), 0, 3)$value
EY <- integrate(function(ty) ty * fy(ty), 0, 3)$value

EX2 <- integrate(function(tx) tx^2 * fx(tx), 0, 3)$value
EY2 <- integrate(function(ty) ty^2 * fy(ty), 0, 3)$value

EXY <- integrate(
  function(tx) sapply(tx, function(xx) {
    integrate(function(ty) xx * ty * fxy(xx, ty), 0, 3)$value
  }),
  0, 3
)$value

VX <- EX2 - EX^2
VY <- EY2 - EY^2
COV <- EXY - EX * EY
CORR <- COV / sqrt(VX * VY)

round(c(
  EX = EX,
  EY = EY,
  EXY = EXY,
  VX = VX,
  VY = VY,
  COV = COV,
  CORR = CORR
), 6)

# Plot marginal density of T_x
tx_grid <- seq(0, 3, by = 0.01)
plot(
  tx_grid, fx(tx_grid), type = "l", lwd = 2,
  xlab = expression(t[x]),
  ylab = expression(f[x](t[x])),
  main = "Example 12.13: Marginal Density of T_x"
)



# R Check for Example 12.15
# Device fails on first failure of its two components

fxy <- function(tx, ty) {
  ifelse(tx > 0 & tx < 3 & ty > 0 & ty < 3, (tx + ty) / 27, 0)
}

p_joint <- integrate(
  function(tx) sapply(tx, function(xx) {
    integrate(function(ty) fxy(xx, ty), 1, 3)$value
  }),
  1, 3
)$value

q_first_month <- 1 - p_joint

round(c(
  p_joint = p_joint,
  q_first_month = q_first_month
), 6)

# Small visual check over n in [0,3]
n_grid <- seq(0, 3, by = 0.05)
q_grid <- sapply(n_grid, function(n) {
  p_joint_n <- integrate(
    function(tx) sapply(tx, function(xx) {
      integrate(function(ty) fxy(xx, ty), n, 3)$value
    }),
    n, 3
  )$value
  1 - p_joint_n
})

plot(
  n_grid, q_grid, type = "l", lwd = 2,
  xlab = "n",
  ylab = expression(n*q[xy]),
  main = "Example 12.15: Joint-Life Failure Probability by Time n"
)
abline(v = 1, lty = 2)
abline(h = q_first_month, lty = 2)



# R Check for Example 12.16
# Compare last-survivor APV under common shock for exponential components

delta <- 0.05
lambda_x <- 0.03
lambda_y <- 0.04

A_last_common_shock <- function(lambda_common) {
  A_x  <- (lambda_x + lambda_common) / (lambda_x + lambda_common + delta)
  A_y  <- (lambda_y + lambda_common) / (lambda_y + lambda_common + delta)
  A_xy <- (lambda_x + lambda_y + lambda_common) /
    (lambda_x + lambda_y + lambda_common + delta)
  A_x + A_y - A_xy
}

A_last_independent <- function() {
  A_x  <- lambda_x / (lambda_x + delta)
  A_y  <- lambda_y / (lambda_y + delta)
  A_xy <- (lambda_x + lambda_y) / (lambda_x + lambda_y + delta)
  A_x + A_y - A_xy
}

lambda_common <- 0.02

base_common <- A_last_common_shock(lambda_common)
base_indep <- A_last_independent()

round(c(
  common_shock = base_common,
  independence = base_indep,
  difference = base_common - base_indep
), 6)

# Sensitivity plot
lambda_grid <- seq(0, 0.08, by = 0.001)
A_grid <- sapply(lambda_grid, A_last_common_shock)

plot(
  lambda_grid, A_grid, type = "l", lwd = 2,
  xlab = expression(lambda[c]),
  ylab = expression(bar(A)[bar(xy)]),
  main = "Example 12.16: Last-Survivor APV under Common Shock"
)
abline(v = lambda_common, lty = 2)
abline(h = base_common, lty = 2)
abline(h = base_indep, lty = 3)






# R Check for Exercise 12-2
# Verify the PDF of T_xy at t = 0.50

qx <- 0.080
qy <- 0.004

p_xy <- function(t) (1 - qx * t^2) * (1 - qy * t^2)
f_xy <- function(t) 2 * 0.084 * t - 0.00128 * t^3

t0 <- 0.50
pdf_exact <- f_xy(t0)

# numerical derivative check
h <- 1e-6
pdf_num <- -(p_xy(t0 + h) - p_xy(t0 - h)) / (2 * h)

round(c(
  pdf_exact = pdf_exact,
  pdf_numeric = pdf_num
), 6)

# plot the joint-life survival function and density
t_grid <- seq(0, 1, by = 0.01)

plot(
  t_grid, p_xy(t_grid), type = "l", lwd = 2,
  xlab = "t",
  ylab = expression(t*p[xy]),
  main = "Exercise 12-2: Joint-Life Survival Function"
)

plot(
  t_grid, f_xy(t_grid), type = "l", lwd = 2,
  xlab = "t",
  ylab = expression(f[xy](t)),
  main = "Exercise 12-2: Joint-Life Density"
)
abline(v = t0, lty = 2)
abline(h = pdf_exact, lty = 2)




# R Check for Exercise 12-4
# Verify e̊_{65:55} for independent nonsmoker age 65 and smoker age 55

p65N <- function(t) ifelse(t >= 0 & t <= 10, 1 - t / 10, 0)
p55S <- function(t) ifelse(t >= 0 & t <= 20, (1 - t / 20)^2, 0)

p_joint <- function(t) p65N(t) * p55S(t)

e_joint <- integrate(p_joint, lower = 0, upper = 10)$value

round(c(e_joint = e_joint), 6)

t_grid <- seq(0, 10, by = 0.05)

plot(
  t_grid, p_joint(t_grid), type = "l", lwd = 2,
  xlab = "t",
  ylab = expression(t*p[65:55]),
  main = "Exercise 12-4: Joint-Life Survival Curve"
)
abline(h = 0, lty = 3)



# R Check for Exercise 12-6
# Verify {}_5 q_{overline{80:85}}

library(mqriskR)

x <- 80
y <- 85
t <- 5

q_last <- tqxybar(
  x = x,
  y = y,
  t = t,
  model = "uniform",
  omega = 100
)

qx <- tqx(x = x, t = t, model = "uniform", omega = 100)
qy <- tqx(x = y, t = t, model = "uniform", omega = 100)

product_check <- qx * qy

round(c(
  q_last = q_last,
  qx = qx,
  qy = qy,
  product_check = product_check
), 6)

# sensitivity to horizon
t_grid <- seq(0, 10, by = 0.1)
q_last_grid <- sapply(t_grid, function(tt) {
  tqxybar(x = x, y = y, t = tt, model = "uniform", omega = 100)
})

plot(
  t_grid, q_last_grid, type = "l", lwd = 2,
  xlab = "t",
  ylab = expression(t*q[bar(80:85)]),
  main = "Exercise 12-6: Last-Survivor Failure Probability"
)
abline(v = t, lty = 2)
abline(h = q_last, lty = 2)



# R Check for Exercise 12-11
# Verify Pr(T_I < T_II)

S_I <- function(x, a = 1.80) ((9 - x) / 9)^a
S_II <- function(x, b = 1.50) ((9 - x) / 9)^b
f_I <- function(x, a = 1.80) ((9 - x) / 9)^a * a / (9 - x)

prob <- integrate(function(x) f_I(x) * S_II(x), lower = 0, upper = 9)$value
round(c(prob = prob), 6)

# sensitivity plot in the coefficient for machine I
a_grid <- seq(0.8, 2.8, by = 0.02)
prob_grid <- sapply(a_grid, function(a) {
  integrate(function(x) f_I(x, a = a) * S_II(x), lower = 0, upper = 9)$value
})

plot(
  a_grid, prob_grid, type = "l", lwd = 2,
  xlab = "Hazard coefficient for Machine I",
  ylab = expression(Pr(T[I] < T[II])),
  main = "Exercise 12-11: Contingent Failure Probability"
)
abline(v = 1.80, lty = 2)
abline(h = prob, lty = 2)



# R Check for Exercise 12-15
# Verify P_xy from the given premium information

d <- 0.06
Px <- 0.10
Pxybar <- 0.06

adotx <- 1 / (Px + d)
Ax <- 1 - d * adotx

# Solve for adotxy from
# (.625 + .625 - 1 + .06*adotxy) / (6.25 + 6.25 - adotxy) = .06

adotxy <- 50 / 12
Axy <- 1 - d * adotxy
Pxy <- Axy / adotxy

round(c(
  adotx = adotx,
  Ax = Ax,
  adotxy = adotxy,
  Axy = Axy,
  Pxy = Pxy,
  Pxybar = Pxybar
), 6)

# compare P_xy as a function of adotxy
adotxy_grid <- seq(3.5, 6.0, by = 0.01)
Pxy_grid <- (1 - d * adotxy_grid) / adotxy_grid

plot(
  adotxy_grid, Pxy_grid, type = "l", lwd = 2,
  xlab = expression(ddot(a)[xy]),
  ylab = expression(P[xy]),
  main = "Exercise 12-15: Joint-Life Premium vs Annuity Value"
)
abline(v = adotxy, lty = 2)
abline(h = Pxy, lty = 2)




# R Check for Exercise 12-32
# Compare independence and common-shock APVs

delta <- 0.05
mu_total <- 0.06
lambda_common <- 0.02

# Independence case
A_indep <- mu_total / (mu_total + delta) +
  mu_total / (mu_total + delta) -
  (2 * mu_total) / (2 * mu_total + delta)

# Common-shock case
mu_star <- mu_total - lambda_common
mu_joint_common <- mu_star + mu_star + lambda_common

A_common <- mu_total / (mu_total + delta) +
  mu_total / (mu_total + delta) -
  mu_joint_common / (mu_joint_common + delta)

increase <- A_common - A_indep

round(c(
  A_indep = A_indep,
  A_common = A_common,
  increase = increase
), 6)

# sensitivity plot in lambda_common
lambda_grid <- seq(0, 0.05, by = 0.001)
A_common_grid <- sapply(lambda_grid, function(lam) {
  mu_star <- mu_total - lam
  mu_joint <- mu_star + mu_star + lam
  mu_total / (mu_total + delta) +
    mu_total / (mu_total + delta) -
    mu_joint / (mu_joint + delta)
})

plot(
  lambda_grid, A_common_grid, type = "l", lwd = 2,
  xlab = expression(lambda),
  ylab = expression(bar(A)[bar(xy)]),
  main = "Exercise 12-32: Common Shock and Last-Survivor APV"
)
abline(v = lambda_common, lty = 2)
abline(h = A_common, lty = 2)
abline(h = A_indep, lty = 3)




# R Check for Exercise 12-33
# Verify {}_5 p_xy under common shock

p_x <- 0.96
p_y <- 0.97
lambda <- 0.01
t <- 5

p_x_star <- p_x * exp(lambda)
p_y_star <- p_y * exp(lambda)

p5_xy <- (p_x_star^t) * (p_y_star^t) * exp(-lambda * t)

round(c(
  p_x_star = p_x_star,
  p_y_star = p_y_star,
  p5_xy = p5_xy
), 6)

# sensitivity plot
lambda_grid <- seq(0, 0.03, by = 0.0005)
p5_grid <- sapply(lambda_grid, function(lam) {
  pxs <- p_x * exp(lam)
  pys <- p_y * exp(lam)
  (pxs^t) * (pys^t) * exp(-lam * t)
})

plot(
  lambda_grid, p5_grid, type = "l", lwd = 2,
  xlab = expression(lambda),
  ylab = expression(5*p[xy]),
  main = "Exercise 12-33: Five-Year Joint Survival under Common Shock"
)
abline(v = lambda, lty = 2)
abline(h = p5_xy, lty = 2)
