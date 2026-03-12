# R Check for Example 8.2
library(mqriskR)

x <- 40
i <- 0.05
omega <- 100

v <- 1 / (1 + i)
d <- i / (1 + i)

Ax_val <- Ax(x = x, i = i, model = "uniform", omega = omega)
ax_val <- ax(x = x, i = i, model = "uniform", omega = omega)

lhs <- Ax_val
rhs <- v - d * ax_val

round(c(lhs = lhs, rhs = rhs, difference = lhs - rhs), 10)



# R Check for Example 8.3
i <- 0.05
ip <- (1 + i)^2 - 1

v  <- 1 / (1 + i)
vp <- 1 / (1 + ip)

d  <- i / (1 + i)

lx_vals <- c(l95 = 100, l96 = 70, l97 = 40, l98 = 20, l99 = 4, l100 = 0)

surv <- c(
  lx_vals["l96"] / lx_vals["l95"],
  lx_vals["l97"] / lx_vals["l95"],
  lx_vals["l98"] / lx_vals["l95"],
  lx_vals["l99"] / lx_vals["l95"]
)

a95  <- sum(v^(1:4)  * surv)
a95_2 <- sum(vp^(1:4) * surv)

A95  <- v  - d * a95
dp   <- ip / (1 + ip)
A95_2 <- vp - dp * a95_2

var_y <- (A95_2 - A95^2) / d^2

round(c(
  a95 = a95,
  a95_2 = a95_2,
  A95 = A95,
  A95_2 = A95_2,
  varY95 = var_y
), 4)


# R Check for Example 8.7
i <- 0.06
n <- 10
coupon <- 40

v <- 1 / (1 + i)
surv <- 0.98^(1:n)

epv <- coupon * sum(v^(1:n) * surv)

round(epv, 2)


# R Check for Example 8.8
i <- 0.05
ip <- 0.1025

d  <- i / (1 + i)
dp <- ip / (1 + ip)

adot_a5  <- 4.13038
adot2_a5 <- 3.79209

Axn  <- 1 - d  * adot_a5
A2xn <- 1 - dp * adot2_a5

var_temp_due <- (A2xn - Axn^2) / d^2

round(c(
  Axn = Axn,
  A2xn = A2xn,
  variance = var_temp_due
), 5)


# R Check for Example 8.14
q80_1990 <- 0.05105
improvement <- 0.005
reduction_factor <- 1 - improvement

q80_2010 <- q80_1990 * reduction_factor^20
p80_2010 <- 1 - q80_2010

round(c(
  q80_2010 = q80_2010,
  p80_2010 = p80_2010
), 5)



# R Check for Example 8.17 b (i)

library(mqriskR)

delta <- 0.04
mu <- 0.01
i <- exp(delta) - 1

# -------------------------------
# Exact values from formulas in the text
# -------------------------------

a_cont_exact <- 1 / (mu + delta)
a_due_exact  <- 1 / (1 - exp(-(mu + delta)))

d <- 1 - exp(-delta)

udd_formula <- (i * d / delta^2) * a_due_exact - (i - delta) / delta^2
w2_formula  <- a_due_exact - 0.5
w3_formula  <- a_due_exact - 0.5 - (mu + delta) / 12


# -------------------------------
# Using package functions
# -------------------------------

abar_pkg <- abarx(
  x = 50,
  i = i,
  model = "exponential",
  lambda = mu
)

adot_pkg <- adotx(
  x = 50,
  i = i,
  model = "exponential",
  lambda = mu
)

udd_pkg <- abarx_udd(
  x = 50,
  i = i,
  model = "exponential",
  lambda = mu
)

w2_pkg <- abarx_woolhouse2(
  x = 50,
  i = i,
  model = "exponential",
  lambda = mu
)

w3_pkg <- abarx_woolhouse3(
  x = 50,
  i = i,
  model = "exponential",
  lambda = mu
)


# -------------------------------
# Compare results
# -------------------------------

round(c(
  abar_exact = a_cont_exact,
  abar_pkg   = abar_pkg,

  adot_exact = a_due_exact,
  adot_pkg   = adot_pkg,

  udd_formula = udd_formula,
  udd_pkg     = udd_pkg,

  woolhouse2_formula = w2_formula,
  woolhouse2_pkg     = w2_pkg,

  woolhouse3_formula = w3_formula,
  woolhouse3_pkg     = w3_pkg
), 5)



# R Check for Exercise 8-4

adotx_val <- 10
A2adotx_val <- 6
i <- 1 / 24

d <- i / (1 + i)
ip <- (1 + i)^2 - 1
dp <- ip / (1 + ip)

Ax_val <- 1 - d * adotx_val
A2x_val <- 1 - dp * A2adotx_val

var_val <- (A2x_val - Ax_val^2) / d^2

round(c(
  d = d,
  ip = ip,
  dp = dp,
  Ax = Ax_val,
  A2x = A2x_val,
  variance = var_val
), 5)



# R Check for Exercise 8-7

i <- 0.05
v <- 1 / (1 + i)

px <- 0.99
px1 <- 0.95
px1_new <- 0.98
adotx1 <- 6.951

adotx2 <- ((adotx1 - 1) * (1 + i)) / px1

adotx_old <- 1 + v * px * (1 + v * px1 * adotx2)
adotx_new <- 1 + v * px * (1 + v * px1_new * adotx2)

increase <- adotx_new - adotx_old

round(c(
  adotx2 = adotx2,
  old_value = adotx_old,
  new_value = adotx_new,
  increase = increase
), 5)



# R Check for Exercise 8-8

annuity_payment <- 12000
d <- 0.08

B_star <- annuity_payment / d

B_star

B <- seq(100000, 200000, by = 1000)
coef_sq <- (B - B_star)^2

plot(B, coef_sq, type = "l",
     xlab = "Death benefit B",
     ylab = "Proportional variance term",
     main = "Variance minimized when B = 150,000")
abline(v = B_star, lty = 2)
grid()



# R Check for Exercise 8-13

abarx <- 10
A2abarx <- 7.375
varY <- 50

# From the algebra in the text:
# varY = (5.25*delta - 100*delta^2)/delta^2 = 5.25/delta - 100
delta <- 5.25 / 150
Abarx <- 1 - delta * abarx

round(c(delta = delta, Abarx = Abarx), 5)




# R Check for Exercise 8-15

delta <- 0.03
mu <- 0.025

threshold_t <- log(1 - 20 * delta) / (-delta)
prob <- exp(-mu * threshold_t)

round(c(
  threshold_t = threshold_t,
  probability = prob
), 5)




# R Check for Exercise 8-18

delta <- 0.10

Abar_m <- 0.15
Abar_f <- 0.09

varY_m <- 5.00
varY_f <- 4.00

A2bar_m <- varY_m * delta^2 + Abar_m^2
A2bar_f <- varY_f * delta^2 + Abar_f^2

Abar_uncond <- 0.5 * (Abar_m + Abar_f)
A2bar_uncond <- 0.5 * (A2bar_m + A2bar_f)

varY_uncond <- (A2bar_uncond - Abar_uncond^2) / delta^2

round(c(
  A2bar_m = A2bar_m,
  A2bar_f = A2bar_f,
  Abar_uncond = Abar_uncond,
  A2bar_uncond = A2bar_uncond,
  varY_uncond = varY_uncond
), 5)



# R Check for Exercise 8-27

y <- c(1.00, 1.87, 2.62)
p <- c(0.10, 0.90 - 0.81, 0.81)

EY <- sum(y * p)
EY2 <- sum(y^2 * p)
VarY <- EY2 - EY^2

round(c(EY = EY, EY2 = EY2, VarY = VarY), 5)


# R Check for Exercise 8-33

adot65 <- 9.90
A35_30 <- 0.21
A35_30_1 <- 0.07

E30_35 <- A35_30 - A35_30_1
NSP <- (E30_35 * adot65) / (1 - A35_30_1)

round(c(E30_35 = E30_35, NSP = NSP), 5)


# R Check for Exercise 8-48

v <- 0.90
y <- c(
  2,
  2 + 3 * v,
  2 + 3 * v + 4 * v^2
)

p <- c(0.20, 0.20, 0.60)

EY <- sum(y * p)
EY2 <- sum(y^2 * p)
VarY <- EY2 - EY^2

round(c(EY = EY, EY2 = EY2, VarY = VarY), 5)


# R Check for Exercise 8-52

qx2000 <- c(0.000475, 0.000514, 0.000554, 0.000598, 0.000648)
AAx <- c(0.011, 0.012, 0.013, 0.014, 0.015)
years_forward <- c(8, 9, 10, 11, 12)

px_proj <- 1 - qx2000 * (1 - AAx)^years_forward
q5 <- 1 - prod(px_proj)

round(c(px_proj, q5 = q5), 6)

