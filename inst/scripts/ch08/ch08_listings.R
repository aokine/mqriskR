## Listing 8.1: Whole life annuity computation from a parametric survival model
library(mqriskR)

x <- 40
i <- 0.05

# Uniform survival model with limiting age omega = 100
ax_val <- ax(x = x, i = i, model = "uniform", omega = 100)

print(ax_val)


## Listing 8.2: Variance of the whole life annuity present value from insurance moments
library(mqriskR)

x <- 95
i <- 0.05
omega <- 100

d <- i / (1 + i)

# First and second insurance moments
A1 <- Ax(x = x, i = i, model = "uniform", omega = omega)
A2 <- A2x(x = x, i = i, model = "uniform", omega = omega)

# Variance using the identity
var_identity <- (A2 - A1^2) / d^2

# Variance using the package variance function
var_pkg <- var_Ax(x = x, i = i, model = "uniform", omega = omega) / d^2

c(identity = var_identity, package = var_pkg)



## Listing 8.3: Temporary life annuity computation
library(mqriskR)

x <- 40
n <- 20
i <- 0.05

axn_val <- axn(x = x, n = n, i = i, model = "uniform", omega = 100)

print(axn_val)



## Listing 4.4: Verification of deferred annuity identities
library(mqriskR)

x <- 40
n <- 20
i <- 0.05

whole <- ax(x = x, i = i, model = "uniform", omega = 100)
temp  <- axn(x = x, n = n, i = i, model = "uniform", omega = 100)

left <- nax(x = x, n = n, i = i, model = "uniform", omega = 100)
right1 <- whole - temp

right2 <- nEx(x = x, n = n, i = i, model = "uniform", omega = 100) *
  ax(x + n, i = i, model = "uniform", omega = 100)

round(c(deferred = left, identity1 = right1, identity2 = right2), 10)




## Listing 8.5: Comparison of annual, m-thly, and continuous whole life annuity values
library(mqriskR)

x <- 40
i <- 0.05
omega <- 100
m_vals <- c(1, 2, 4, 12)

a_annual <- ax(x = x, i = i, model = "uniform", omega = omega)
a_m <- sapply(m_vals, function(m) {
  ax_m(x = x, m = m, i = i, model = "uniform", omega = omega)
})
a_cont <- abarx(x = x, i = i, model = "uniform", omega = omega)

data.frame(
  m = c(m_vals, Inf),
  value = c(a_m, a_cont)
)


## Listing 8.6:Comparison of exact, UDD, and Woolhouse approximations for an \(m\)-thly annuity
library(mqriskR)

x <- 40
i <- 0.05
m <- 12
omega <- 100

exact <- ax_m(x = x, m = m, i = i, model = "uniform", omega = omega)

udd_based <- ax_m_udd(
  x = x,
  m = m,
  i = i,
  model = "uniform",
  omega = omega
)

woolhouse_2 <- ax_m_woolhouse2(
  x = x,
  m = m,
  i = i,
  model = "uniform",
  omega = omega
)

woolhouse_3 <- ax_m_woolhouse3(
  x = x,
  m = m,
  i = i,
  model = "uniform",
  omega = omega
)

c(exact = exact,
  udd = udd_based,
  woolhouse_2 = woolhouse_2,
  woolhouse_3 = woolhouse_3)


## Listing 8.7: Temporary annuity value with and without mortality improvement
library(mqriskR)

x0 <- 60; n <- 5; i <- 0.06
base_year <- 2000
issue_year <- 2008

qx_base <- c(0.0100, 0.0110, 0.0122, 0.0135, 0.0150)
AAx <- c(0.01, 0.01, 0.01, 0.01, 0.01)

a_improved <- axn_improved(
  x0 = x0, n = n, i = i,
  qx_base_vec = qx_base,
  AAx_vec = AAx,
  base_year = base_year,
  issue_year = issue_year
)

p_no_improve <- cumprod(1 - qx_base)
a_no_improve <- sum((1 + i)^(-(1:n)) * p_no_improve)

c(no_improvement = a_no_improve, improvement = a_improved)
