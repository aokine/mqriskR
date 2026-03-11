## Listing 7.1: Basic discrete insurance functions from a parametric survival model
library(mqriskR)

x <- 40
n <- 20
i <- 0.05

# Use a simple uniform survival model with limiting age omega = 100
Ax_val    <- Ax(x = x, i = i, model = "uniform", omega = 100)
Axn1_val  <- Axn1(x = x, n = n, i = i, model = "uniform", omega = 100)
nAx_val   <- nAx(x = x, n = n, i = i, model = "uniform", omega = 100)
nEx_val   <- nEx(x = x, n = n, i = i, model = "uniform", omega = 100)
Axn_val   <- Axn(x = x, n = n, i = i, model = "uniform", omega = 100)

out <- c(
  Ax   = Ax_val,
  Axn1 = Axn1_val,
  nAx  = nAx_val,
  nEx  = nEx_val,
  Axn  = Axn_val
)

print(round(out, 6))

# Check the endowment decomposition
round(Axn_val - (Axn1_val + nEx_val), 10)




## Listing 7.2: Comparison of annual, continuous, and \(m\)-thly insurance APVs
library(mqriskR)

x <- 40
i <- 0.05

A_annual <- Ax(x = x, i = i, model = "uniform", omega = 100)
A_cont   <- Abarx(x = x, i = i, model = "uniform", omega = 100)

m_grid <- c(2, 4, 12, 24)
A_m <- sapply(m_grid, function(mm) {
  Ax_m(x = x, i = i, m = mm, model = "uniform", omega = 100)
})

vals <- c(annual = A_annual, A_m, continuous = A_cont)
print(round(vals, 6))

