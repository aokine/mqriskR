library(mqriskR)

# Example cashflows (net from investor perspective)
# -X at t=0,1,2,3; +Y at t=4,...,8
X <- 100
Y <- 120

cf <- c(-X, -X, -X, -X, rep(Y, 5))
t  <- c(0, 1, 2, 3, 4:8)

# NPV as a function of i
npv <- function(i) pv_cashflows(cf, t, i)

# Solve yield rate
i_hat <- solve_yield(cf, t, interval = c(-0.5, 1))
i_hat

# Visualize NPV(i)
grid <- seq(-0.2, 0.6, by = 0.005)
vals <- vapply(grid, npv, numeric(1))

plot(grid, vals, xlab = "i", ylab = "NPV(i)", main = "Equation of Value: NPV vs i")
abline(h = 0)
abline(v = i_hat, lty = 2)
