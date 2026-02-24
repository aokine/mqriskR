## Listing 2.1: Discrete distribution utilities (PMF -> CDF, mean, variance)
x <- c(0, 1, 2)
p <- c(0.20, 0.30, 0.50)  # probabilities must sum to 1

# CDF values at support points (after ordering by x)
o <- order(x)
x <- x[o]; p <- p[o]
F <- cumsum(p)

# Mean and variance
EX  <- sum(x * p)
EX2 <- sum((x^2) * p)
VarX <- EX2 - EX^2

EX
VarX
F





## Listing 2.2: R distribution functions: density (d*), CDF (p*), quantile (q*), simulation (r*)}, label={lst:ch2_r_distribution_convention}]
# Normal distribution: X ~ N(mu, sigma^2)
mu <- 10
sigma <- 2

# Probability: P(X <= 12)
pnorm(12, mean = mu, sd = sigma)

# 95th percentile of X
qnorm(0.95, mean = mu, sd = sigma)

# Simulate 5 observations
rnorm(5, mean = mu, sd = sigma)


## Listing 2.3: Poisson probabilities and simulation
lambda <- 3

# P(X = 2) and P(X <= 2) for X ~ Pois(lambda)
dpois(2, lambda)
ppois(2, lambda)

# Simulate and compare sample mean/variance to lambda
set.seed(1)
x <- rpois(10000, lambda)
mean(x)
var(x)


## Listing 2.4: Normal probabilities without tables using pnorm and qnorm
# Standard normal probabilities
pnorm(1.96)                 # P(Z <= 1.96)
1 - pnorm(1.96)             # P(Z > 1.96)

# Two-sided probability P(|Z| <= 1.96)
pnorm(1.96) - pnorm(-1.96)

# Find z such that P(Z <= z) = 0.975
qnorm(0.975)


## Listing 2.5 Verification of the double expectation theorem (discrete example)
# Joint pmf table for Y in {1,2} and X in {0,1,2}
x_vals <- c(0, 1, 2)
y_vals <- c(1, 2)

p <- matrix(
  c(0.10, 0.20, 0.30,
    0.10, 0.10, 0.20),
  nrow = length(y_vals), byrow = TRUE
)

# Marginal of X
px <- colSums(p)
EX <- sum(x_vals * px)
VarX <- sum((x_vals^2) * px) - EX^2

# Conditional moments of X given Y
py <- rowSums(p)
cond_px_given_y <- sweep(p, 1, py, "/")

EX_given_y  <- as.numeric(cond_px_given_y %*% x_vals)
EX2_given_y <- as.numeric(cond_px_given_y %*% (x_vals^2))
Var_given_y <- EX2_given_y - EX_given_y^2

# Double expectation theorem
EX_alt <- sum(EX_given_y * py)
Var_alt <- sum(Var_given_y * py) + (sum((EX_given_y^2) * py) - EX_alt^2)

c(EX = EX, EX_alt = EX_alt, VarX = VarX, Var_alt = Var_alt)



## Listing 2.6: Verification of the double expectation theorem (continuous example)
# X ~ Unif(0,12), Y | X=x ~ Unif(0,x)

EX <- 6
VarX <- 12
EX2 <- VarX + EX^2  # E[X^2]

EY <- 0.5 * EX
VarY <- (1/12) * EX2 + (1/4) * VarX

c(EY = EY, VarY = VarY)
\end{lstlisting}
