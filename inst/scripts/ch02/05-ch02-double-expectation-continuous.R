# Chapter 2 Script 05: Verification of the double expectation theorem (continuous example)

# Setup:
# X ~ Unif(0, 12)
# Y | X=x ~ Unif(0, x)

EX <- 6
VarX <- 12
EX2 <- VarX + EX^2  # E[X^2]

# From conditional moments:
# E[Y|X] = X/2
# Var(Y|X) = X^2/12
EY <- 0.5 * EX
VarY <- (1/12) * EX2 + (1/4) * VarX

cat("E[Y] =", EY, "\n")
cat("Var(Y) =", VarY, "\n")

# Optional Monte Carlo check (kept simple, base R)
set.seed(1)
n <- 200000
X <- runif(n, min = 0, max = 12)
Y <- runif(n, min = 0, max = X)  # vectorized: max can be a vector

cat("Simulated E[Y] ~", mean(Y), "\n")
cat("Simulated Var(Y) ~", var(Y), "\n")
