# Chapter 2 Script 02: Poisson probabilities and simulation in R

lambda <- 3

# Probabilities for X ~ Pois(lambda)
p_X_eq_2  <- dpois(2, lambda)
p_X_le_2  <- ppois(2, lambda)

cat("lambda =", lambda, "\n")
cat("P(X = 2) =", p_X_eq_2, "\n")
cat("P(X <= 2) =", p_X_le_2, "\n")

# Simulation check: mean and variance should be close to lambda
set.seed(1)
x <- rpois(10000, lambda)

cat("Simulated mean ~", mean(x), "\n")
cat("Simulated var  ~", var(x), "\n")
