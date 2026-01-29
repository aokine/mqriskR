# Chapter 2 Script 01: Discrete distribution utilities (PMF -> CDF, mean, variance)

x <- c(0, 1, 2)
p <- c(0.20, 0.30, 0.50)

# Basic validation
stopifnot(length(x) == length(p))
stopifnot(all(p >= 0))
stopifnot(abs(sum(p) - 1) < 1e-12)

# Order support (recommended best practice)
o <- order(x)
x <- x[o]
p <- p[o]

# CDF at support points
F <- cumsum(p)

# Mean and variance
EX  <- sum(x * p)
EX2 <- sum((x^2) * p)
VarX <- EX2 - EX^2

cat("Support x:\n"); print(x)
cat("PMF p:\n"); print(p)
cat("CDF F:\n"); print(F)
cat("E[X] =", EX, "\n")
cat("Var(X) =", VarX, "\n")
