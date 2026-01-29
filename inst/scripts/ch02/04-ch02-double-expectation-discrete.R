# Chapter 2 Script 04: Verification of the double expectation theorem (discrete example)

# X in {0,1,2}, Y in {1,2}
x_vals <- c(0, 1, 2)
y_vals <- c(1, 2)

# Joint pmf table p(y, x) with rows = y, cols = x
p <- matrix(
  c(0.10, 0.20, 0.30,
    0.10, 0.10, 0.20),
  nrow = length(y_vals),
  byrow = TRUE
)

stopifnot(abs(sum(p) - 1) < 1e-12)

# Marginal distribution of X
px <- colSums(p)
EX <- sum(x_vals * px)
VarX <- sum((x_vals^2) * px) - EX^2

# Marginal distribution of Y
py <- rowSums(p)

# Conditional pmf of X given Y
cond_px_given_y <- sweep(p, 1, py, "/")

# Conditional moments of X given Y
EX_given_y  <- as.numeric(cond_px_given_y %*% x_vals)
EX2_given_y <- as.numeric(cond_px_given_y %*% (x_vals^2))
Var_given_y <- EX2_given_y - EX_given_y^2

# Double expectation theorem
EX_alt <- sum(EX_given_y * py)
Var_alt <- sum(Var_given_y * py) + (sum((EX_given_y^2) * py) - EX_alt^2)

cat("E[X] direct =", EX, "\n")
cat("E[X] via double expectation =", EX_alt, "\n")
cat("Var(X) direct =", VarX, "\n")
cat("Var(X) via law of total variance =", Var_alt, "\n")
