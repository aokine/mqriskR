# Chapter 3 — Discrete-time Markov chains (homogeneous examples)

P <- matrix(c(
  0.80, 0.20, 0.00,
  0.30, 0.60, 0.10,
  0.00, 0.00, 1.00
), nrow = 3, byrow = TRUE)

# State labels are 0,1,2 in the text; in R indices are 1,2,3.
# If process known to be in State 1 at time n, then pi_n = (0,1,0)
pi_n <- c(0, 1, 0)

pi_n1 <- pi_n %*% P
pi_n2 <- pi_n %*% (P %*% P)

pi_n1
pi_n2

# r-step transition probabilities are in P^r.
# Example: probability from state i to j in r steps:
r <- 2
i <- 1  # State 1 in text
j <- 0  # State 0 in text
P_r <- P %*% P
p_ij_r <- P_r[i + 1, j + 1]
p_ij_r
