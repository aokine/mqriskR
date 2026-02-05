# Chapter 3 — Discrete-time Markov chains (non-homogeneous example)

P0 <- matrix(c(
  0.60, 0.40,
  0.70, 0.30
), nrow = 2, byrow = TRUE)

P1 <- matrix(c(
  0.50, 0.50,
  0.80, 0.20
), nrow = 2, byrow = TRUE)

# Process begins in State 1 at time 0 => pi0 = (0,1)
pi0 <- c(0, 1)

pi1 <- pi0 %*% P0
pi2 <- pi1 %*% P1

pi1
pi2

# Probability in State 0 at time 2:
pi2[1]
