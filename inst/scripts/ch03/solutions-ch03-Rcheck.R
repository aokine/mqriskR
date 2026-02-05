library(expm)

# -------------------------
# Exercise 3-2 R-check
# -------------------------
P0 <- matrix(c(
  0.85, 0.15, 0.00,
  0.00, 0.70, 0.30,
  0.00, 0.00, 1.00
), 3, byrow = TRUE)

P1 <- matrix(c(
  0.90, 0.10, 0.00,
  0.10, 0.70, 0.20,
  0.00, 0.00, 1.00
), 3, byrow = TRUE)

P2 <- matrix(c(
  0.95, 0.05, 0.00,
  0.20, 0.70, 0.10,
  0.00, 0.00, 1.00
), 3, byrow = TRUE)

# For k >= 3
Pk <- matrix(c(
  0.95, 0.05, 0.00,
  0.50, 0.50, 0.00,
  0.00, 0.00, 1.00
), 3, byrow = TRUE)

# start endangered at t=0 => State 1 in book's numbering (0,1,2)
pi0 <- c(0, 1, 0)

# extinction by t=1: pi0 %*% P0 -> prob in State 2
pi1 <- as.numeric(pi0 %*% P0)
p_extinct_by_1 <- pi1[3]

# extinction by t=2
pi2 <- as.numeric(pi1 %*% P1)
p_extinct_by_2 <- pi2[3]

# extinction by t=3
pi3 <- as.numeric(pi2 %*% P2)
p_extinct_by_3 <- pi3[3]

c(by1 = p_extinct_by_1, by2 = p_extinct_by_2, by3 = p_extinct_by_3)

# For this problem extinction can't happen after t=3 (given model structure),
# so P(ever extinct) = P(extinct by t=3)
p_ever_extinct <- p_extinct_by_3
p_ever_extinct


# -------------------------
# Exercise 3-3 R-check (CTMC)
# -------------------------
Q <- matrix(0, 4, 4)
Q[1,2] <- 0.30
Q[1,3] <- 0.50
Q[1,4] <- 0.70
Q[1,1] <- -(0.30 + 0.50 + 0.70)  # -1.50

t <- 1
P_t <- expm(Q * t)

# Probability in State 2 at t=1 given start State 0:
P_t[1, 3]


# -------------------------
# Exercise 3-4 R-check (expected payments)
# -------------------------
P0 <- matrix(c(
  0.60, 0.30, 0.10,
  0.00, 0.00, 1.00,
  0.00, 0.00, 1.00
), 3, byrow = TRUE)

# same as P0 for k=1
P1 <- P0

P2plus <- matrix(c(
  0.00, 0.30, 0.70,
  0.00, 0.00, 1.00,
  0.00, 0.00, 1.00
), 3, byrow = TRUE)

pi <- c(1, 0, 0)

# (a) pay 1 at t=0,1,2,... if in State 0 or 1
pay01 <- c(1, 1, 0)

Epay_a <- 0
Epay_a <- Epay_a + sum(pi * pay01)      # t=0

pi <- as.numeric(pi %*% P0)             # t=1
Epay_a <- Epay_a + sum(pi * pay01)

pi <- as.numeric(pi %*% P1)             # t=2
Epay_a <- Epay_a + sum(pi * pay01)

pi <- as.numeric(pi %*% P2plus)         # t=3
Epay_a <- Epay_a + sum(pi * pay01)

# After t=3, probability in State 0 or 1 is 0, so stop.
Epay_a

# (b) pay 4 at t=1,2,3,... if in State 1
pi <- c(1, 0, 0)
Epay_b <- 0

pi <- as.numeric(pi %*% P0)             # t=1
Epay_b <- Epay_b + 4 * pi[2]

pi <- as.numeric(pi %*% P1)             # t=2
Epay_b <- Epay_b + 4 * pi[2]

pi <- as.numeric(pi %*% P2plus)         # t=3
Epay_b <- Epay_b + 4 * pi[2]

Epay_b
