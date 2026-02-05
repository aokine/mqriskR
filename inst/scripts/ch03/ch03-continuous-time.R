# Chapter 3 — Continuous-time Markov chains with expm

# Install expm if needed:
# install.packages("expm")

library(expm)

# Example structure: 4-state multiple decrement from State 0 -> {1,2,3}
# with constant forces mu01=.30, mu02=.50, mu03=.70 for t>=0.
# Generator Q has:
# q_0j = mu_0j, q_00 = -sum(mu_0j)
# absorbing states 1,2,3: q_jj=0, others 0

mu01 <- 0.30
mu02 <- 0.50
mu03 <- 0.70
mu0  <- mu01 + mu02 + mu03

Q <- matrix(0, nrow = 4, ncol = 4)
Q[1, 2] <- mu01
Q[1, 3] <- mu02
Q[1, 4] <- mu03
Q[1, 1] <- -mu0
# states 1,2,3 absorbing => rows 2:4 already zero

Q

t <- 1
P_t <- expm(t * Q)
P_t

# Probability P[X(1)=2 | X(0)=0] corresponds to State 2 (text) = index 3 in R
# because R indices: 1->State0, 2->State1, 3->State2, 4->State3
P_t[1, 3]
