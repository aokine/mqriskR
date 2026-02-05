# Chapter 3: Payment expectations for a DTMC example (Exercise 3-4 style)

P0 <- matrix(c(
  0.60, 0.30, 0.10,
  0.00, 0.00, 1.00,
  0.00, 0.00, 1.00
), nrow = 3, byrow = TRUE)

P1 <- P0

Pk <- matrix(c(
  0.00, 0.30, 0.70,
  0.00, 0.00, 1.00,
  0.00, 0.00, 1.00
), nrow = 3, byrow = TRUE)

# Sequence: P^(0)=P0, P^(1)=P1, P^(k)=Pk for k>=2
Ps <- function(k) if (k == 0) P0 else if (k == 1) P1 else Pk

# Start in State 0 at time 0
pi <- c(1, 0, 0)

# (a) payment 1 at times t=0,1,2,... if in state 0 or 1
# We'll compute expected payments up to some horizon N and show convergence
N <- 15
pay_a <- numeric(N + 1)
pay_b <- numeric(N + 1)

for (t in 0:N) {
  # payment at time t depends on state at time t
  pay_a[t + 1] <- 1 * (pi[1] + pi[2])     # in state 0 or 1
  pay_b[t + 1] <- if (t >= 1) 4 * pi[2] else 0  # in state 1, for t>=1

  # advance state if not at horizon
  if (t < N) {
    pi <- drop(pi %*% Ps(t))
  }
}

list(
  expected_payments_a_partial_sum = sum(pay_a),
  expected_payments_b_partial_sum = sum(pay_b),
  pay_a = pay_a,
  pay_b = pay_b
)
