## Listing 3.1: Homogeneous DTMC: state vectors and multi-step probabilities via matrix powers
P <- matrix(c(
  0.80, 0.20, 0.00,
  0.30, 0.60, 0.10,
  0.00, 0.00, 1.00
), nrow = 3, byrow = TRUE)

# Known to be in State 1 at time n:
pi_n <- c(0, 1, 0)

# One-step state vector:
pi_np1 <- drop(pi_n %*% P)

# Two-step state vector using P^2:
P2 <- P %*% P
pi_np2 <- drop(pi_n %*% P2)

pi_np1
pi_np2
P2


## Listing 3.2: Non-homogeneous DTMC: sequential transition matrices
P0 <- matrix(c(
  0.60, 0.40,
  0.70, 0.30
), nrow = 2, byrow = TRUE)

P1 <- matrix(c(
  0.50, 0.50,
  0.80, 0.20
), nrow = 2, byrow = TRUE)

# Begin in State 1 at time 0:
pi0 <- c(0, 1)

pi1 <- drop(pi0 %*% P0)
pi2 <- drop(pi1 %*% P1)

pi1
pi2
pi2[1]  # P(State 0 at time 2 | State 1 at time 0)


## Listing 3.3 CTMC transition probabilities via the matrix exponential expm(Q*t)
# Rates from State 0 to States 1,2,3
mu01 <- 0.30
mu02 <- 0.50
mu03 <- 0.70

Q <- matrix(0, nrow = 4, ncol = 4)
Q[1, 2] <- mu01
Q[1, 3] <- mu02
Q[1, 4] <- mu03
Q[1, 1] <- -(mu01 + mu02 + mu03)

t <- 1
Pt <- expm::expm(Q * t)

Pt
Pt[1, 3]  # P(X(1)=2 | X(0)=0)



## Listing 3.4: Expected payments in a non-homogeneous DTMC (state-probability recursion)
# State order is: 0, 1, 2  (so pi[1]=State0, pi[2]=State1, pi[3]=State2)

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

P_step <- function(k) if (k == 0) P0 else if (k == 1) P1 else Pk

state0 <- 1
state1 <- 2
state2 <- 3

pi <- c(1, 0, 0)  # start in State 0 at time 0
N <- 15           # compute partial sums through time N

pay_a <- numeric(N + 1)  # payments for part (a)
pay_b <- numeric(N + 1)  # payments for part (b)

for (t in 0:N) {
  # (a) Pay 1 at time t if in State 0 or State 1
  pay_a[t + 1] <- 1 * (pi[state0] + pi[state1])

  # (b) Pay 4 at time t>=1 if in State 1
  pay_b[t + 1] <- if (t >= 1) 4 * pi[state1] else 0

  # update pi(t) -> pi(t+1) using the matrix for interval (t, t+1]
  if (t < N) pi <- drop(pi %*% P_step(t))
}

sum(pay_a)  # partial expected value of payments in (a) through time N
sum(pay_b)  # partial expected value of payments in (b) through time N


