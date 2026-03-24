##Listing 14.1: Cause-specific multiple-decrement APV
library(mqriskR)

q1 <- c(.02, .02, .02, .02, .02)   # default
q2 <- c(.03, .04, .05, .06, .00)   # call
q3 <- c(.00, .00, .00, .00, .98)   # maturity

qtau <- q1 + q2 + q3

ptau <- numeric(length(qtau))
ptau[1] <- 1
for (k in 2:length(qtau)) {
  ptau[k] <- prod(1 - qtau[1:(k - 1)])
}

Axj_md(qj = q1, ptau = ptau, i = 0.06, benefit = 1000)


##Listing 14.2: Projected asset share path
library(mqriskR)

r <- c(.05, .05, .05, .05, .05)
e <- c(30, 30, 30, 30, 30)

b1 <- c(1000, 1000, 1000, 1000, 1000)   # death benefit
b2 <- c(50, 100, 300, 600, 0)           # withdrawal benefit
b3 <- c(0, 0, 0, 0, 1000)               # endowment benefit

q1 <- c(.02, .03, .04, .05, .06)
q2 <- c(.30, .20, .20, .10, .00)
p_tau <- c(.68, .77, .76, .85, .94)

AS_path(AS0 = 50, G = 200, r = r, e = e, b1 = b1, b2 = b2, q1 = q1, q2 = q2, p_tau = p_tau, i = 0.06, b3 = b3)



##Listing 14.3: Euler approximation of state probabilities
library(mqriskR)

mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

out <- tp00_tp01_euler(
  h = 0.10,
  n = 2.00,
  mu01 = mu01,
  mu02 = mu02,
  mu10 = mu10,
  mu12 = mu12
)

out
subset(out, t %in% c(1.0, 2.0))



##Listing 4.4: Three-step transition probabilities in the CCRC model
library(mqriskR)

P <- matrix(
  c(0.94, 0.03, 0.02, 0.01,
    0.50, 0.30, 0.18, 0.02,
    0.00, 0.00, 0.93, 0.07,
    0.00, 0.00, 0.00, 1.00),
  nrow = 4, byrow = TRUE
)

c(
  p00_3 = markov_nstep_prob(P, n = 3, i = 1, j = 1),
  p02_3 = markov_nstep_prob(P, n = 3, i = 1, j = 3)
)


## Listing 14.5: Backward reserve path for the recovery model
library(mqriskR)

mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

out <- thiele_path_01(h = 0.10, n = 2.00, delta = 0.04, Pbar = 446.95, B = 1000, R = 1000, mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12)

out
subset(out, t %in% c(1.90, 1.00, 0.00))




## Listing 14.6: Continuous premium approximation in the disability model with recovery
library(mqriskR)

mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20

ex14_9 <- tp00_tp01_euler(
  h = 0.10,
  n = 2.00,
  mu01 = mu01,
  mu02 = mu02,
  mu10 = mu10,
  mu12 = mu12
)

Pbar_trapz_ms(
  t = ex14_9$t,
  tp00 = ex14_9$tp00,
  tp01 = ex14_9$tp01,
  delta = 0.04,
  mu02 = mu02,
  mu12 = mu12,
  B02 = 1000,
  B12 = 1000,
  R = 1000
)

