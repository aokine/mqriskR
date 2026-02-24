## Listing 5.1: Computing $S_0(t)$, $F_0(t)$, $f_0(t)$, $\lambda_0(t)$, and $\Lambda_0(t)$
library(mqriskR)

# --- Uniform (de Moivre) model ---
omega <- 100
t <- c(0, 10, 25, 50, 99, 100)

S0(t, model = "uniform", omega = omega)
F0(t, model = "uniform", omega = omega)
f0(t, model = "uniform", omega = omega)
hazard0(t, model = "uniform", omega = omega)
cumhaz0(t, model = "uniform", omega = omega)

# --- Exponential (constant force) model ---
lambda <- 0.02
t2 <- c(0, 5, 10, 25, 50)

S0(t2, model = "exponential", lambda = lambda)
F0(t2, model = "exponential", lambda = lambda)
f0(t2, model = "exponential", lambda = lambda)
hazard0(t2, model = "exponential", lambda = lambda)
cumhaz0(t2, model = "exponential", lambda = lambda)



## Listing 5.2: Computing ${}_t p_x$, ${}_t q_x$, and $f_x(t)$

library(mqriskR)

omega <- 120
x <- 30
t <- c(0, 1, 5, 10, 20)

# Conditional survival and failure probabilities
tpx(t, x, model = "uniform", omega = omega)  # {}_t p_x
tqx(t, x, model = "uniform", omega = omega)  # {}_t q_x = 1 - {}_t p_x

# Conditional density f_x(t)
fx(t, x, model = "uniform", omega = omega)

# Verify the identity {}_t p_x = S0(x+t)/S0(x)
lhs <- tpx(t, x, model = "uniform", omega = omega)
rhs <- S0(x + t, model = "uniform", omega = omega) / S0(x, model = "uniform", omega = omega)
max(abs(lhs - rhs))



## Listing 5.3: Computing $\overset{\circ}{e}_x$ and $e_x$
library(mqriskR)

# Uniform model: omega is the limiting age
omega <- 120
x <- c(0, 30, 60, 100)

# Complete expectation of life (continuous-time)
ex_complete(x, model = "uniform", omega = omega)

# Curtate expectation of life (integer part of future lifetime)
ex_curtate(x, model = "uniform", omega = omega)

# Exponential model: lambda is the constant force
lambda <- 0.02
ex_complete(x, model = "exponential", lambda = lambda)
ex_curtate(x, model = "exponential", lambda = lambda)

