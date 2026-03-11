## R check for Example 7.1
k <- 0:4
pk <- c(0.20, 0.30, 0.20, 0.15, 0.15)
i <- 0.01

z <- 10 * (1 + i)^(-(k + 1))

EZ  <- sum(z * pk)
EZ2 <- sum(z^2 * pk)
VarZ <- EZ2 - EZ^2

cdf <- cumsum(pk[order(z)])
z_sorted <- sort(z)
median_z <- z_sorted[min(which(cdf >= 0.50))]

round(c(EZ = EZ, VarZ = VarZ, median = median_z), 5)

plot(z_sorted, cdf, type = "s",
     xlab = "z", ylab = "F_Z(z)",
     main = "Example 7.1: CDF of the present value random variable")
grid()




## R check for Example 7.5
# Example 7.1 distribution reused
k <- 0:4
pk <- c(0.20, 0.30, 0.20, 0.15, 0.15)
i <- 0.01
v <- 1 / (1 + i)

# 3-year term insurance part
Z_term <- ifelse(k < 3, v^(k + 1), 0)

# 3-year pure endowment part
Z_pe <- ifelse(k >= 3, v^3, 0)

# 3-year endowment insurance
Z_endow <- Z_term + Z_pe

A_term <- sum(Z_term * pk)
A_pe   <- sum(Z_pe * pk)
A_end  <- sum(Z_endow * pk)

A2_term <- sum(Z_term^2 * pk)
A2_pe   <- sum(Z_pe^2 * pk)
A2_end  <- sum(Z_endow^2 * pk)

cov_tp <- sum(Z_term * Z_pe * pk) - A_term * A_pe

round(c(A_term = A_term, A_pe = A_pe, A_end = A_end,
        A2_term = A2_term, A2_pe = A2_pe, A2_end = A2_end,
        covariance = cov_tp), 5)



## R check for Example 7.7
library(mqriskR)

i <- 0.05
delta <- interest_convert(i = i)$delta
lambda <- 0.25

Abar <- Abarx(x = 0, i = i, model = "exponential", lambda = lambda)
VarAbar <- var_Abarx(x = 0, i = i, model = "exponential", lambda = lambda)

# 90th percentile of Z = exp(-delta T)
p90 <- 0.90^(delta / lambda)

round(c(Abar = Abar, Var = VarAbar, p90 = p90), 5)

t_grid <- seq(0, 20, by = 0.01)
z_grid <- exp(-delta * t_grid)

plot(t_grid, z_grid, type = "l", lwd = 2,
     xlab = "t", ylab = "z = exp(-delta t)",
     main = "Example 7.7: Present value as a function of failure time")
grid()




## R check for Example 7.9
k <- 0:4
pk <- c(0.20, 0.30, 0.20, 0.15, 0.15)
i <- 0.01
v <- 1 / (1 + i)

IAx <- sum((k + 1) * v^(k + 1) * pk)

n <- 3
DAxn1 <- sum((n - k[1:n]) * v^(k[1:n] + 1) * pk[1:n])

round(c(IAx = IAx, DAxn1 = DAxn1), 5)




# Exercise 7-3
# Cov(Z_{x:n}^1, {}_{n|}Z_x) for the Example 7.1 model

k <- 0:4
pk <- c(0.20, 0.30, 0.20, 0.15, 0.15)
i <- 0.01
v <- 1 / (1 + i)
n <- 3

U <- ifelse(k < n, v^(k + 1), 0)   # term insurance part
V <- ifelse(k >= n, v^(k + 1), 0)  # deferred insurance part

EU <- sum(U * pk)
EV <- sum(V * pk)
EUV <- sum(U * V * pk)
cov_uv <- EUV - EU * EV

round(c(EU = EU, EV = EV, EUV = EUV, covariance = cov_uv), 6)

barplot(rbind(U * pk, V * pk), beside = TRUE,
        names.arg = k,
        legend.text = c("Term part", "Deferred part"),
        xlab = "k", ylab = "Probability-weighted PV",
        main = "Exercise 7-3: Term and deferred components")
grid()



# Exercise 7-6
# Check the recursion A_x = v q_x + v p_x A_{x+1}

library(mqriskR)

x <- 40
i <- 0.05

lhs <- Ax(x, i, model = "uniform", omega = 100)
rhs <- discount(i, 1) * (
  tqx(1, x = x, model = "uniform", omega = 100) +
    tpx(1, x = x, model = "uniform", omega = 100) *
    Ax(x + 1, i, model = "uniform", omega = 100)
)

round(c(lhs = lhs, rhs = rhs, difference = lhs - rhs), 10)

ages <- 40:70
lhs_vec <- sapply(ages, function(a) Ax(a, i, model = "uniform", omega = 100))
rhs_vec <- sapply(ages, function(a) {
  discount(i, 1) * (
    tqx(1, x = a, model = "uniform", omega = 100) +
      tpx(1, x = a, model = "uniform", omega = 100) *
      Ax(a + 1, i, model = "uniform", omega = 100)
  )
})

plot(ages, lhs_vec - rhs_vec, type = "h",
     xlab = "x", ylab = "LHS - RHS",
     main = "Exercise 7-6: Recursion check error")
abline(h = 0, lty = 2)
grid()




# Exercise 7-11
# Check {}_{n|}A_x = {}_nE_x * A_{x+n}

library(mqriskR)

x <- 40
n <- 10
i <- 0.05

lhs <- nAx(x, n, i, model = "uniform", omega = 100)
rhs <- nEx(x, n, i, model = "uniform", omega = 100) *
  Ax(x + n, i, model = "uniform", omega = 100)

round(c(lhs = lhs, rhs = rhs, difference = lhs - rhs), 10)

n_grid <- 1:25
lhs_vec <- sapply(n_grid, function(nn) nAx(x, nn, i, model = "uniform", omega = 100))
rhs_vec <- sapply(n_grid, function(nn) {
  nEx(x, nn, i, model = "uniform", omega = 100) *
    Ax(x + nn, i, model = "uniform", omega = 100)
})

plot(n_grid, lhs_vec, type = "l", lwd = 2,
     xlab = "n", ylab = "Value",
     main = "Exercise 7-11: Deferred insurance identity")
lines(n_grid, rhs_vec, lty = 2, lwd = 2)
legend("topright", legend = c("nAx", "nEx * A_{x+n}"),
       lty = c(1, 2), lwd = 2, bty = "n")
grid()



# Exercise 7-15
# Continuous APV with given density

delta <- 0.10

f <- function(t) ifelse(t > 0 & t <= 100, t / 5000, 0)
integrand <- function(t) 50 * exp(-delta * t) * f(t)

apv_num <- integrate(integrand, lower = 0, upper = 100)$value
apv_formula <- 1 - 11 * exp(-10)

round(c(numerical = apv_num, closed_form = apv_formula,
        difference = apv_num - apv_formula), 8)

t_grid <- seq(0, 100, by = 0.1)
plot(t_grid, integrand(t_grid), type = "l", lwd = 2,
     xlab = "t", ylab = "Discounted density contribution",
     main = "Exercise 7-15: Integrand for the APV")
grid()



# Exercise 7-21
# Deferred continuous insurance under exponential lifetime

delta <- 0.05
lambda <- 0.25
n <- 10

integrand <- function(t) exp(-delta * t) * exp(-lambda * t) * lambda

lhs_num <- integrate(integrand, lower = n, upper = Inf)$value
rhs_formula <- (lambda / (delta + lambda)) * exp(-(delta + lambda) * n)

round(c(numerical = lhs_num, closed_form = rhs_formula,
        difference = lhs_num - rhs_formula), 8)

t_grid <- seq(n, 40, by = 0.05)
plot(t_grid, integrand(t_grid), type = "l", lwd = 2,
     xlab = "t", ylab = "Integrand",
     main = "Exercise 7-21: Deferred continuous insurance integrand")
grid()



# Exercise 7-25
# Check (IA)_x = A_x + {}_1E_x (IA)_{x+1}

IA35 <- 3.711
A35  <- 0.1300
p35  <- 0.9964
v    <- 0.9434

E1_35 <- v * p35
IA36 <- (IA35 - A35) / E1_35

round(c(E1_35 = E1_35, IA36 = IA36), 6)

lhs <- IA35
rhs <- A35 + E1_35 * IA36

round(c(lhs = lhs, rhs = rhs, difference = lhs - rhs), 10)

barplot(c(A35, E1_35 * IA36),
        names.arg = c("A35", "E1_35 * IA36"),
        ylab = "Contribution",
        main = "Exercise 7-25: Decomposition of (IA)35")
grid()



# Exercise 7-28
# Recover A_35 from UDD information

i <- 0.05
delta <- log(1 + i)
q35 <- 0.01
Abar36 <- 0.185
v <- 1 / (1 + i)
p35 <- 1 - q35

Abar35_1 <- (i / delta) * v * q35
Abar35 <- Abar35_1 + v * p35 * Abar36
A35 <- (delta / i) * Abar35

round(c(Abar35_1 = Abar35_1, Abar35 = Abar35, A35 = A35), 6)

barplot(c(Abar35_1, v * p35 * Abar36),
        names.arg = c("1-year term part", "survival continuation"),
        ylab = "Contribution",
        main = "Exercise 7-28: Components of Abar35")
grid()





# Exercise 7-29
# UDD approximation for 2-year term insurance

i <- 0.10
delta <- log(1 + i)
q_x <- 0.05
q_x1 <- 0.08
p_x <- 1 - q_x
v <- 1 / (1 + i)

A_disc <- v * q_x + v^2 * p_x * q_x1
A_cont_udd <- (i / delta) * A_disc

round(c(A_disc = A_disc, Abar_udd = A_cont_udd), 6)

barplot(c(v * q_x, v^2 * p_x * q_x1),
        names.arg = c("year 1", "year 2"),
        ylab = "Contribution",
        main = "Exercise 7-29: Annual term-insurance contributions")
grid()





# Exercise 7-31
# Visualize piecewise-constant vs fully continuous increasing benefit

t_grid <- seq(0, 5, by = 0.001)
b_piece <- floor(t_grid + 1)
b_cont <- t_grid

plot(t_grid, b_piece, type = "s", lwd = 2,
     xlab = "t", ylab = "Benefit amount",
     main = "Exercise 7-31: Piecewise vs continuous increasing benefit")
lines(t_grid, b_cont, lwd = 2, lty = 2)
legend("topleft",
       legend = c("piecewise increasing", "continuous increasing"),
       lty = c(1, 2), lwd = 2, bty = "n")
grid()

