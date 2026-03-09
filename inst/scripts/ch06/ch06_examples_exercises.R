# ------------------------------------------------------------
# Example 6.6
# l_x = 10000 / (x + 1)^3
# Verify E[T_x] and Var(T_x | T0 > x)
# ------------------------------------------------------------

lx_fun <- function(x) 10000 / (x + 1)^3
tpx_fun <- function(t, x) ((x + 1)^3) / (x + t + 1)^3

check_ex6_6 <- function(x = 2, tmax = 50) {
  EX_num <- integrate(function(t) tpx_fun(t, x), 0, Inf)$value
  EX2_num <- 2 * integrate(function(t) t * tpx_fun(t, x), 0, Inf)$value
  Var_num <- EX2_num - EX_num^2

  EX_theory <- (x + 1) / 2
  Var_theory <- 3 * (x + 1)^2 / 4

  list(
    x = x,
    EX_numeric = EX_num,
    EX_theory = EX_theory,
    Var_numeric = Var_num,
    Var_theory = Var_theory
  )
}

print(check_ex6_6(x = 0))
print(check_ex6_6(x = 2))
print(check_ex6_6(x = 5))

# Plot conditional survival function

curve(tpx_fun(x, 2), from = 0, to = 20,
      xlab = "t",
      ylab = expression(p[2](t)),
      main = "Example 6.6: Conditional survival function for x = 2",
      lwd = 2)
grid()


# ------------------------------------------------------------
# Example 6.10
# Compare UDD, constant force, and Balducci values of mu_{x+t}
# when l_x = 1000 and l_{x+1} = 900
# ------------------------------------------------------------

check_ex6_10 <- function(lx = 1000, lx1 = 900, t0 = 0.25) {
  qx <- (lx - lx1) / lx
  px <- lx1 / lx

  mu_udd <- qx / (1 - t0 * qx)
  mu_cf <- -log(px)
  mu_bal <- qx / (1 - (1 - t0) * qx)

  list(
    qx = qx,
    px = px,
    mu_udd = mu_udd,
    mu_cf = mu_cf,
    mu_balducci = mu_bal
  )
}

print(check_ex6_10())

# Plot mu_{x+t} across t in [0,1]
qx <- 0.10
px <- 0.90
t_grid <- seq(0, 1, by = 0.001)

mu_udd <- qx / (1 - t_grid * qx)
mu_cf <- rep(-log(px), length(t_grid))
mu_bal <- qx / (1 - (1 - t_grid) * qx)

matplot(t_grid, cbind(mu_udd, mu_cf, mu_bal), type = "l", lwd = 2, lty = c(1:3), col = c(1:3),
        xlab = "t", ylab = expression(mu[x+t]),
        main = expression("Example 6.10: Comparison of " * mu[x+t] * " under three assumptions"))
legend("topright",
       legend = c("UDD", "Constant force", "Balducci"),
       lty = c(1:3),, col = c(1:3), lwd = 2, bty = "n")
grid()




# ------------------------------------------------------------
# Example 6.13
# S0(t) = (9000 - 10 t - t^2)/9000, 0 <= t <= 90
# ------------------------------------------------------------

S0_ex613 <- function(t) (9000 - 10 * t - t^2) / 9000
mu_exact_ex613 <- function(t) (10 + 2 * t) / (9000 - 10 * t - t^2)

check_ex6_13 <- function() {
  l60 <- round(100000 * S0_ex613(60))
  l61 <- round(100000 * S0_ex613(61))

  q60 <- (l60 - l61) / l60
  p60 <- 1 - q60

  # part (b)
  tp_udd <- 1 - 0.25 * q60
  tp_cf  <- p60^0.25
  tp_bal <- (1 - q60) / (1 - 0.75 * q60)

  # exact and approximate mu at 60.75
  mu_exact <- mu_exact_ex613(60.75)
  mu_udd <- q60 / (1 - 0.75 * q60)
  mu_cf  <- -log(p60)
  mu_bal <- q60 / (1 - 0.25 * q60)

  list(
    l60 = l60,
    l61 = l61,
    q60 = q60,
    p60 = p60,
    tpx_udd = tp_udd,
    tpx_cf = tp_cf,
    tpx_bal = tp_bal,
    mu_exact = mu_exact,
    mu_udd = mu_udd,
    mu_cf = mu_cf,
    mu_bal = mu_bal
  )
}

print(check_ex6_13())


# ------------------------------------------------------------
# Plot exact and approximate mu_{x+t} for 60 <= x+t <= 61
# ------------------------------------------------------------

age_grid <- seq(60, 61, by = 0.001)
t_grid <- age_grid - 60

mu_exact_vals <- mu_exact_ex613(age_grid)

out <- check_ex6_13()

qx <- out$q60
px <- out$p60

# approximations as functions of t
mu_udd <- qx / (1 - t_grid * qx)
mu_cf <- rep(-log(px), length(t_grid))
mu_bal <- qx / (1 - (1 - t_grid) * qx)

matplot(age_grid,
        cbind(mu_exact_vals, mu_udd, mu_cf, mu_bal),
        type = "l",
        lwd = 2,
        lty = c(1,2,3,4),
        col = 1,
        xlab = "Age",
        ylab = expression(mu[x]),
        main = expression("Example 6.13: Exact and approximate " * mu[x+t]))

legend("topleft",
       legend = c("Exact", "UDD", "Constant force", "Balducci"),
       lty = c(1,2,3,4),
       lwd = 2,
       bty = "n")

grid()




# ------------------------------------------------------------
# Example 6.15
# Makeham: A = 0.002, B = 1e-5, c = 1.10
# ------------------------------------------------------------

check_ex6_15 <- function(A = 0.002, B = 1e-5, c = 1.10) {
  mu_exact <- hazard0(40, model = "makeham", A = A, B = B, c = c)

  S39 <- S0(39, model = "makeham", A = A, B = B, c = c)
  S40 <- S0(40, model = "makeham", A = A, B = B, c = c)
  S41 <- S0(41, model = "makeham", A = A, B = B, c = c)

  p39 <- S40 / S39
  p40 <- S41 / S40

  q39 <- 1 - p39
  q40 <- 1 - p40

  # UDD approximations
  mu_left  <- q40
  mu_right <- q39 / (1 - q39)

  # midpoint rule approximation
  mu_mid <- -0.5 * (log(p39) + log(p40))

  data.frame(
    method = c("Exact", "UDD left endpoint", "UDD right endpoint", "Midpoint"),
    value = c(mu_exact, mu_left, mu_right, mu_mid),
    error = c(0, mu_left - mu_exact, mu_right - mu_exact, mu_mid - mu_exact)
  )
}

print(check_ex6_15())

# Plot exact mu_x near age 40 with approximation lines
A <- 0.002; B <- 1e-5; c <- 1.10
age_grid <- seq(39.5, 40.5, by = 0.001)
mu_vals <- hazard0(age_grid, model = "makeham", A = A, B = B, c = c)
out <- check_ex6_15(A, B, c)

plot(age_grid, mu_vals, type = "l", lwd = 2,
     xlab = "Age", ylab = expression(mu[x]),
     main = expression("Example 6.15: Exact " * mu[x] * " near age 40 and approximations to " * mu[40]))

abline(h = out$value[2], lty = 2)
abline(h = out$value[3], lty = 3)
abline(h = out$value[4], lty = 4)
abline(v = 40, lty = 5)

points(40, out$value[1], pch = 19)

legend("topleft",
       legend = c("Exact curve", "UDD left endpoint", "UDD right endpoint", "Midpoint", "Exact value at age 40"),
       lty = c(1, 2, 3, 4, NA),
       lwd = c(2, 2, 2, 2, NA),
       pch = c(NA, NA, NA, NA, 19),
       bty = "n")
grid()


out <- check_ex6_15()

barplot(out$error[-1],
        names.arg = c("UDD left", "UDD right", "Midpoint"),
        ylab = "Approximation error",
        main = "Example 6.15: Errors in approximating mu_40")
abline(h = 0, lty = 2)
grid()



# ------------------------------------------------------------
# Exercise 6-3
# Y = min(K0, 4)
# ------------------------------------------------------------

check_ex6_3 <- function() {
  qx <- c(0.10, 0.20, 0.30, 0.40, 0.50)
  px <- 1 - qx

  pY <- c(
    qx[1],
    px[1] * qx[2],
    px[1] * px[2] * qx[3],
    px[1] * px[2] * px[3] * qx[4],
    px[1] * px[2] * px[3] * px[4]
  )

  y <- 0:4
  EY <- sum(y * pY)
  EY2 <- sum((y^2) * pY)
  VarY <- EY2 - EY^2

  list(
    distribution = data.frame(Y = y, prob = pY),
    prob_sum = sum(pY),
    EY = EY,
    VarY = VarY
  )
}

out <- check_ex6_3()
print(out$distribution)
print(out$prob_sum)
print(out$EY)
print(out$VarY)

barplot(out$distribution$prob,
        names.arg = out$distribution$Y,
        xlab = "Y",
        ylab = "Probability",
        main = "Exercise 6-3: Distribution of Y = min(K0, 4)")
grid()



# ------------------------------------------------------------
# Exercise 6-19
# Verify pi_4 from Table 6.1a
# ------------------------------------------------------------

check_ex6_19 <- function(make_plot = TRUE) {
  lx_vals <- c(100000, 97408, 97259, 97160, 97082)
  names(lx_vals) <- 0:4

  px <- lx_vals[-1] / lx_vals[-length(lx_vals)]
  qx <- 1 - px

  P_list <- lapply(seq_along(px), function(i) {
    matrix(c(px[i], qx[i],
             0,     1),
           nrow = 2, byrow = TRUE)
  })

  pi0 <- c(1, 0)
  pi1 <- drop(pi0 %*% P_list[[1]])
  pi2 <- drop(pi1 %*% P_list[[2]])
  pi3 <- drop(pi2 %*% P_list[[3]])
  pi4 <- drop(pi3 %*% P_list[[4]])

  out <- list(
    px = px,
    qx = qx,
    pi1 = pi1,
    pi2 = pi2,
    pi3 = pi3,
    pi4 = pi4
  )

  if (make_plot) {
    probs_alive_dead <- rbind(
      c(1, 0),
      pi1, pi2, pi3, pi4
    )
    matplot(0:4, probs_alive_dead, type = "l", lwd = 2, lty = 1:2,col = 1:2,
            xlab = "Time",
            ylab = "State probability",
            main = "Exercise 6-19: State probabilities over time")
    legend("right", legend = c("Alive (State 0)", "Dead (State 1)"),
           lty = 1:2,col = 1:2, lwd = 2, bty = "n")
    grid()
  }

  out
}

out <- check_ex6_19()
print(out$pi4)


# ------------------------------------------------------------
# Exercise 6-27
# l_x = 1000 * (100 - x)^(1/2)
# ------------------------------------------------------------

check_ex6_27 <- function(make_plot = TRUE) {
  lx_fun <- function(x) 1000 * sqrt(100 - x)
  mu_exact_fun <- function(x) 0.5 / (100 - x)

  l36 <- lx_fun(36)
  l37 <- round(lx_fun(37))   # as instructed
  p36 <- l37 / l36
  q36 <- 1 - p36

  t0 <- 0.25
  mu_exact <- mu_exact_fun(36.25)
  mu_udd <- q36 / (1 - t0 * q36)
  mu_cf <- -log(p36)
  mu_bal <- q36 / (1 - (1 - t0) * q36)

  out <- data.frame(
    method = c("Exact", "UDD", "Constant force", "Balducci"),
    value = c(mu_exact, mu_udd, mu_cf, mu_bal)
  )

  if (make_plot) {
    age_grid <- seq(36, 37, by = 0.001)
    t_grid <- age_grid - 36
    mu_exact_vals <- mu_exact_fun(age_grid)
    mu_udd_vals <- q36 / (1 - t_grid * q36)
    mu_cf_vals <- rep(-log(p36), length(age_grid))
    mu_bal_vals <- q36 / (1 - (1 - t_grid) * q36)

    matplot(age_grid,
            cbind(mu_exact_vals, mu_udd_vals, mu_cf_vals, mu_bal_vals),
            type = "l", lwd = 2, lty = c(1,2,3,4), col = 1,
            xlab = "Age", ylab = expression(mu[x]),
            main = expression("Exercise 6-27: Exact and approximate " * mu[x] * " on (36,37]"))
    legend("topleft",
           legend = c("Exact", "UDD", "Constant force", "Balducci"),
           lty = c(1,2,3,4), lwd = 2, bty = "n")
    grid()
  }

  out
}

out <- check_ex6_27()
print(out)



# ------------------------------------------------------------
# Exercise 6-37
# Select table with 2-year select period, UDD
# ------------------------------------------------------------

check_ex6_37 <- function() {
  l_60   <- 80625
  l_60_1 <- 79954
  l_62   <- 78839

  q_60_sel <- 1 - l_60_1 / l_60
  q_60_1_sel <- 1 - l_62 / l_60_1

  # from age [60]+0.60 to [60]+1.00 : 0.40 years left in current year
  q_part1 <- (0.40 * q_60_sel) / (1 - 0.60 * q_60_sel)

  # then 0.50 years in next year under UDD
  q_part2 <- 0.50 * q_60_1_sel

  total_q <- q_part1 + (1 - q_part1) * q_part2

  list(
    q_sel_60 = q_60_sel,
    q_sel_60_plus_1 = q_60_1_sel,
    q_part1 = q_part1,
    q_part2 = q_part2,
    total_q = total_q
  )
}

out <- check_ex6_37()
print(out)
