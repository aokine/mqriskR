
## In-Text Examples
## Example 5.2 (PDF and HRF)
# Example 5.2 verification:
# F0(t) = 1 - 0.10 * sqrt(100 - t),  0 <= t <= 100

F0 <- function(t) {
  out <- 1 - 0.10 * sqrt(pmax(100 - t, 0))
  # enforce support
  out[t < 0] <- 0
  out[t > 100] <- 1
  out
}

S0 <- function(t) 1 - F0(t)

# Closed-form results from the solution:
f0_closed <- function(t) ifelse(t >= 0 & t < 100, 0.05 * (100 - t)^(-1/2), 0)
haz_closed <- function(t) ifelse(t >= 0 & t < 100, 0.50 * (100 - t)^(-1), Inf)

# Numerical derivative check for f0(t) = d/dt F0(t)
f0_num <- function(t, h = 1e-6) (F0(t + h) - F0(t - h)) / (2*h)

t_test <- c(0, 1, 10, 30, 75, 99)

cbind(
  t = t_test,
  f0_closed = f0_closed(t_test),
  f0_num = f0_num(t_test),
  haz_closed = haz_closed(t_test),
  haz_from_def = f0_closed(t_test) / S0(t_test)
)

# max absolute differences (ignoring endpoints issues)
max(abs(f0_closed(t_test) - f0_num(t_test)), na.rm = TRUE)
max(abs(haz_closed(t_test) - f0_closed(t_test) / S0(t_test)), na.rm = TRUE)


## Example 5.3 (mean and median)
# Example 5.3 verification for the same survival model as Example 5.2
S0 <- function(t) {
  out <- 0.10 * sqrt(pmax(100 - t, 0))
  out[t < 0] <- 1
  out[t > 100] <- 0
  out
}

# Verify E[T0] = integral_0^infty S0(t) dt = integral_0^100 S0(t) dt
E_num <- integrate(function(t) S0(t), lower = 0, upper = 100)$value
E_closed <- 200/3
c(E_num = E_num, E_closed = E_closed, abs_diff = abs(E_num - E_closed))

# Verify median y solves S0(y) = 0.5 -> y = 75
median_closed <- 75
S0(median_closed)

# Numerical solve for median as a check
median_num <- uniroot(function(y) S0(y) - 0.5, interval = c(0, 100))$root
c(median_num = median_num, median_closed = median_closed, abs_diff = abs(median_num - median_closed))



## Example 5.4 (conditional survival, density, and expectation)
# Example 5.4 verification for the model from Example 5.2
S0 <- function(t) {
  out <- 0.10 * sqrt(pmax(100 - t, 0))
  out[t < 0] <- 1
  out[t > 100] <- 0
  out
}

f0 <- function(t) ifelse(t >= 0 & t < 100, 0.05 * (100 - t)^(-1/2), 0)

# (a) _20 p_36 = S0(56)/S0(36)
p20_36_num <- S0(56) / S0(36)
p20_36_closed <- 0.82916  # from the solution (rounded)
c(p20_36_num = p20_36_num, p20_36_closed = p20_36_closed, abs_diff = abs(p20_36_num - p20_36_closed))

# (b) fx(t) = f0(36+t)/S0(36). Compare to closed-form 0.0625 / sqrt(64 - t)
fx_num <- function(t) f0(36 + t) / S0(36)
fx_closed <- function(t) 0.0625 / sqrt(64 - t)

t_test <- c(0, 1, 10, 20, 40, 63)
cbind(t = t_test, fx_num = fx_num(t_test), fx_closed = fx_closed(t_test), abs_diff = abs(fx_num(t_test) - fx_closed(t_test)))

# (c) e36^o = integral_0^infty _t p_36 dt = integral_0^64 S0(36+t)/S0(36) dt
tp36 <- function(t) S0(36 + t) / S0(36)
e36_num <- integrate(tp36, lower = 0, upper = 64)$value
e36_closed <- 128/3
c(e36_num = e36_num, e36_closed = e36_closed, abs_diff = abs(e36_num - e36_closed))



## Example 5.7 (validity checks for candidate survival functions)
# Example 5.7 verification: check S(0)=1, limit->0 (numerically), and non-increasing.

S1 <- function(t) exp(t - 0.70 * (2^t - 1))
S2 <- function(t) (1 + t)^(-2)
S3 <- function(t) exp(-t^2)

grid <- seq(0, 20, by = 0.01)

check_S <- function(Sfun, name) {
  S0_val <- Sfun(0)
  # monotonic: all differences <= 0 (allow tiny numerical slack)
  diffs <- diff(Sfun(grid))
  non_increasing <- all(diffs <= 1e-10)
  # "limit" check via large t
  tail_val <- Sfun(max(grid))
  list(
    name = name,
    S_at_0 = S0_val,
    approx_limit = tail_val,
    non_increasing = non_increasing,
    max_positive_diff = max(diffs)
  )
}

check_S(S1, "S1")
check_S(S2, "S2")
check_S(S3, "S3")

# Optional: visualize the functions
plot(grid, S1(grid), type="l", col="red", ylim=c(0,1.5),
     main="Candidate Survival Functions", ylab="S(t)", xlab="t")
lines(grid, S2(grid), col="blue")
lines(grid, S3(grid), col="darkgreen")
legend("topright", legend=c("S1","S2","S3"),
       col=c("red","blue","darkgreen"), lty=1)

# Optional: check derivative sign near 0 for S1 (as in the written solution)
h <- 1e-6
S1_prime0_num <- (S1(h) - S1(0)) / h
S1_prime0_num



## Example 5.8 (polynomial survival model)}, label={lst:ex5_8_verify}]
# Example 5.8:
# S0(t) = (9000 - 10 t - t^2)/9000 on 0 <= t <= omega

S0 <- function(t) (9000 - 10*t - t^2) / 9000

# (a) Find omega from S0(omega)=0
omega_num <- uniroot(function(u) S0(u), interval = c(0, 200))$root
omega_num

# (b) f0(30) = -S0'(30)
S0_prime <- function(t) -(10 + 2*t) / 9000
f0 <- function(t) -S0_prime(t)
f0(30)

# (c) f20(10) = f0(30)/S0(20)
f20_10 <- f0(30) / S0(20)
f20_10

# (d) lambda0(30) = f0(30)/S0(30)
haz0_30 <- f0(30) / S0(30)
haz0_30

# (e) lambda20(10 | T0>20) = hazard at age 30 (should match part d)
# Build conditional survival: S20(t) = S0(20+t)/S0(20)
S20 <- function(t) S0(20 + t) / S0(20)
# hazard for T20 at duration t: -S20'(t)/S20(t)
S20_prime <- function(t) S0_prime(20 + t) / S0(20)
haz20 <- function(t) -S20_prime(t) / S20(t)

haz20_10 <- haz20(10)
c(haz0_30 = haz0_30, haz20_10 = haz20_10, abs_diff = abs(haz0_30 - haz20_10))


## End of Chapter Exercises

############################################################
## Chapter 5: Optional R checks for selected end-of-chapter exercises
############################################################

## Helper: safe numeric integration wrapper
num_int <- function(f, lower, upper, ...) {
  integrate(f, lower = lower, upper = upper, ..., subdivisions = 2000L, rel.tol = 1e-10)$value
}

## Helper: finite-difference derivative
fd_deriv <- function(f, x, h = 1e-6) (f(x + h) - f(x - h)) / (2 * h)

############################################################
## Exercise 5-1: lambda0(t) = a + b t
############################################################
check_ex5_1 <- function(a = 0.02, b = 0.001) {

  lambda0 <- function(t) a + b * t

  # Survival reconstructed from hazard
  S0 <- function(t) exp(-(a * t + 0.5 * b * t^2))

  f0 <- function(t) S0(t) * lambda0(t)

  # Choose plotting grid
  grid <- seq(0, 150, length.out = 20001)

  # Numerical integration check (truncated)
  mass <- sum(f0(grid)) * (grid[2] - grid[1])

  # Numerical mode
  mode_num <- grid[which.max(f0(grid))]

  # Analytical mode
  mode_exact <- (-a + sqrt(b)) / b

  # Plot
  plot(grid, f0(grid), type = "l", lwd = 2,
       main = "Density f0(t) with Numerical and Analytical Mode",
       xlab = "t", ylab = "f0(t)")

  abline(v = mode_num, col = "red", lty = 2, lwd = 2)
  abline(v = mode_exact, col = "blue", lty = 3, lwd = 2)

  legend("topright",
         legend = c("Numerical Mode", "Analytical Mode"),
         col = c("red", "blue"),
         lty = c(2,3),
         lwd = 2)

  list(
    a = a,
    b = b,
    approx_mass = mass,
    mode_numeric = mode_num,
    mode_analytical = mode_exact
  )
}

print(check_ex5_1())




############################################################
## Exercise 5-2: S0(t) = a t^2 + b with E[T0]=60 -> omega=90, median
############################################################
check_ex5_2 <- function(omega = 90) {

  S0 <- function(t) 1 - (t^2) / (omega^2)

  # E[T0] = integral_0^omega S0(t) dt
  ET0 <- num_int(S0, 0, omega)

  # Median m: S0(m) = 0.5
  med <- uniroot(function(m) S0(m) - 0.5, lower = 0, upper = omega)$root

  list(
    omega = omega,
    ET0_numeric = ET0,
    median_numeric = med
  )
}

print(check_ex5_2())



############################################################
## Exercise 5-3: lambda0(t) = exp(-r t) is not valid
## because S0(t) -> exp(-1/r) > 0
############################################################

check_ex5_3 <- function(r = 0.5, t_max = 20) {

  S0 <- function(t) exp((exp(-r * t) - 1) / r)

  tgrid <- seq(0, t_max, length.out = 2000)
  Svals <- S0(tgrid)

  limit_theory <- exp(-1 / r)

  # Plot
  plot(tgrid, Svals, type = "l", lwd = 2,
       main = expression(paste("Survival Function ", S[0](t),
                               " with ", lambda[0](t) == e^{-r*t})),
       xlab = "t", ylab = expression(S[0](t)))

  abline(h = limit_theory, lty = 2)
  text(t_max * 0.6, limit_theory,
       labels = paste0("Limit = ", round(limit_theory, 4)),
       pos = 3)

  list(
    r = r,
    theoretical_limit = limit_theory,
    S_at_large_t = S0(t_max)
  )
}

print(check_ex5_3())



############################################################
## Exercise 5-4: modified de Moivre with gamma=4, e0 = omega/5
############################################################
check_ex5_4 <- function(omega = 80, gamma = 4) {

  S0 <- function(t) (1 - t / omega)^gamma

  e0_num <- num_int(S0, 0, omega)
  e0_theory <- omega / (gamma + 1)

  list(
    omega = omega, gamma = gamma,
    e0_numeric = e0_num,
    e0_theory = e0_theory,
    diff = e0_num - e0_theory
  )
}

print(check_ex5_4())




############################################################
## Exercise 5-7: exponential X1,X2; solve for lambdas from S_Y(2), S_Z(2)
############################################################
check_ex5_7 <- function(SY2 = 0.24, SZ2 = 0.86) {

  # From theory: SY(2) = exp(-2(l1+l2)) => l1+l2 = -(1/2) log(SY2)
  sum_l <- -(1/2) * log(SY2)

  # SZ(2) = exp(-2l1) + exp(-2l2) - exp(-2(l1+l2))
  # Let a=exp(-2l1), b=exp(-2l2), then:
  # ab = SY2,  a+b-ab = SZ2 => a+b = SZ2 + SY2
  s <- SZ2 + SY2

  # Solve a( s - a ) = SY2 => a^2 - s a + SY2 = 0
  disc <- s^2 - 4 * SY2
  a1 <- (s + sqrt(disc)) / 2
  a2 <- (s - sqrt(disc)) / 2

  # Convert back to lambdas
  lam_from_a <- function(a) -(1/2) * log(a)
  cand <- data.frame(
    a = c(a1, a2),
    lambda = c(lam_from_a(a1), lam_from_a(a2))
  )

  # identify which corresponds to lambda1 > lambda2:
  # smaller a means larger lambda
  a_small <- min(a1, a2)
  a_large <- max(a1, a2)
  lambda1 <- lam_from_a(a_small)
  lambda2 <- lam_from_a(a_large)

  list(
    sum_lambda = sum_l,
    candidates_for_exp_neg2lambda = cand,
    lambda1 = lambda1,
    lambda2 = lambda2,
    check_order = (lambda1 > lambda2)
  )
}



print(check_ex5_7())



############################################################
## Exercise 5-8: deferred probability m|q0 = S(m) - S(m+1)
## Plot + correct constancy/monotonicity checks (exclude boundary m=omega)
############################################################
check_ex5_8 <- function(omega = 40, lambda = 0.2, make_plot = TRUE, include_endpoint_in_preview = TRUE) {

  # For behavior checks we only use m = 0,...,omega-1 (intervals fully inside support)
  m_check <- 0:(omega - 1)

  # For display/plot we can optionally include m=omega (shows boundary artifact clearly)
  m_grid <- if (include_endpoint_in_preview) 0:omega else m_check

  # (a) Uniform on (0, omega): S(t) = 1 - t/omega on [0,omega], 0 after omega
  S_unif <- function(t) ifelse(t < 0, 1, ifelse(t <= omega, 1 - t / omega, 0))

  dq_unif <- S_unif(m_grid) - S_unif(m_grid + 1)
  dq_unif_check <- S_unif(m_check) - S_unif(m_check + 1)

  # (b) Exponential with rate lambda: S(t) = exp(-lambda t), t>=0
  S_exp <- function(t) exp(-lambda * pmax(t, 0))

  dq_exp <- S_exp(m_grid) - S_exp(m_grid + 1)
  dq_exp_check <- S_exp(m_check) - S_exp(m_check + 1)

  # (c) f0(t)=0.00125 t on [0,omega] => S(t)=1 - 0.000625 t^2 on [0,omega]
  # Clamp t to [0,omega] so S is well-defined at endpoints.
  S_linpdf <- function(t) {
    tt <- pmax(pmin(t, omega), 0)
    1 - 0.000625 * tt^2
  }

  dq_lin <- S_linpdf(m_grid) - S_linpdf(m_grid + 1)
  dq_lin_check <- S_linpdf(m_check) - S_linpdf(m_check + 1)

  preview <- data.frame(
    m = m_grid,
    dq_unif = dq_unif,
    dq_exp  = dq_exp,
    dq_lin  = dq_lin
  )

  if (make_plot) {
    matplot(
      preview$m,
      cbind(preview$dq_unif, preview$dq_exp, preview$dq_lin),
      type = "l", lty = 1:3, lwd = 2,
      xlab = "m",
      ylab = "m|q0 = S0(m) - S0(m+1)",
      main = "Exercise 5-8: Behavior of m|q0 across m"
    )
    legend(
      "topright",
      legend = c("Uniform", "Exponential", "f0(t)=0.00125 t on [0,omega]"),
      lty = 1:3, lwd = 2, bty = "n"
    )
    grid()
    if (include_endpoint_in_preview) {
      abline(v = omega, lty = 3)  # highlight boundary
    }
  }

  list(
    uniform_is_constant        = max(abs(diff(dq_unif_check))) < 1e-10,
    exponential_is_decreasing  = all(diff(dq_exp_check) <= 1e-10),
    linpdf_is_increasing       = all(diff(dq_lin_check) >= -1e-10),
    preview = preview
  )
}

# Example run:
out <- check_ex5_8(make_plot = TRUE)
out$uniform_is_constant
out$exponential_is_decreasing
out$linpdf_is_increasing
tail(out$preview)




############################################################
## Exercise 5-12: S0(t) = (9000 - 10t - t^2)/9000, compute q50 - mu50
############################################################
check_ex5_12 <- function() {

  S0 <- function(t) (9000 - 10 * t - t^2) / 9000
  Sprime <- function(t) (-10 - 2 * t) / 9000

  S50 <- S0(50); S51 <- S0(51)

  q50 <- (S50 - S51) / S50
  mu50 <- (-Sprime(50)) / S50

  list(
    S50 = S50, S51 = S51,
    q50 = q50,
    mu50 = mu50,
    diff = q50 - mu50
  )
}

print(check_ex5_12())

############################################################
## Exercise 5-15: mu_y = (80 - y)^(-1/2), median of T20 solves p20(t)=0.5
############################################################
check_ex5_15 <- function() {

  mu <- function(y) (80 - y)^(-1/2)

  # survival from age 20 for duration t:
  p20 <- function(t) {
    # integrate mu from 20 to 20+t
    val <- num_int(mu, 20, 20 + t)
    exp(-val)
  }

  t_med <- uniroot(function(t) p20(t) - 0.5, lower = 0, upper = 59.999)$root

  list(
    median_T20_numeric = t_med,
    p20_at_median = p20(t_med)
  )
}

print(check_ex5_15())


############################################################
## Exercise 5-18: select model S_[x](t;x)=1 - t/(40-x)
############################################################
check_ex5_18 <- function() {

  S_sel <- function(t, x) 1 - t / (40 - x)

  # (a) _4 p_[30]
  p4_30 <- S_sel(4, 30)

  # (b) e_[30] = integral_0^(40-30) S(t;30) dt
  e30 <- num_int(function(t) S_sel(t, 30), 0, 10)

  # (c) mu_[20]+t = -d/dt log S(t;20) = 1/(20 - t)
  mu20t <- function(t) 1 / (20 - t)

  tgrid <- seq(0, 19, by = 0.5)
  # numeric hazard from S:
  S20 <- function(t) S_sel(t, 20)
  mu_num <- -fd_deriv(function(u) log(S20(u)), tgrid)

  list(
    p4_30 = p4_30,
    e_30_numeric = e30,
    mu_identity_check = data.frame(t = tgrid, mu_theory = mu20t(tgrid), mu_numeric = mu_num)
  )
}

out518 <- check_ex5_18()
print(out518$p4_30)
print(out518$e_30_numeric)
head(out518$mu_identity_check)
