# R Check for Example 9.1
# Special 10-year contract with refund of net annual premiums on survival

A_endow <- 0.60
A_pure  <- 0.47
d <- 0.05

A_term <- A_endow - A_pure
adot_30_10 <- (1 - A_endow) / d
P <- A_term / (adot_30_10 - 10 * A_pure)

round(c(
  A_term = A_term,
  adot_30_10 = adot_30_10,
  premium = P
), 5)



# R Check for Example 9.4
# Variance of present-value-of-loss random variable for a 2-year term insurance

qx <- 0.10
qx1 <- 0.20
px <- 1 - qx
v <- 0.90

A_x_2_term <- v * qx + v^2 * px * qx1
adot_x_2 <- 1 + v * px
P <- A_x_2_term / adot_x_2

L1 <- v - P
L2 <- v^2 - P * (1 + v)
L3 <- -P * (1 + v)

pr1 <- qx
pr2 <- px * qx1
pr3 <- px * (1 - qx1)

EL <- pr1 * L1 + pr2 * L2 + pr3 * L3
VarL <- pr1 * L1^2 + pr2 * L2^2 + pr3 * L3^2 - EL^2

round(c(
  premium = P,
  loss_year1 = L1,
  loss_year2 = L2,
  loss_no_claim = L3,
  expected_loss = EL,
  variance = VarL
), 5)



# R Check for Example 9.6
# Expected value of loss under elevated first-year mortality

i <- 0.06
d <- i / (1 + i)

A60 <- 0.36933
A61 <- 0.38300
q60 <- 0.01376

adot60 <- (1 - A60) / d
P60 <- 1000 * A60 / adot60

q60_imp <- 10 * q60
p60_imp <- 1 - q60_imp

A60_imp <- (q60_imp / (1 + i)) + (p60_imp / (1 + i)) * A61
adot60_imp <- (1 - A60_imp) / d

EL_imp <- 1000 * A60_imp - P60 * adot60_imp

round(c(
  adot60 = adot60,
  premium = P60,
  A60_imp = A60_imp,
  adot60_imp = adot60_imp,
  expected_loss = EL_imp
), 5)





# R Check for Example 9.11
# Standard and impaired-life temporary annuity and term insurance values

i <- 0.06
v <- 1 / (1 + i)
d <- i / (1 + i)

tp <- c(1.00000, 0.99105, 0.98143, 0.97107, 0.95994, 0.94800)

# Part (a)
adot_temp5 <- sum(v^(0:4) * tp[1:5])
nEx5 <- v^5 * tp[6]
A_term5 <- 1 - d * adot_temp5 - nEx5
P_term5 <- A_term5 / adot_temp5

round(c(
  adot_temp5 = adot_temp5,
  A_term5 = A_term5,
  P_term5 = P_term5
), 5)

# Part (b): impaired life with mu' = mu + 0.002
delta <- log(1 + i)
delta_prime <- delta + 0.002
i_prime <- exp(delta_prime) - 1

adot_temp5_prime <- sum((1 / (1 + i_prime))^(0:4) * tp[1:5])
nEx5_prime <- (1 / (1 + i_prime))^5 * tp[6]
A_term5_prime <- 1 - d * adot_temp5_prime - nEx5_prime
P_term5_prime <- A_term5_prime / adot_temp5_prime

round(c(
  i_prime = i_prime,
  adot_temp5_prime = adot_temp5_prime,
  A_term5_prime = A_term5_prime,
  P_term5_prime = P_term5_prime
), 5)


# Exercise 9.4 R Check: whole life premium with premiums increasing at interest

i <- 0.05
v <- 1 / (1 + i)
omega <- 105
x <- 75

# 1 + e_75 under uniform future lifetime over 30 years
t_vals <- 0:(omega - x - 1)
surv_vals <- (omega - x - t_vals) / (omega - x)
one_plus_ex <- sum(surv_vals)

# Whole life insurance APV under the same model
k_vals <- 0:(omega - x - 1)
kq_vals <- rep(1 / (omega - x), length(k_vals))
Ax_75 <- sum(v^(k_vals + 1) * kq_vals)

P1 <- 1000 * Ax_75 / one_plus_ex

round(c(
  one_plus_ex = one_plus_ex,
  Ax_75 = Ax_75,
  P1 = P1
), 5)


# Exercise 9.5 R Check: mortgage-protection decreasing term premium

loan <- 100000
i <- 0.05
n <- 25
adot_40_25 <- 14
p25 <- 0.80

# Annual mortgage payment
mortgage_payment <- loan / annuity_certain(n = n, i = i, due = FALSE)

# Approximate insurance APV using the given closed-form solution
v <- 1 / (1 + i)
q25 <- 1 - p25

apv_benefit <- 149000 * (
  1 - (i / (1 + i)) * adot_40_25 - p25 * v^25 - q25 * v^26
)

premium <- apv_benefit / adot_40_25

round(c(
  mortgage_payment = mortgage_payment,
  apv_benefit = apv_benefit,
  premium = premium
), 5)



# Exercise 9.10 R Check: expected present value of loss from given variance and moments

Ax <- 0.29224
A2x <- 0.11723
VarL <- 0.10
i <- 0.05
d <- i / (1 + i)

VarZ <- A2x - Ax^2
P_over_d <- sqrt(VarL / VarZ) - 1

EL <- (1 + P_over_d) * Ax - P_over_d

round(c(
  VarZ = VarZ,
  P_over_d = P_over_d,
  expected_loss = EL
), 5)



# Exercise 9.13 R Check: aggregate-loss premium target and portfolio sensitivity

Ax <- 0.24905
A2x <- 0.09476
i <- 0.06
d <- i / (1 + i)
z95 <- qnorm(0.95)

prob_pos_loss <- function(P, n, Ax, A2x, d) {
  mean_i <- (1 + P / d) * Ax - P / d
  var_i <- (1 + P / d)^2 * (A2x - Ax^2)
  mean_agg <- n * mean_i
  sd_agg <- sqrt(n * var_i)
  1 - pnorm((0 - mean_agg) / sd_agg)
}

# Premium for n = 100 with target probability about 0.05
objective <- function(P) prob_pos_loss(P, n = 100, Ax = Ax, A2x = A2x, d = d) - 0.05
P_star_100 <- uniroot(objective, interval = c(0.015, 0.03))$root

# Limiting premium (equivalence principle)
P_limit <- d * Ax / (1 - Ax)

# Probability for P = 0.024 and n = 100
prob_at_024 <- prob_pos_loss(0.024, n = 100, Ax = Ax, A2x = A2x, d = d)

round(c(
  P_star_100 = P_star_100,
  P_limit = P_limit,
  prob_at_024 = prob_at_024
), 5)

# Plot probability of positive aggregate loss against n for a few premium values
n_vals <- 10:300
p1 <- sapply(n_vals, prob_pos_loss, P = 0.02188, Ax = Ax, A2x = A2x, d = d)
p2 <- sapply(n_vals, prob_pos_loss, P = 0.02400, Ax = Ax, A2x = A2x, d = d)
p3 <- sapply(n_vals, prob_pos_loss, P = 0.02500, Ax = Ax, A2x = A2x, d = d)

plot(n_vals, p1, type = "l", lwd = 2,
     xlab = "Number of contracts",
     ylab = "Pr(Aggregate loss > 0)",
     main = "Aggregate loss probability versus portfolio size")
lines(n_vals, p2, lwd = 2, lty = 2)
lines(n_vals, p3, lwd = 2, lty = 3)
abline(h = 0.05, lty = 4)
legend("topright",
       legend = c("P = 0.02188", "P = 0.02400", "P = 0.02500", "Target 0.05"),
       lty = c(1, 2, 3, 4), lwd = c(2, 2, 2, 1), bty = "n")


# Exercise 9.25 R Check: monthly versus semiannual true fractional premiums under UDD

x# R Check: monthly versus semiannual true fractional premiums under UDD

A_term <- 0.015
adot_temp <- 8
nEx <- 0.604
i <- 0.05

adot_temp_m <- function(m) {
  ic <- interest_convert(i = i, m = m)
  im <- ic$im                      # nominal annual rate convertible m-thly
  d <- i / (1 + i)
  dm_small <- (im / m) / (1 + im / m)
  dm <- m * dm_small              # this is d^(m)

  alpha <- (i * d) / (im * dm)
  beta  <- (i - im) / (im * dm)

  alpha * adot_temp - beta * (1 - nEx)
}

adot2 <- adot_temp_m(2)
adot12 <- adot_temp_m(12)

prem2 <- 10000 * A_term / adot2
prem12 <- 10000 * A_term / adot12

round(c(
  adot2 = adot2,
  adot12 = adot12,
  prem2 = prem2,
  prem12 = prem12,
  difference = prem12 - prem2
), 5)



# Exercise 9.27 R Check: gross premium with expenses

adot_x_20 <- 10
adot_x <- 20
d <- 0.04

Ax <- 1 - d * adot_x

G <- (Ax + 0.05 + 0.02 * adot_x) / ((1 - 0.03) * adot_x_20)
P_net_20pay <- Ax / adot_x_20

round(c(
  Ax = Ax,
  net_premium = P_net_20pay,
  gross_premium = G
), 5)

