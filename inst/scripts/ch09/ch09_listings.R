## Listing 9.1: Net annual premiums for standard contingent contracts
library(mqriskR)

x <- 40
n <- 20
i <- 0.05

Px_val   <- Px(x = x, i = i, model = "uniform", omega = 100)
Pxn1_val <- Pxn1(x = x, n = n, i = i, model = "uniform", omega = 100)
PnEx_val <- PnEx(x = x, n = n, i = i, model = "uniform", omega = 100)
Pxn_val  <- Pxn(x = x, n = n, i = i, model = "uniform", omega = 100)

out <- c(
  Px   = Px_val,
  Pxn1 = Pxn1_val,
  PnEx = PnEx_val,
  Pxn  = Pxn_val
)

print(round(out, 6))

# Check the endowment premium decomposition
round(Pxn_val - (Pxn1_val + PnEx_val), 10)




## Listing 9.2: Comparison of full-payment and limited-payment net premiums
library(mqriskR)

x <- 40; n <- 20; t <- 10; i <- 0.05

whole_full   <- Px(x = x, i = i, model = "uniform", omega = 100)
whole_10pay  <- tPx(x = x, t = t, i = i, model = "uniform", omega = 100)

term_full    <- Pxn1(x = x, n = n, i = i, model = "uniform", omega = 100)
term_10pay   <- tPxn1(x = x, n = n, t = t, i = i, model = "uniform", omega = 100)

endow_full   <- Pxn(x = x, n = n, i = i, model = "uniform", omega = 100)
endow_10pay  <- tPxn(x = x, n = n, t = t, i = i, model = "uniform", omega = 100)

out <- c(
  whole_full  = whole_full,
  whole_10pay = whole_10pay,
  term_full   = term_full,
  term_10pay  = term_10pay,
  endow_full  = endow_full,
  endow_10pay = endow_10pay
)

print(round(out, 6))



## Listing 9.3: Present value of loss for a whole life insurance
library(mqriskR)

x <- 40
i <- 0.05

premium <- Px(x = x, i = i, model = "uniform", omega = 100)

mean_loss <- EL0x(x = x, i = i, P = premium,
  model = "uniform", omega = 100
)

var_loss <- varL0x(
  x = x, i = i, P = premium,
  model = "uniform", omega = 100
)

c(
  premium = premium,
  mean_loss = mean_loss,
  variance_loss = var_loss
)



## Listing 9.4: Normal approximation for aggregate present value of loss
Ax <- 0.24905
A2x <- 0.09476
i <- 0.06
P <- 0.025

d <- i / (1 + i)
c_factor <- 1 + P / d

mean_L <- c_factor * Ax - P / d
var_L  <- c_factor^2 * (A2x - Ax^2)

z_95 <- qnorm(0.95)
n_required <- (z_95 * sqrt(var_L) / abs(mean_L))^2

c(
  mean_single = mean_L,
  var_single = var_L,
  n_required = ceiling(n_required)
)



## Listing 9.5: Fully continuous premium rate for a whole life insurance

x <- 50
i <- 0.05

premium_rate <- PbarAbarx(x = x, i = i, model = "exponential", lambda = 0.02)
abar_val <- abarx(x = x, i = i, model = "exponential", lambda = 0.02)
delta <- interest_convert(i = i)$delta

check_val <- 1 / abar_val - delta

c(
  premium_rate = premium_rate,
  identity_value = check_val
)



## Listing 9.6:True fractional premium rates for a term insurance
library(mqriskR)

x <- 30
n <- 10
i <- 0.05

P_annual <- Pxn1(x = x, n = n, i = i, model = "uniform", omega = 100)
P_semi   <- Pxn1_m(x = x, n = n, m = 2, i = i, model = "uniform", omega = 100)
P_month  <- Pxn1_m(x = x, n = n, m = 12, i = i, model = "uniform", omega = 100)

c(
  annual = P_annual,
  semiannual = P_semi,
  monthly = P_month
)


## Listing 9.7: Net and gross premiums for a whole life insurance with expenses
library(mqriskR)

x <- 40
i <- 0.05

# Net premium per unit benefit
net_unit <- Px(
  x = x,
  i = i,
  model = "uniform",
  omega = 100
)

# Net premium for a benefit of 1000
net_premium <- 1000 * net_unit

# Gross premium including expenses
gross_premium <- Gx(
  x = x,
  i = i,
  benefit = 1000,
  first_premium_pct = 0.75,
  renewal_premium_pct = 0.10,
  first_policy_exp = 10,
  renewal_policy_exp = 2,
  settlement_exp = 20,
  model = "uniform",
  omega = 100
)

c(
  net_premium = net_premium,
  gross_premium = gross_premium,
  loading = gross_premium - net_premium
)


