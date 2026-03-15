## Listing 10.1 Prospective reserves for several insurance contracts
library(mqriskR)

x <- 40
n <- 20
t <- 10
i <- 0.05

whole_life <- tVx(
  x = x,
  t = t,
  i = i,
  model = "uniform",
  omega = 100
)

term_insurance <- tVxn1(
  x = x,
  n = n,
  t = t,
  i = i,
  model = "uniform",
  omega = 100
)

endowment_insurance <- tVxn(
  x = x,
  n = n,
  t = t,
  i = i,
  model = "uniform",
  omega = 100
)

c(
  whole_life = whole_life,
  term_insurance = term_insurance,
  endowment_insurance = endowment_insurance
)


## Listing 10.2: Comparison of prospective and retrospective whole life reserves
library(mqriskR)

x <- 40; t <- 10; i <- 0.05

prospective <- tVx(x = x,t = t,i = i,model = "uniform",omega = 100)

retrospective <- tVx_ret(x = x,t = t,i = i,model = "uniform", omega = 100)

c(prospective = prospective,
  retrospective = retrospective,
  difference = prospective - retrospective)




## Listing 10.3: Equivalent formulas for the whole life reserve
library(mqriskR)

x <- 40; t <- 10; i <- 0.05

reserve_direct <- tVx(x = x,t = t,i = i,model = "uniform", omega = 100)

reserve_annuity_ratio <- 1 - adotx(x = x + t,i = i,model = "uniform",omega = 100) / adotx(x = x,i = i,model = "uniform", omega = 100)

reserve_premium_diff <- (Px(x = x + t,i = i,model = "uniform",omega = 100) -Px(x = x,i = i,model = "uniform",omega = 100)) * adotx(x = x + t,i = i,model = "uniform",omega = 100)

c(
  direct = reserve_direct,
  annuity_ratio = reserve_annuity_ratio,
  premium_difference = reserve_premium_diff
)


## Listing 10.4: Conditional mean and variance of loss at duration \(t\)
library(mqriskR)

x <- 40
t <- 10
i <- 0.05

premium <- Px(
  x = x,
  i = i,
  model = "uniform",
  omega = 100
)

reserve <- tVx(
  x = x,
  t = t,
  i = i,
  model = "uniform",
  omega = 100
)

mean_loss <- ELtx(
  x = x,
  t = t,
  i = i,
  P = premium,
  model = "uniform",
  omega = 100
)

var_loss <- varLtx(
  x = x,
  t = t,
  i = i,
  P = premium,
  model = "uniform",
  omega = 100
)

c(
  reserve = reserve,
  expected_loss = mean_loss,
  variance_loss = var_loss
)



## Listing 10.5: Reserve for a whole life insurance with immediate payment of claims
library(mqriskR)

x <- 40; t <- 10; i <- 0.05

reserve_val <- tVbarx(x = x,t = t,i = i,model = "uniform", omega = 100)

reserve_val



## Listing 10.6: Comparison of annual and monthly whole life reserves
library(mqriskR)

x <- 40; t <- 10; i <- 0.05

reserve_annual <- tVx(x = x,t = t,i = i,model = "uniform",omega = 100)

reserve_monthly <- tVx_m(x = x,t = t,m = 12,i = i, model = "uniform",omega = 100)

c(annual = reserve_annual,
  monthly = reserve_monthly)


## Listing 10.7: Gain decomposition for a discrete insurance contract
library(mqriskR)

Vt <- 0.085697; Vt1 <- 0.096115
P <- 0.008013; B <- 1
i_assumed <- 0.06; i_actual <- 0.065
q_assumed <- 0.00356; q_actual <- 0.00300


gt <- GT_disc(Vt = Vt, Vt1 = Vt1, P = P, i_actual = i_actual, q_actual = q_actual, B = B )

gm <- GM_disc(Vt = Vt, Vt1 = Vt1, P = P, i_assumed = i_assumed, q_actual = q_actual, B = B)

gi <- GI_disc( Vt = Vt, Vt1 = Vt1, P = P, i_actual = i_actual, q_assumed = q_assumed, B = B)

round(c(
  total_gain = gt,
  mortality_gain = gm,
  interest_gain = gi,
  check = gm + gi
), 7)


## Listing 10.8: Backward Euler approximation of a fully continuous reserve path
library(mqriskR)

times <- seq(19, 20, by = 0.25)

reserve_path <- thiele_backward_path(
  times = times,
  V_terminal = 1000,
  P = 26.96,
  delta = 0.0582689,
  mu = 0.00197,
  benefit = 1000
)

data.frame(
  time = times,
  reserve = reserve_path
)


