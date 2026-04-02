## Listing 15.1: Pure endowment under variable interest scenarios
library(mqriskR)

qx <- c(.03, .04, .05, .06, .07)

scen1 <- c(.06, .07, .08, .09, .10)
scen2 <- c(.06, .06, .06, .06, .06)
scen3 <- c(.06, .05, .04, .03, .02)

c(Scenario1 = nEx_var(qx = qx, i = scen1, benefit = 1000),
  Scenario2 = nEx_var(qx = qx, i = scen2, benefit = 1000),
  Scenario3 = nEx_var(qx = qx, i = scen3, benefit = 1000))



## Listing 5.2: Term insurance under variable interest scenarios
library(mqriskR)

qx <- c(.03, .04, .05, .06, .07)

scen1 <- c(.06, .07, .08, .09, .10)
scen2 <- c(.06, .06, .06, .06, .06)
scen3 <- c(.06, .05, .04, .03, .02)

c(Scenario1 = Axn1_var(qx = qx, i = scen1),
  Scenario2 = Axn1_var(qx = qx, i = scen2),
  Scenario3 = Axn1_var(qx = qx, i = scen3))



##Listing 15.3: Deterministic interest-rate stress test for a temporary annuity
library(mqriskR)

qx <- rep(.02, 5)

level <- c(.06, .06, .06, .06, .06)
rise  <- c(.06, .07, .08, .09, .10)
fall  <- c(.06, .05, .04, .03, .03)

apv_level <- axn_var(qx = qx, i = level, type = "immediate")
apv_rise  <- axn_var(qx = qx, i = rise,  type = "immediate")
apv_fall  <- axn_var(qx = qx, i = fall,  type = "immediate")

c(level = apv_level,
  rise  = apv_rise,
  fall  = apv_fall,
  ratio_fall_to_level = apv_fall / apv_level)



##Listing 15.4: Bootstrapping semiannual spot rates
library(mqriskR)

maturity <- c(0.5, 1.0, 1.5, 2.0)
coupon_yield <- c(0.0244, 0.0260, 0.0276, 0.0293)

z_from_coupon_semi(maturity = maturity, coupon_yield = coupon_yield)



##Listing 15.5: Present value using spot rates
library(mqriskR)

pv_spot_cashflows(
  amounts = c(200000, 50000, 50000, 100000),
  times   = c(0, 0.5, 1.0, 2.0),
  spot    = c(0, 0.02440, 0.02601, 0.02936),
  compounding = "semiannual")



##Listing 5.6: Net premium using spot rates
library(mqriskR)

qx <- c(.02, .03, .04, .05, .06)
z  <- c(.03, .04, .05, .06, .07)

apvb <- Axn1_spot(qx = qx, z = z)
apvp <- axn_spot(qx = qx, z = z, type = "due")

premium <- 1000000 * apvb / apvp

c(APV_benefit = apvb,
  APV_premium_factor = apvp,
  Premium = premium)



##Listing 15.7 Forward rates implied by spot rates}, label={lst:ch15_forward_rates}]
library(mqriskR)

z <- c(0.03, 0.04, 0.05, 0.06, 0.07)

c(f_1_4 = fnk_from_z(z, n = 1, k = 4),
  f_2_2 = fnk_from_z(z, n = 2, k = 2))

forward_matrix_from_z(z)


##Listing 15.8: Spot rates recovered from forward one-year rates
library(mqriskR)

fn1 <- c(0.04, 0.05, 0.06, 0.07, 0.08)

z_from_fn1(fn1)


