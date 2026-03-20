## Listing 11.1: Full preliminary term premiums and their benchmark values
library(mqriskR)

x <- 40; i <- 0.05

alpha_f <- alphaF(x = x,i = i,model = "uniform",omega = 100)

beta_f <- betaF(x = x,i = i,model = "uniform",omega = 100)

first_year_cost <- Axn1(x = x,n = 1,i = i,model = "uniform",omega = 100)

attained_age_premium <- Px(x = x + 1,i = i,model = "uniform",omega = 100)

c(
  alpha_F = alpha_f,
  first_year_cost = first_year_cost,
  beta_F = beta_f,
  attained_age_premium = attained_age_premium
)



## Listing 11.2: Comparison of FPT and net level premium reserves
library(mqriskR)

x <- 40
i <- 0.05
t <- 5

fpt_reserve <- tVFx(
  x = x,
  t = t,
  i = i,
  model = "uniform",
  omega = 100
)

shifted_nlp <- tVx(
  x = x + 1,
  t = t - 1,
  i = i,
  model = "uniform",
  omega = 100
)

c(
  fpt_reserve = fpt_reserve,
  shifted_nlp = shifted_nlp
)


## Listing 11.3:Fractional-duration reserve and mean reserve
library(mqriskR)

x <- 40; t <- 10; s <- 0.5; i <- 0.05

# Package value
fractional_reserve <- tsVx(x = x, t = t, s = s, i = i, model = "uniform", omega = 100 )

# Direct calculation from Equation (11.12)
Vt <- tVx(x = x, t = t, i = i, model = "uniform", omega = 100 )

Vt1 <- tVx( x = x, t = t + 1, i = i, model = "uniform", omega = 100 )

P <- Px(x = x,i = i, model = "uniform", omega = 100)

fractional_formula <- Vt * (1 - s) + Vt1 * s + P * (1 - s)

# Mean reserve
mean_reserve <- meanVx( x = x, t = t, i = i, model = "uniform", omega = 100)

c(reserve_t_plus_s = fractional_reserve,
  formula_11_12 = fractional_formula,
  mean_reserve = mean_reserve)


## Listing 11.4: Gross premium reserve and expense reserve
library(mqriskR)

x <- 40; t <- 10; i <- 0.05; G <- 0.04

gross_reserve <- tVGx( x = x, t = t, i = i, G = G, benefit = 1, renewal_premium_pct = 0.10, renewal_policy_exp = 0.002, settlement_exp = 0.02, model = "uniform", omega = 100 )

expense_reserve <- tVEx( x = x, t = t, i = i, G = G, benefit = 1, renewal_premium_pct = 0.10, renewal_policy_exp = 0.002, settlement_exp = 0.02, model = "uniform", omega = 100)

net_reserve <- tVx( x = x, t = t, i = i, model = "uniform", omega = 100)

c(gross_reserve = gross_reserve,
  net_reserve = net_reserve,
  expense_reserve = expense_reserve,
  check = net_reserve + expense_reserve)




## Listing 11.5: Ordered gross gain decomposition
library(mqriskR)

out <- decompGg_disc(VtG = 3950.73, Vt1G = 4607.07, G = 685.00,
                     i_assumed = 0.06, q_assumed = 0.00592,
                     r_assumed = 0.05, e_assumed = 0,
                     s_assumed = 300, i_actual = 0.065,
                     q_actual = 0.005, r_actual = 0.06,
                     e_actual = 0, s_actual = 100,
                     b = 50000, order = c("interest", "mortality", "expense"))

out



