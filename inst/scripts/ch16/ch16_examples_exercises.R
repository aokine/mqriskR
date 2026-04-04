# Chapter 16 R Check: Example 16.7
# Catch-up premium under stipulated premium guarantee

library(mqriskR)

stip_premium <- 10000
years_elapsed <- 10
cum_paid <- 90000

required_to_date <- stip_premium * years_elapsed
shortfall <- required_to_date - cum_paid
payment_needed_now <- shortfall + stip_premium

out <- data.frame(
  years_elapsed = years_elapsed,
  required_to_date = required_to_date,
  cumulative_paid = cum_paid,
  shortfall = shortfall,
  payment_needed_now = payment_needed_now
)

print(out)

# Simple plot
barplot(
  height = c(required_to_date, cum_paid),
  names.arg = c("Required cumulative", "Actual cumulative"),
  main = "Example 16.7: Stipulated Premium Catch-Up",
  ylab = "Amount"
)
abline(h = required_to_date, lty = 2)



# Chapter 16 R Check: Example 16.8
# Lapse timing with and without shadow fund guarantee

library(mqriskR)

yr <- 1:6
contract_av <- c(1000, 800, 400, 100, 0, 0)
surr_charge <- c(900, 800, 600, 400, 200, 100)
cash_value <- c(100, 0, 0, 0, 0, 0)
shadow_av <- c(1100, 900, 700, 500, 200, 0)

lapse_no_guarantee <- yr[min(which(contract_av == 0))]
lapse_shadow <- yr[min(which(shadow_av == 0))]

out <- data.frame(
  year = yr,
  contract_av = contract_av,
  surrender_charge = surr_charge,
  cash_value = cash_value,
  shadow_av = shadow_av
)

print(out)
cat("Lapse year without secondary guarantee:", lapse_no_guarantee, "\n")
cat("Lapse year with shadow-fund secondary guarantee:", lapse_shadow, "\n")




# Chapter 16 R Check: Example 16.9
# Basic UL reserving quantities using package helpers where available

library(mqriskR)

face_per_1000 <- 100
AV9_per_1000 <- 57.60
GMF9_per_1000 <- 140.40
gmp_per_1000 <- 14.49
expense_rate <- 0.04
guar_charge10_per_1000 <- 11.80
curr_charge10_per_1000 <- 10.76
curr_i <- 0.05
guar_i <- 0.03
pvfb_minus_pvfp_per_1000 <- 70

# Account value roll-forward (computed directly from the example)
AV10_per_1000 <- (AV9_per_1000 - curr_charge10_per_1000) * (1 + curr_i)

# Package helpers
GMF10_per_1000 <- GMF_rollforward_ul(
  GMF_prev = GMF9_per_1000,
  GMP = gmp_per_1000,
  r = expense_rate,
  policy_charge = guar_charge10_per_1000,
  i = guar_i
)

AV10 <- AV10_per_1000 * face_per_1000
GMF10 <- GMF10_per_1000 * face_per_1000

r10 <- rt_ul(AV = AV10, GMF = GMF10)

pre_floor_CRVM_per_1000 <- Vprefloor_crvm_ul(
  r = r10,
  pvfb_minus_pvfp = pvfb_minus_pvfp_per_1000
)

pre_floor_CRVM <- pre_floor_CRVM_per_1000 * face_per_1000

out <- data.frame(
  AV10_per_1000 = AV10_per_1000,
  GMF10_per_1000 = GMF10_per_1000,
  AV10 = AV10,
  GMF10 = GMF10,
  r10 = r10,
  pre_floor_CRVM_per_1000 = pre_floor_CRVM_per_1000,
  pre_floor_CRVM = pre_floor_CRVM
)

print(round(out, 5))

barplot(
  height = c(AV10, GMF10, pre_floor_CRVM),
  names.arg = c("AV10", "GMF10", "Pre-floor CRVM"),
  main = "Example 16.9: UL Reserving Quantities",
  ylab = "Amount"
)




# Chapter 16 R Check: Example 16.10
# Implied guaranteed rates for IUL

library(mqriskR)

valuation_rate <- 0.04
option_cost_initial <- 0.05
option_cost_renewal <- 0.02
guaranteed_floor <- 0.00

igr_initial <- guaranteed_floor + option_cost_initial * (1 + valuation_rate)
igr_renewal <- guaranteed_floor + option_cost_renewal * (1 + valuation_rate)

out <- data.frame(
  implied_guaranteed_rate_initial = igr_initial,
  implied_guaranteed_rate_renewal = igr_renewal
)

print(out)

barplot(
  height = c(igr_initial, igr_renewal),
  names.arg = c("Initial term", "Beyond initial term"),
  main = "Example 16.10: Implied Guaranteed Rates",
  ylab = "Rate"
)




# Chapter 16 R Check: Example 16.11
# AG 38 reserve under shadow fund method using package helpers

library(mqriskR)

shadow_fund <- 60000
nsp_full_funding <- 100000
valuation_nsp <- 150000
surrender_charge <- 5000
basic_reserve <- 10000
deficiency_reserve <- 0

prefunding_ratio <- ag38_prefunding_ratio(
  excess_payment = shadow_fund,
  nsp_required = nsp_full_funding
)

ag38_out <- ag38_reserve_ul(
  basic_reserve = basic_reserve,
  deficiency_reserve = deficiency_reserve,
  excess_payment = shadow_fund,
  nsp_required = nsp_full_funding,
  valuation_nsp = valuation_nsp,
  surrender_charge = surrender_charge
)

out <- data.frame(
  prefunding_ratio = prefunding_ratio,
  net_amount_additional = ag38_out$net_amount_additional,
  reduced_deficiency_reserve = ag38_out$reduced_deficiency_reserve,
  reserve_step8 = ag38_out$step8_reserve,
  ag38_reserve = ag38_out$final_reserve
)

print(out)

barplot(
  height = c(shadow_fund, nsp_full_funding, ag38_out$final_reserve),
  names.arg = c("Shadow fund", "Full-funding NSP", "AG38 reserve"),
  main = "Example 16.11: AG 38 Reserve Components",
  ylab = "Amount"
)


# Chapter 16 R Check: Exercise 16-3
# Corridor-factor test for Type A universal life

library(mqriskR)

B <- 100000
G <- 1500
r <- 0.025
e <- 50
qx <- 0.00278
iq <- 0.06
ic <- 0.07
f <- 2.50
AV_prev <- 40000

fund_before_coi <- AV_prev + G * (1 - r) - e

# Ordinary Type A COI from the formula in Exercise 16-3(a)
coi_typeA <- (
  (qx / (1 + iq)) * (B - fund_before_coi * (1 + ic))
) / (
  1 - (qx / (1 + iq)) * (1 + ic)
)

# Corridor-factor COI from the formula in Exercise 16-3(a)
coi_cf <- (
  (qx / (1 + iq)) * fund_before_coi * (1 + ic) * (f - 1)
) / (
  1 + (qx / (1 + iq)) * (1 + ic) * (f - 1)
)

coi_used <- max(coi_typeA, coi_cf)
corridor_applies <- coi_cf > coi_typeA

AV_t <- (fund_before_coi - coi_used) * (1 + ic)

out <- data.frame(
  coi_typeA = coi_typeA,
  coi_cf = coi_cf,
  coi_used = coi_used,
  corridor_applies = corridor_applies,
  AV_t = AV_t
)

print(round(out, 5))

barplot(
  height = c(coi_typeA, coi_cf),
  names.arg = c("Ordinary COI", "Corridor COI"),
  main = "Exercise 16-3: Competing COI Values",
  ylab = "Cost of insurance"
)




# Chapter 16 R Check: Exercise 16-8
# Indexed UL credited rates and cash values using package helpers

library(mqriskR)

B <- 100000
G <- rep(1000, 3)
premium_expense <- 0.04
admin_fee <- 50
index <- c(1000, 1080, 1200, 1100)
coi_per_1000 <- c(2.0, 3.0, 4.0)
surr_per_1000 <- c(5.0, 4.0, 3.0)

raw <- iP_eiul(index)
credited <- i_credit_eiul(
  i_raw = raw,
  part = 1.00,
  floor = 0.01,
  cap = 0.10
)

AV <- numeric(3)
cash_value <- numeric(3)

for (t in 1:3) {
  prior <- if (t == 1) 0 else AV[t - 1]
  net_premium <- G[t] * (1 - premium_expense)
  amount_at_risk <- B - (prior + net_premium)
  coi <- (coi_per_1000[t] / 1000) * amount_at_risk
  AV[t] <- (prior + net_premium - coi - admin_fee) * (1 + credited[t])
  cash_value[t] <- AV[t] - surr_per_1000[t] * (B / 1000)
}

out <- data.frame(
  year = 1:3,
  raw = raw,
  credited = credited,
  AV = AV,
  cash_value = cash_value
)

print(round(out, 5))

plot(
  1:3, raw, type = "b", pch = 16,
  ylim = range(c(raw, credited)),
  xlab = "Year", ylab = "Rate",
  main = "Exercise 16-8: Raw vs Credited Rates"
)
lines(1:3, credited, type = "b", pch = 17)
abline(h = 0.10, lty = 2)
abline(h = 0.01, lty = 2)
legend(
  "topright",
  legend = c("Raw index growth", "Credited rate"),
  lty = 1, pch = c(16, 17), bty = "n"
)



# Chapter 16 R Check: Exercise 16-9
# APV of premiums and pure endowment under mortality and withdrawal

library(mqriskR)

qd <- c(.001, .002, .003, .004, .005)
qw <- c(.02, .02, .03, .04, .05)
G <- c(20000, 25000, 25000, 30000, 20000)
i <- 0.05
B <- 100000

p_tau <- pxtau_ul(qd = qd, qw = qw, year_end_withdrawal = TRUE)
tp_tau <- tpxtau_ul(qd = qd, qw = qw, year_end_withdrawal = TRUE)

apvp <- G[1] +
  G[2] * tp_tau[1] / (1 + i)^1 +
  G[3] * tp_tau[2] / (1 + i)^2 +
  G[4] * tp_tau[3] / (1 + i)^3 +
  G[5] * tp_tau[4] / (1 + i)^4

apvb <- B * tp_tau[5] / (1 + i)^5

out <- data.frame(
  year = 1:5,
  p_tau = p_tau,
  tp_tau = tp_tau
)

print(round(out, 5))
cat("APV of premiums =", round(apvp, 2), "\n")
cat("APV of pure endowment =", round(apvb, 2), "\n")

plot(
  1:5, tp_tau, type = "b", pch = 16,
  xlab = "Policy year", ylab = "Cumulative in-force probability",
  main = "Exercise 16-9: Cumulative Persistency"
)


# Chapter 16 R Check: Exercise 16-11
# Revised AV10, GMF10, and r10 after extra premium using package helpers

library(mqriskR)

face_per_1000 <- 100
AV10_old <- 4918.20
extra_premium <- 1000
expense_rate <- 0.04
curr_i <- 0.05

GMF9_per_1000 <- 140.40
GMP_per_1000 <- 14.49
guar_expense_rate <- 0.04
guar_charge10_per_1000 <- 11.80
guar_i <- 0.03

AV10_new <- AV10_old + extra_premium * (1 - expense_rate) * (1 + curr_i)

GMF10_per_1000 <- GMF_rollforward_ul(
  GMF_prev = GMF9_per_1000,
  GMP = GMP_per_1000,
  r = guar_expense_rate,
  policy_charge = guar_charge10_per_1000,
  i = guar_i
)

GMF10 <- GMF10_per_1000 * face_per_1000

r10 <- rt_ul(
  AV = AV10_new,
  GMF = GMF10
)

out <- data.frame(
  AV10_new = AV10_new,
  GMF10 = GMF10,
  r10 = r10
)

print(round(out, 5))

barplot(
  height = c(AV10_new, GMF10),
  names.arg = c("AV10", "GMF10"),
  main = "Exercise 16-11: Revised AV and GMF",
  ylab = "Amount"
)



