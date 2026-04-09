# Chapter 18 - R Check for Example 18.3
# Early-retirement reduction factors and immediate pension amounts

library(mqriskR)

# Example inputs
x <- 35
z <- 65
CASx <- 60000
p <- 2
g <- 0.04

# Retirement ages y = 60, 61, 62, 63, 64
y <- 60:64

# Salary scale helper
salary_at_age <- function(age, base_age = x, base_salary = CASx, growth = g) {
  base_salary * (1 + growth)^(age - base_age)
}

# Approximate projected FAS at age y + 1/2
# using the same weighting pattern as the text
FAS_half <- sapply(y, function(yy) {
  s_y3 <- salary_at_age(yy - 3)
  s_y2 <- salary_at_age(yy - 2)
  s_y1 <- salary_at_age(yy - 1)
  s_y  <- salary_at_age(yy)
  (0.5 * s_y3 + s_y2 + s_y1 + 0.5 * s_y) / 3
})

# Projected benefit payable at NRA
PAB_half <- 0.02 * (y + 0.5 - x) * FAS_half

# Immediate benefit using 5% per year early reduction
reduction_factor <- 1 - 0.05 * (z - (y + 0.5))
immediate_benefit <- PAB_half * reduction_factor

out_18_3 <- data.frame(
  y = y,
  retirement_age = y + 0.5,
  FAS_half = FAS_half,
  PAB_half = PAB_half,
  reduction_factor = reduction_factor,
  immediate_benefit = immediate_benefit
)

print(round(out_18_3, 2))

# Plot: reduction factor by retirement age
plot(
  out_18_3$retirement_age, out_18_3$reduction_factor,
  type = "b",
  xlab = "Retirement age",
  ylab = "Reduction factor",
  main = "Example 18.3: Early-retirement reduction factor"
)

# Plot: immediate benefit by retirement age
plot(
  out_18_3$retirement_age, out_18_3$immediate_benefit,
  type = "b",
  xlab = "Retirement age",
  ylab = "Immediate annual benefit",
  main = "Example 18.3: Immediate pension by retirement age"
)



# Chapter 18 - R Check for Example 18.4
# Withdrawal-benefit APV decomposition by withdrawal age

library(mqriskR)

# Illustrative setup consistent with Example 18.4 structure
x <- 35
z <- 65
CASx <- 60000
p <- 2
g <- 0.04
adue_65 <- 12

# Withdrawal ages in Example 18.4
y <- 40:59

# Illustrative active-service survival probabilities to age y
# (replace with service-table values if available)
p_active <- cumprod(c(1, rep(0.97, length(y) - 1)))

# Illustrative withdrawal probabilities by attained age
q_withdrawal <- seq(0.03, 0.01, length.out = length(y))

# Illustrative post-withdrawal survival to NRA
p_withdrawn_to_65 <- 0.98^(z - y - 0.5)

# Salary helper
salary_at_age <- function(age, base_age = x, base_salary = CASx, growth = g) {
  base_salary * (1 + growth)^(age - base_age)
}

# Approximate projected benefit at withdrawal age y + 1/2, deferred to NRA
FAS_half <- sapply(y, function(yy) {
  s_y3 <- salary_at_age(yy - 3)
  s_y2 <- salary_at_age(yy - 2)
  s_y1 <- salary_at_age(yy - 1)
  s_y  <- salary_at_age(yy)
  (0.5 * s_y3 + s_y2 + s_y1 + 0.5 * s_y) / 3
})

PAB_half <- 0.02 * (y + 0.5 - x) * FAS_half

# Discount factor from age 35 to age 65 as in Eq. (18.12)
v30 <- 1 / (1.06^30)

# Contribution by withdrawal age
apv_contrib <- PAB_half * v30 * p_active * q_withdrawal * p_withdrawn_to_65 * adue_65

out_18_4 <- data.frame(
  y = y,
  PAB_half = PAB_half,
  p_active = p_active,
  q_withdrawal = q_withdrawal,
  p_withdrawn_to_65 = p_withdrawn_to_65,
  apv_contribution = apv_contrib
)

print(round(out_18_4, 4))
cat("Total illustrative APV =", round(sum(apv_contrib), 4), "\n")

# Plot: APV contribution by withdrawal age
plot(
  y, apv_contrib,
  type = "h",
  lwd = 2,
  xlab = "Withdrawal age y",
  ylab = "APV contribution",
  main = "Example 18.4: Withdrawal-benefit APV contribution by age"
)



# Chapter 18 - R Check for Example 18.5
# TUC versus PUC normal cost and accrued liability

library(mqriskR)

projected_benefit <- 50000
past_service <- 26
total_service <- 35
accrued_benefit <- 20000
current_year_accrual <- 1500

i <- 0.06
p_surv <- 0.92
adue_65 <- 12
v_to_ret <- 1 / (1 + i)^9

NC_TUC <- NC_TUC_db(
  accrual_benefit = current_year_accrual,
  v_to_ret = v_to_ret,
  p_surv = p_surv,
  adue_ret = adue_65
)

AAL_TUC <- AAL_TUC_db(
  accrued_benefit = accrued_benefit,
  v_to_ret = v_to_ret,
  p_surv = p_surv,
  adue_ret = adue_65
)

NC_PUC <- NC_PUC_db(
  projected_benefit = projected_benefit,
  total_service = total_service,
  v_to_ret = v_to_ret,
  p_surv = p_surv,
  adue_ret = adue_65
)

AAL_PUC <- AAL_PUC_db(
  projected_benefit = projected_benefit,
  past_service = past_service,
  total_service = total_service,
  v_to_ret = v_to_ret,
  p_surv = p_surv,
  adue_ret = adue_65
)

out_18_5 <- c(
  NC_TUC = NC_TUC,
  AAL_TUC = AAL_TUC,
  NC_PUC = NC_PUC,
  AAL_PUC = AAL_PUC
)

print(round(out_18_5, 2))

# Barplot comparison
barplot(
  height = out_18_5,
  beside = TRUE,
  main = "Example 18.5: TUC versus PUC",
  ylab = "Amount"
)


# Chapter 18 - R Check for Exercise 18-1

library(mqriskR)

# Base calculation
AV_65 <- AVz_dc(
  x = 25, z = 65, Sx = 40000,
  c = 0.08, i = 0.05, g = 0.035
)

Income_65 <- Income_dc(
  AVz = AV_65,
  adue_z = 13
)

RR_65 <- replacement_ratio_dc(
  x = 25, z = 65, Sx = 40000,
  c = 0.08, i = 0.05,
  adue_z = 13, g = 0.035
)

out_18_1 <- c(
  AV_65 = AV_65,
  Income_65 = Income_65,
  RR_65 = RR_65
)

print(round(out_18_1, 5))

# Sensitivity of replacement ratio to salary growth
g_grid <- seq(0.01, 0.06, by = 0.005)

RR_grid <- sapply(g_grid, function(g) {
  replacement_ratio_dc(
    x = 25, z = 65, Sx = 40000,
    c = 0.08, i = 0.05,
    adue_z = 13, g = g
  )
})

plot(
  g_grid, RR_grid,
  type = "b",
  xlab = "Salary growth rate",
  ylab = "Replacement ratio",
  main = "Exercise 18-1: Replacement ratio vs salary growth"
)



# Chapter 18 - R Check for Exercise 18-2

library(mqriskR)

c_target <- contribution_rate_target(
  x = 30, z = 65, Sx = 60000,
  RR_target = 0.50,
  i = 0.06, adue_z = 11,
  g = 0.04
)

RR_check <- replacement_ratio_dc(
  x = 30, z = 65, Sx = 60000,
  c = c_target, i = 0.06,
  adue_z = 11, g = 0.04
)

out_18_2 <- c(
  c_target = c_target,
  RR_check = RR_check
)

print(round(out_18_2, 5))

# Sensitivity to interest rate
i_grid <- seq(0.03, 0.08, by = 0.005)

c_grid <- sapply(i_grid, function(ii) {
  contribution_rate_target(
    x = 30, z = 65, Sx = 60000,
    RR_target = 0.50,
    i = ii, adue_z = 11,
    g = 0.04
  )
})

plot(
  i_grid, c_grid,
  type = "b",
  xlab = "Fund interest rate",
  ylab = "Required contribution rate",
  main = "Exercise 18-2: Contribution rate vs interest rate"
)


# Chapter 18 - R Check for Exercise 18-3

library(mqriskR)

# Salary path with 4% regular increases and 6% merit increases
ages <- 30:64
salary <- numeric(length(ages))
salary[1] <- 100000

for (j in 2:length(ages)) {
  salary[j] <- salary[j - 1] * 1.04
  if (ages[j] %in% c(31, 32, 33)) {
    salary[j] <- salary[j] * 1.06
  }
}

names(salary) <- ages

FAS_65 <- mean(tail(salary, 5))
PAB_fas <- 0.01 * FAS_65 * 35
final_salary <- salary[length(salary)]
RR_fas <- PAB_fas / final_salary

CAP_65 <- 0.01 * sum(salary)
ratio_cae_to_fas <- CAP_65 / PAB_fas

out_18_3 <- c(
  FAS_65 = FAS_65,
  PAB_fas = PAB_fas,
  final_salary = final_salary,
  RR_fas = RR_fas,
  CAP_65 = CAP_65,
  ratio_cae_to_fas = ratio_cae_to_fas
)

print(round(out_18_3, 5))

# Plot salary path
plot(
  ages, salary,
  type = "b",
  xlab = "Age",
  ylab = "Salary",
  main = "Exercise 18-3: Salary path"
)

# Plot comparison of benefits
barplot(
  c(FAS_plan_benefit = PAB_fas, CAE_plan_benefit = CAP_65),
  main = "Exercise 18-3: Pension formula comparison",
  ylab = "Annual benefit"
)



# Chapter 18 - R Check for Exercise 18-4

library(mqriskR)

i <- 0.06

# Inputs from the written solution
ten_E_55 <- 0.51147
adue_65 <- (1 - 0.36778) / (i / (1 + i))
adue_55 <- (1 - 0.24766) / (i / (1 + i))

ERRF_55 <- ten_E_55 * adue_65 / adue_55
informal_factor <- 1 - 0.05 * 10
relative_diff <- (informal_factor - ERRF_55) / ERRF_55

out_18_4 <- c(
  adue_65 = adue_65,
  adue_55 = adue_55,
  ERRF_55 = ERRF_55,
  informal_factor = informal_factor,
  relative_diff = relative_diff
)

print(round(out_18_4, 5))

# Plot actuarial vs informal factor
barplot(
  c(actuarial = ERRF_55, informal = informal_factor),
  main = "Exercise 18-4: Early-retirement reduction factors",
  ylab = "Reduction factor"
)




# Chapter 18 - R Check for Exercise 18-6
# Illustrative numerical version using reasonable assumptions
# Part (a) matches the book exactly.
# Parts (b) and (c) use illustrative service-table and annuity assumptions.

library(mqriskR)

# ------------------------------------------------------------
# Part (a): accrued benefit and benefit accrual
# ------------------------------------------------------------

AB_56 <- 0.01 * 5 * 150000
AB_57 <- 0.01 * 6 * 156000
BA <- AB_57 - AB_56

out_part_a <- c(
  AB_56 = AB_56,
  AB_57 = AB_57,
  BA = BA
)

print(out_part_a)

# ------------------------------------------------------------
# Exercise 18-5 setup carried into 18-6
# ------------------------------------------------------------

x <- 51
z <- 65
i <- 0.06
v <- 1 / (1 + i)

# Decrement ages
y_er <- 61:64   # early retirement
y_w  <- 56:60   # withdrawal
y_i  <- 56:64   # disability
y_d  <- 61:64   # death with spouse annuity

# ------------------------------------------------------------
# Reasonable illustrative assumptions
# ------------------------------------------------------------

# Active-service one-year survival probabilities from age 51 onward
# (reflecting small combined decrement rates)
p_active_1yr <- c(
  "51" = 0.975,
  "52" = 0.974,
  "53" = 0.973,
  "54" = 0.972,
  "55" = 0.971,
  "56" = 0.970,
  "57" = 0.969,
  "58" = 0.968,
  "59" = 0.967,
  "60" = 0.966,
  "61" = 0.965,
  "62" = 0.964,
  "63" = 0.963,
  "64" = 0.962
)

# Build {}_{y-51}p_{51}^{(tau)} for y = 52,...,65
p_tau_to_y <- numeric(14)
names(p_tau_to_y) <- 52:65

cum_surv <- 1
for (age in 51:64) {
  cum_surv <- cum_surv * p_active_1yr[as.character(age)]
  p_tau_to_y[as.character(age + 1)] <- cum_surv
}

# Retirement decrements q_y^(r), ages 61:64
q_r <- c("61" = 0.08, "62" = 0.12, "63" = 0.18, "64" = 0.25)

# Withdrawal decrements q_y^(w), ages 56:60
q_w <- c("56" = 0.04, "57" = 0.035, "58" = 0.03, "59" = 0.025, "60" = 0.02)

# Disability decrements q_y^(i), ages 56:64
q_i <- c(
  "56" = 0.010, "57" = 0.010, "58" = 0.009, "59" = 0.009, "60" = 0.008,
  "61" = 0.008, "62" = 0.007, "63" = 0.007, "64" = 0.006
)

# Death decrements q_y^(d), ages 61:64
q_d <- c("61" = 0.003, "62" = 0.0035, "63" = 0.004, "64" = 0.0045)

# Survival from withdrawal age y+1/2 to 65, using a simple 98.5% annual survival
p_w_to_65 <- sapply(y_w, function(y) 0.985^(65 - y - 0.5))
names(p_w_to_65) <- y_w

# Retirement annuity-due values
adue_r_65 <- 12.0
adue_r_half <- c("61" = 12.8, "62" = 12.6, "63" = 12.4, "64" = 12.2)

# Disability annuity-due values, somewhat higher than retirement annuity
adue_i_half <- c(
  "56" = 14.5, "57" = 14.3, "58" = 14.1, "59" = 13.9, "60" = 13.7,
  "61" = 13.5, "62" = 13.3, "63" = 13.1, "64" = 12.9
)

# Spouse annuity-due values, spouse assumed younger
adue_spouse_half <- c("61" = 14.0, "62" = 13.8, "63" = 13.6, "64" = 13.4)

# ------------------------------------------------------------
# Helper for {}_{y-51}p_{51}^{(tau)}
# ------------------------------------------------------------
p_tau_from_51 <- function(y) {
  if (y == 51) return(1)
  p_tau_to_y[as.character(y)]
}

# ------------------------------------------------------------
# APV component calculator
# benefit_amount = 1860 for normal cost
# benefit_amount = 7500 for accrued liability
# ------------------------------------------------------------
apv_components_18_6 <- function(benefit_amount) {

  # (a) normal retirement
  APV_NR <- benefit_amount *
    v^(65 - 51) *
    p_tau_from_51(65) *
    adue_r_65

  # (b) early retirement
  APV_ER <- sum(sapply(y_er, function(y) {
    reduction <- 1 - 0.03 * (65 - y - 0.5)
    reduction *
      benefit_amount *
      v^(y + 0.5 - 51) *
      p_tau_from_51(y) *
      q_r[as.character(y)] *
      adue_r_half[as.character(y)]
  }))

  # (c) withdrawal
  APV_W <- sum(sapply(y_w, function(y) {
    benefit_amount *
      v^(65 - 51) *
      p_tau_from_51(y) *
      q_w[as.character(y)] *
      p_w_to_65[as.character(y)] *
      adue_r_65
  }))

  # (d) disability
  APV_I <- sum(sapply(y_i, function(y) {
    benefit_amount *
      v^(y + 0.5 - 51) *
      p_tau_from_51(y) *
      q_i[as.character(y)] *
      adue_i_half[as.character(y)]
  }))

  # (e) death
  APV_D <- sum(sapply(y_d, function(y) {
    0.50 *
      (1 - 0.03 * (65 - y - 0.5)) *
      benefit_amount *
      v^(y + 0.5 - 51) *
      p_tau_from_51(y) *
      q_d[as.character(y)] *
      adue_spouse_half[as.character(y)]
  }))

  c(
    APV_NR = APV_NR,
    APV_ER = APV_ER,
    APV_W  = APV_W,
    APV_I  = APV_I,
    APV_D  = APV_D,
    TOTAL  = APV_NR + APV_ER + APV_W + APV_I + APV_D
  )
}

# ------------------------------------------------------------
# Part (b): TUC normal cost using benefit accrual BA = 1860
# ------------------------------------------------------------
NC_components <- apv_components_18_6(benefit_amount = BA)

# ------------------------------------------------------------
# Part (c): TUC accrued liability using accrued benefit AB_56 = 7500
# ------------------------------------------------------------
AAL_components <- apv_components_18_6(benefit_amount = AB_56)

cat("\nPart (b): Illustrative TUC normal cost components\n")
print(round(NC_components, 2))

cat("\nPart (c): Illustrative TUC accrued liability components\n")
print(round(AAL_components, 2))

# ------------------------------------------------------------
# Plot component comparison
# ------------------------------------------------------------
comp_mat <- rbind(
  Normal_Cost = NC_components[1:5],
  Accrued_Liability = AAL_components[1:5]
)

barplot(
  t(comp_mat),
  beside = TRUE,
  legend.text = TRUE,
  main = "Exercise 18-6: Illustrative component APVs",
  ylab = "Actuarial present value"
)

