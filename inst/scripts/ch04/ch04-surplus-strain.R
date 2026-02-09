# Chapter 4 — Surplus strain intuition (Section 4.2.3)

# Toy model:
# - first-year commission/expense load is high
# - renewal expenses are lower
# This is NOT a reserving model; it's just to illustrate capital strain.

new_policies_per_year <- 500
years <- 10

premium <- 500
first_year_acq_cost <- 800   # commission + underwriting + issue costs
renewal_cost <- 50

capital0 <- 1e6
capital <- numeric(years + 1)
capital[1] <- capital0

inforce <- 0

for (y in 1:years) {
  newbiz <- new_policies_per_year
  # cash flow this year:
  # premium from inforce + newbiz, minus renewal costs on inforce, minus acquisition costs on newbiz
  cash_in <- premium * (inforce + newbiz)
  cash_out <- renewal_cost * inforce + first_year_acq_cost * newbiz

  capital[y + 1] <- capital[y] + (cash_in - cash_out)

  # next year's inforce (ignore lapses/mortality here—just a demo)
  inforce <- inforce + newbiz
}

data.frame(year = 0:years, capital = capital)
