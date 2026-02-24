## Listing 4.1: Risk pooling and one-year capital stress (Exercises 4-4 and 4-8)
set.seed(1)
simulate_one_year <- function(n_policies, p_claim, benefit, premium, capital, n_sims = 10000) {
  claims <- rbinom(n_sims, size = n_policies, prob = p_claim)
  profit <- n_policies * premium - claims * benefit
  end_capital <- capital + profit
  list(
    ruin_prob = mean(end_capital < 0),
    mean_end_capital = mean(end_capital),
    q05_end_capital = unname(quantile(end_capital, 0.05))
  )
}
benefit <- 1000
premium <- 50
capital <- 2000
p_claim <- 0.02

n_grid <- c(10, 20, 50, 100, 200, 500, 1000)

res <- lapply(n_grid, function(n) {
  out <- simulate_one_year(n, p_claim, benefit, premium, capital)
  c(n_policies = n, out)
})

print(do.call(rbind, res))


## Listing 4.2: Surplus strain intuition from high first-year acquisition costs
new_policies_per_year <- 500
years <- 10

premium <- 500
first_year_acq_cost <- 800
renewal_cost <- 50

capital0 <- 1e6
capital <- numeric(years + 1)
capital[1] <- capital0

inforce <- 0

for (y in 1:years) {
  newbiz <- new_policies_per_year
  cash_in <- premium * (inforce + newbiz)
  cash_out <- renewal_cost * inforce + first_year_acq_cost * newbiz
  capital[y + 1] <- capital[y] + (cash_in - cash_out)
  inforce <- inforce + newbiz
}

data.frame(year = 0:years, capital = capital)


## Listing 4.3: Defined benefit pension formula illustration
db_pension <- function(final_avg_salary, years_service, accrual = 1 + 2/3) {
  rate <- accrual / 100
  rate * final_avg_salary * years_service
}

db_pension(final_avg_salary = 90000, years_service = 25)



