# Chapter 4 — Risk pooling and capital stress (Exercises 4-4, 4-8)

set.seed(1)

simulate_one_year <- function(n_policies, p_claim, benefit, premium, capital, n_sims = 10000) {
  claims <- rbinom(n_sims, size = n_policies, prob = p_claim)
  profit <- n_policies * premium - claims * benefit
  end_capital <- capital + profit

  list(
    ruin_prob = mean(end_capital < 0),
    mean_end_capital = mean(end_capital),
    q05_end_capital = unname(quantile(end_capital, 0.05)),
    q01_end_capital = unname(quantile(end_capital, 0.01))
  )
}

# Example parameters (change these for the exercises)
benefit <- 1000
premium <- 50
capital <- 2000

# Try a range of portfolio sizes and claim probabilities
n_grid <- c(10, 20, 50, 100, 200, 500, 1000)
p_claim <- 0.02

res <- lapply(n_grid, function(n) {
  out <- simulate_one_year(n, p_claim, benefit, premium, capital)
  c(n_policies = n, out)
})

print(do.call(rbind, res))


