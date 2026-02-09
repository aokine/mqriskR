# Chapter 4 — R solution support (selected exercises)
# Exercises: 4-4 and 4-8
#
# Goal:
# - Simulate/compute the one-year claim risk for a small term insurer.
# - Explore "How many policies should we sell?" and "How much capital is needed?"
#
# Model (simple, explicit assumptions):
# - One-year renewable term, benefit = 1000
# - Premium per policy = 50
# - Initial paid-in capital (given by exercise)
# - Number of insured lives = n
# - Each life has independent probability q of dying within 1 year
# - Claim count ~ Binomial(n, q), aggregate claims = 1000 * ClaimCount
#
# NOTE: Exercise 4-4 and 4-8 do NOT specify q.
# So we treat q as a parameter and do sensitivity analysis.
#
# --- Helpers ---------------------------------------------------------------

simulate_one_year <- function(n, q, benefit = 1000, prem = 50, capital = 2000,
                              n_sims = 200000, seed = 1) {
  set.seed(seed)
  claims_count <- rbinom(n_sims, size = n, prob = q)
  claims <- benefit * claims_count
  premium_income <- prem * n

  ending_surplus <- capital + premium_income - claims
  ruined <- ending_surplus < 0

  list(
    n = n, q = q, benefit = benefit, prem = prem, capital = capital,
    prob_ruin = mean(ruined),
    expected_claims = benefit * n * q,
    expected_profit = capital + premium_income - benefit * n * q,
    ending_surplus = ending_surplus
  )
}

# Exact (non-simulation) probability of ruin under Binomial claims:
# Ruin occurs if capital + prem*n - benefit*X < 0  =>  X > (capital + prem*n)/benefit
prob_ruin_exact <- function(n, q, benefit = 1000, prem = 50, capital = 2000) {
  threshold <- (capital + prem * n) / benefit
  k0 <- floor(threshold)  # ruin if X >= k0 + 1
  if (k0 >= n) return(0)
  if (k0 < 0) return(1)
  pbinom(k0, size = n, prob = q, lower.tail = FALSE)  # P[X > k0]
}

# Capital required for a target ruin probability alpha:
# Find smallest capital such that P(ruin) <= alpha.
required_capital_for_ruin <- function(n, q, alpha = 0.01, benefit = 1000, prem = 50,
                                      cap_grid = seq(0, 200000, by = 100)) {
  pr <- vapply(cap_grid, function(c0) prob_ruin_exact(n, q, benefit, prem, c0), numeric(1))
  idx <- which(pr <= alpha)
  if (length(idx) == 0) return(NA_real_)
  cap_grid[min(idx)]
}

# Risk measures (VaR / TVaR) of aggregate claims, for a given n,q:
claims_risk_measures <- function(n, q, p = 0.99, benefit = 1000) {
  # Claim count distribution is Binomial(n,q)
  # VaR of claims is benefit * quantile of Binomial
  k_var <- qbinom(p, size = n, prob = q)
  var_claims <- benefit * k_var

  # TVaR (conditional tail mean) for discrete distribution:
  ks <- 0:n
  probs <- dbinom(ks, size = n, prob = q)
  claims <- benefit * ks

  tail <- claims >= var_claims
  tvar_claims <- sum(claims[tail] * probs[tail]) / sum(probs[tail])

  list(VaR = var_claims, TVaR = tvar_claims, k_var = k_var)
}

# --- Exercise 4-4 ----------------------------------------------------------
# Risk Mitigation, Inc.:
# capital = 2000, benefit = 1000, prem = 50
# Question asks how many policies to sell and how much capital to hold.
#
# We'll explore ruin probabilities for a range of n under plausible q values.

ex4_4_analysis <- function(q_values = c(0.001, 0.002, 0.005, 0.01),
                           n_values = c(10, 20, 30, 40, 50, 75, 100, 150, 200),
                           capital = 2000, benefit = 1000, prem = 50,
                           alpha_targets = c(0.10, 0.05, 0.01)) {

  out <- list()

  for (q in q_values) {
    tab <- data.frame(
      n = n_values,
      q = q,
      prob_ruin = NA_real_,
      expected_claims = NA_real_,
      expected_profit_ex_capital = NA_real_  # premium - expected claims
    )

    for (i in seq_along(n_values)) {
      n <- n_values[i]
      tab$prob_ruin[i] <- prob_ruin_exact(n, q, benefit, prem, capital)
      tab$expected_claims[i] <- benefit * n * q
      tab$expected_profit_ex_capital[i] <- prem * n - benefit * n * q
    }

    # Capital requirement for fixed n and target ruin probs
    cap_req <- expand.grid(n = n_values, alpha = alpha_targets)
    cap_req$q <- q
    cap_req$cap_required <- mapply(
      function(n, alpha) required_capital_for_ruin(n, q, alpha, benefit, prem),
      cap_req$n, cap_req$alpha
    )

    out[[paste0("q=", q)]] <- list(table = tab, capital_requirements = cap_req)
  }

  out
}

# Run and print a compact view:
res_4_4 <- ex4_4_analysis()

# Print a quick summary for one q (you can change which one):
print(res_4_4[["q=0.005"]]$table)
print(res_4_4[["q=0.005"]]$capital_requirements)

# --- Exercise 4-8 ----------------------------------------------------------
# New Beginnings, Inc.:
# capital = 3000, and a client wants 1,000 one-year policies.
#
# We quantify the concentration risk via ruin probability and capital needed.

ex4_8_analysis <- function(q = 0.005, n = 1000, capital = 3000,
                           benefit = 1000, prem = 50,
                           alpha_targets = c(0.10, 0.05, 0.01, 0.001)) {

  pr <- prob_ruin_exact(n, q, benefit, prem, capital)
  rm <- claims_risk_measures(n, q, p = 0.99, benefit = benefit)

  cap_req <- data.frame(
    alpha = alpha_targets,
    cap_required = vapply(alpha_targets, function(a)
      required_capital_for_ruin(n, q, a, benefit, prem), numeric(1)
    )
  )

  list(
    inputs = list(n = n, q = q, capital = capital, benefit = benefit, prem = prem),
    prob_ruin = pr,
    claims_VaR_99 = rm$VaR,
    claims_TVaR_99 = rm$TVaR,
    capital_requirements = cap_req
  )
}

res_4_8 <- ex4_8_analysis(q = 0.005, n = 1000, capital = 3000)
print(res_4_8)

# Optional: Monte Carlo sanity check for 4-8 (can be slow at large n_sims)
# sim_4_8 <- simulate_one_year(n = 1000, q = 0.005, capital = 3000, n_sims = 200000)
# sim_4_8$prob_ruin

