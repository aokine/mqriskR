test_that("interest conversions are consistent", {
  out <- interest_convert(i = 0.05)
  expect_equal(out$d, 0.05/1.05, tolerance = 1e-12)
  expect_equal(out$delta, log(1.05), tolerance = 1e-12)
})

test_that("annuity due relationship holds (annual)", {
  a  <- annuity_certain(n = 10, i = 0.05, due = FALSE)
  ad <- annuity_certain(n = 10, i = 0.05, due = TRUE)
  expect_equal(ad, (1 + 0.05) * a, tolerance = 1e-10)
})

test_that("pv_cashflows matches manual discounting", {
  cf <- c(-100, 60, 60)
  t <- c(0, 1, 2)
  i <- 0.10
  expect_equal(pv_cashflows(cf, t, i), sum(cf * (1+i)^(-t)), tolerance = 1e-12)
})

test_that("solve_yield finds a root", {
  cf <- c(-100, 60, 60)
  t <- c(0, 1, 2)
  i_hat <- solve_yield(cf, t, interval = c(-0.5, 1))
  expect_lt(abs(pv_cashflows(cf, t, i_hat)), 1e-8)
})
