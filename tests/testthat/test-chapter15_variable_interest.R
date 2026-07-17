test_that("nEx_var reproduces Example 15.1", {
  qx <- c(.03, .04, .05, .06, .07)

  val1 <- nEx_var(qx = qx, i = c(.06, .07, .08, .09, .10), benefit = 1000)
  val2 <- nEx_var(qx = qx, i = c(.06, .06, .06, .06, .06), benefit = 1000)
  val3 <- nEx_var(qx = qx, i = c(.06, .05, .04, .03, .02), benefit = 1000)

  expect_equal(val1, 526.61, tolerance = 1e-2)
  expect_equal(val2, 577.93, tolerance = 1e-2)
  expect_equal(val3, 635.97, tolerance = 1e-2)
})

test_that("Axn1_var reproduces Example 15.2", {
  qx <- c(.03, .04, .05, .06, .07)

  val1 <- Axn1_var(qx = qx, i = c(.06, .07, .08, .09, .10))
  val2 <- Axn1_var(qx = qx, i = c(.06, .06, .06, .06, .06))
  val3 <- Axn1_var(qx = qx, i = c(.06, .05, .04, .03, .02))

  expect_equal(val1, 0.1799, tolerance = 1e-4)
  expect_equal(val2, 0.1875, tolerance = 1e-3)
  expect_equal(val3, 0.1958, tolerance = 1e-3)
})

test_that("Axn_var equals Axn1_var plus nEx_var", {
  qx <- c(.03, .04, .05, .06, .07)
  i  <- c(.06, .07, .08, .09, .10)

  term_val  <- Axn1_var(qx = qx, i = i)
  pe_val    <- nEx_var(qx = qx, i = i)
  endow_val <- Axn_var(qx = qx, i = i)

  expect_equal(endow_val, term_val + pe_val, tolerance = 1e-12)
})

test_that("axn_var reproduces Example 15.3 annuity values", {
  qx <- rep(.02, 5)

  level <- axn_var(qx = qx, i = c(.06, .06, .06, .06, .06), type = "immediate")
  fall  <- axn_var(qx = qx, i = c(.06, .05, .04, .03, .03), type = "immediate")
  rise  <- axn_var(qx = qx, i = c(.06, .07, .08, .09, .10), type = "immediate")

  expect_equal(level, 3.9756, tolerance = 1e-4)
  expect_equal(fall, 4.1102, tolerance = 1e-4)
  expect_equal(rise, 3.8461, tolerance = 1e-4)
})

test_that("z_from_coupon_semi reproduces Table 15.6 approximately", {
  maturity <- c(0.5, 1.0, 1.5, 2.0)
  coupon_yield <- c(0.0244, 0.0260, 0.0276, 0.0293)

  z <- z_from_coupon_semi(maturity, coupon_yield)

  expect_equal(100 * z[1], 2.440, tolerance = 1e-3)
  expect_equal(100 * z[2], 2.601, tolerance = 1e-3)
  expect_equal(100 * z[3], 2.763, tolerance = 1e-3)
  expect_equal(100 * z[4], 2.936, tolerance = 1e-3)
})

test_that("pv_spot_cashflows reproduces Example 15.4", {
  val <- pv_spot_cashflows(
    amounts = c(200000, 50000, 50000, 100000),
    times   = c(0, 0.5, 1, 2),
    spot    = c(0, 0.02440, 0.02601, 0.02936),
    compounding = "semiannual"
  )

  expect_equal(val, 392459.12, tolerance = 1e-2)
})

test_that("spot-rate Chapter 15 functions reproduce Example 15.5 premium pieces", {
  qx <- c(.02, .03, .04, .05, .06)
  z  <- c(.03, .04, .05, .06, .07)

  apvb <- Axn1_spot(qx = qx, z = z)
  apvp_factor <- axn_spot(qx = qx, z = z, type = "due")

  premium <- 1000000 * apvb / apvp_factor

  expect_equal(apvb, 0.1527, tolerance = 5e-4)
  expect_equal(apvp_factor, 4.3054, tolerance = 1e-4)
  expect_equal(premium, 35467.09, tolerance = 1e-2)
})




test_that("spot-rate identities hold", {
  qx <- c(.02, .03, .04, .05, .06)
  z  <- c(.03, .04, .05, .06, .07)

  expect_equal(
    Axn_spot(qx = qx, z = z),
    Axn1_spot(qx = qx, z = z) + nEx_spot(qx = qx, z = z),
    tolerance = 1e-12
  )
})

test_that("fnk_from_z reproduces Example 15.6 values", {
  z <- c(0.03, 0.04, 0.05, 0.06, 0.07)

  expect_equal(100 * fnk_from_z(z, n = 1, k = 4), 8.024, tolerance = 1e-3)
  expect_equal(100 * fnk_from_z(z, n = 2, k = 2), 8.038, tolerance = 1e-3)
})

test_that("forward_matrix_from_z gives expected selected entries", {
  z <- c(0.03, 0.04, 0.05, 0.06, 0.07)
  fmat <- forward_matrix_from_z(z)

  expect_equal(100 * fmat["n=1", "k=1"], 5.01, tolerance = 1e-2)
  expect_equal(100 * fmat["n=1", "k=2"], 6.01, tolerance = 1e-2)
  expect_equal(100 * fmat["n=2", "k=1"], 7.03, tolerance = 1e-2)
  expect_true(is.na(fmat["n=4", "k=2"]))
})

test_that("z_from_fn1 converts one-year forward rates to spot rates", {
  fn1 <- c(0.04, 0.05, 0.06, 0.07, 0.08)
  z <- z_from_fn1(fn1)

  expected <- c(
    0.0400000000,
    ((1.04 * 1.05)^(1 / 2)) - 1,
    ((1.04 * 1.05 * 1.06)^(1 / 3)) - 1,
    ((1.04 * 1.05 * 1.06 * 1.07)^(1 / 4)) - 1,
    ((1.04 * 1.05 * 1.06 * 1.07 * 1.08)^(1 / 5)) - 1
  )

  expect_equal(z, expected, tolerance = 1e-12)
})

test_that("z_from_coupon_annual works on a simple annual table", {
  maturity <- 1:4
  coupon_yield <- c(0.02, 0.04, 0.06, 0.08)

  z <- z_from_coupon_annual(maturity, coupon_yield)

  expect_equal(z[1], 0.02, tolerance = 1e-12)
  expect_true(all(diff(z) > 0))
})

test_that("vt_var returns cumulative variable-rate discount factors", {
  v <- vt_var(c(0.06, 0.07, 0.08))

  expect_equal(v[1], 1 / 1.06, tolerance = 1e-12)
  expect_equal(v[2], 1 / (1.06 * 1.07), tolerance = 1e-12)
  expect_equal(v[3], 1 / (1.06 * 1.07 * 1.08), tolerance = 1e-12)
})




testthat::test_that("variable discount factors allow valid negative rates", {
  rates <- c(-0.01, 0.02, 0.03)

  expected <- cumprod(1 / (1 + rates))

  testthat::expect_equal(
    vt_var(rates),
    expected,
    tolerance = 1e-12
  )
})


testthat::test_that("variable-interest APVs satisfy standard identities", {
  qx <- c(0.02, 0.03, 0.04, 0.05)
  rates <- c(0.03, 0.04, 0.05, 0.06)

  testthat::expect_equal(
    Axn_var(qx, rates),
    Axn1_var(qx, rates) + nEx_var(qx, rates),
    tolerance = 1e-12
  )
})


testthat::test_that("spot-rate APVs satisfy standard identities", {
  qx <- c(0.02, 0.03, 0.04, 0.05)
  spot <- c(0.03, 0.04, 0.05, 0.06)

  testthat::expect_equal(
    Axn_spot(qx, spot),
    Axn1_spot(qx, spot) + nEx_spot(qx, spot),
    tolerance = 1e-12
  )
})


testthat::test_that("annuity-due includes the time-zero payment", {
  qx <- rep(0.02, 5)
  rates <- rep(0.05, 5)

  immediate <- axn_var(qx, rates, type = "immediate")
  due <- axn_var(qx, rates, type = "due")

  testthat::expect_true(due > immediate)
  testthat::expect_equal(
    due,
    sum(
      c(1, head(vt_var(rates), -1)) *
        c(1, cumprod(1 - qx))[seq_along(qx)]
    ),
    tolerance = 1e-12
  )
})


testthat::test_that("constant variable rates agree with direct discounting", {
  qx <- c(0.02, 0.03, 0.04)
  rate <- 0.05
  rates <- rep(rate, length(qx))

  survival_start <- c(1, cumprod(1 - qx))[seq_along(qx)]
  expected <- sum(
    (1 + rate)^(-seq_along(qx)) *
      survival_start *
      qx
  )

  testthat::expect_equal(
    Axn1_var(qx, rates),
    expected,
    tolerance = 1e-12
  )
})


testthat::test_that("spot cash-flow valuation supports valid negative rates", {
  out <- pv_spot_cashflows(
    amounts = c(100, 100),
    times = c(1, 2),
    spot = c(-0.01, 0.02),
    compounding = "annual"
  )

  expected <- 100 / 0.99 + 100 / 1.02^2

  testthat::expect_equal(out, expected, tolerance = 1e-12)
})


testthat::test_that("annual spot bootstrapping reprices par bonds", {
  maturity <- 1:4
  coupon_yield <- c(0.02, 0.04, 0.06, 0.08)
  par <- 1000

  spot <- z_from_coupon_annual(
    maturity = maturity,
    coupon_yield = coupon_yield,
    par = par
  )

  for (j in seq_along(maturity)) {
    coupon <- par * coupon_yield[j]
    cash_flows <- c(
      rep(coupon, maturity[j] - 1L),
      par + coupon
    )

    price <- sum(
      cash_flows /
        (1 + spot[seq_len(maturity[j])])^seq_len(maturity[j])
    )

    testthat::expect_equal(price, par, tolerance = 1e-8)
  }
})


testthat::test_that("semiannual spot bootstrapping reprices par bonds", {
  maturity <- c(0.5, 1.0, 1.5, 2.0)
  coupon_yield <- c(0.0244, 0.0260, 0.0276, 0.0293)
  par <- 1000

  spot <- z_from_coupon_semi(
    maturity = maturity,
    coupon_yield = coupon_yield,
    par = par
  )

  periods <- seq_along(maturity)

  for (j in seq_along(maturity)) {
    coupon <- par * coupon_yield[j] / 2
    cash_flows <- c(
      rep(coupon, periods[j] - 1L),
      par + coupon
    )

    price <- sum(
      cash_flows /
        (1 + spot[seq_len(periods[j])] / 2)^seq_len(periods[j])
    )

    testthat::expect_equal(price, par, tolerance = 1e-8)
  }
})


testthat::test_that("forward and spot conversions are internally consistent", {
  forward <- c(0.04, 0.05, 0.06, 0.07, 0.08)
  spot <- z_from_fn1(forward)

  recovered <- c(
    spot[1],
    vapply(
      1:(length(spot) - 1L),
      function(n) fnk_from_z(spot, n = n, k = 1),
      numeric(1)
    )
  )

  testthat::expect_equal(
    recovered,
    forward,
    tolerance = 1e-12
  )
})


testthat::test_that("forward rate matrix has expected structure", {
  spot <- c(0.03, 0.04, 0.05, 0.06)

  out <- forward_matrix_from_z(spot)

  testthat::expect_true(is.matrix(out))
  testthat::expect_equal(dim(out), c(3, 3))
  testthat::expect_equal(out[1, 1], fnk_from_z(spot, 1, 1))
  testthat::expect_true(is.na(out[3, 2]))
})


testthat::test_that("variable-interest functions reject invalid inputs", {
  testthat::expect_error(
    vt_var(numeric(0)),
    "positive length"
  )

  testthat::expect_error(
    vt_var(c(0.05, -1)),
    "greater than -1"
  )

  testthat::expect_error(
    Axn1_var(
      qx = c(0.02, 0.03),
      i = c(0.04, 0.05, 0.06)
    ),
    "same length"
  )

  testthat::expect_error(
    nEx_var(
      qx = c(0.02, 0.03),
      i = c(0.04, 0.05),
      benefit = c(1, 2)
    ),
    "benefit must be a numeric scalar"
  )

  testthat::expect_error(
    forward_matrix_from_z(0.05),
    "at least two"
  )
})


testthat::test_that("bootstrapping requires consecutive maturities", {
  testthat::expect_error(
    z_from_coupon_annual(
      maturity = c(1, 3),
      coupon_yield = c(0.03, 0.04)
    ),
    "consecutive integer maturities"
  )

  testthat::expect_error(
    z_from_coupon_semi(
      maturity = c(0.5, 1.5),
      coupon_yield = c(0.03, 0.04)
    ),
    "consecutive half-year maturities"
  )
})
