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
