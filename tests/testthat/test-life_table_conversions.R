library(testthat)
library(mqriskR)

# ------------------------------------------------------------
# S0_to_lx
# ------------------------------------------------------------

test_that("S0_to_lx converts survival probabilities correctly", {

  S0_vals <- c(1.00, 0.90, 0.72, 0.432)

  expect_equal(
    S0_to_lx(S0_vals, radix = 10000),
    c(10000, 9000, 7200, 4320)
  )
})

test_that("S0_to_lx rejects invalid inputs", {

  expect_error(S0_to_lx(c(1, 1.1, 0.8)))
  expect_error(S0_to_lx(c(1, -0.1, 0.8)))
  expect_error(S0_to_lx(c(1, 0.9, 0.95)))   # not nonincreasing
  expect_error(S0_to_lx(c(1, 0.9), radix = -1))
})

# ------------------------------------------------------------
# lx_to_S0
# ------------------------------------------------------------

test_that("lx_to_S0 converts life-table values correctly", {

  lx_vals <- c(10000, 9000, 7200, 4320)

  expect_equal(
    lx_to_S0(lx_vals),
    c(1.00, 0.90, 0.72, 0.432)
  )
})

test_that("lx_to_S0 rejects invalid inputs", {

  expect_error(lx_to_S0(c(10000, 9000, 9500)))  # not nonincreasing
  expect_error(lx_to_S0(c(10000, -1, 500)))
  expect_error(lx_to_S0(c(0, 0, 0)))
})

# ------------------------------------------------------------
# px_to_lx
# ------------------------------------------------------------

test_that("px_to_lx constructs life-table values correctly", {

  px_vals <- c(0.90, 0.80, 0.60, 0.30, 0.00)

  expect_equal(
    px_to_lx(px_vals, radix = 10000),
    c(10000, 9000, 7200, 4320, 1296, 0)
  )
})

test_that("px_to_lx rejects invalid inputs", {

  expect_error(px_to_lx(c(0.9, 1.2, 0.5)))
  expect_error(px_to_lx(c(0.9, -0.1, 0.5)))
  expect_error(px_to_lx(c(0.9, 0.8), radix = 0))
})

# ------------------------------------------------------------
# qx_to_lx
# ------------------------------------------------------------

test_that("qx_to_lx constructs life-table values correctly", {

  qx_vals <- c(0.10, 0.20, 0.40, 0.70, 1.00)

  expect_equal(
    qx_to_lx(qx_vals, radix = 10000),
    c(10000, 9000, 7200, 4320, 1296, 0)
  )
})

test_that("qx_to_lx rejects invalid inputs", {

  expect_error(qx_to_lx(c(0.1, 1.2, 0.5)))
  expect_error(qx_to_lx(c(0.1, -0.1, 0.5)))
  expect_error(qx_to_lx(c(0.1, 0.2), radix = -100))
})

# ------------------------------------------------------------
# consistency checks
# ------------------------------------------------------------

test_that("px_to_lx and qx_to_lx are consistent", {

  px_vals <- c(0.90, 0.80, 0.60, 0.30, 0.00)
  qx_vals <- 1 - px_vals

  expect_equal(
    px_to_lx(px_vals, radix = 10000),
    qx_to_lx(qx_vals, radix = 10000)
  )
})

test_that("lx_to_S0 and S0_to_lx are inverse up to radix", {

  lx_vals <- c(10000, 9000, 7200, 4320)
  S0_vals <- lx_to_S0(lx_vals)

  expect_equal(
    S0_to_lx(S0_vals, radix = 10000),
    lx_vals
  )
})
