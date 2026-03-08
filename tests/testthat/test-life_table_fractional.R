library(testthat)
library(mqriskR)

# ------------------------------------------------------------
# Test table used throughout
# ------------------------------------------------------------
tbl <- life_table(
  x = 0:4,
  px = c(0.90, 0.80, 0.60, 0.30, 0.00),
  radix = 10000
)

# For x = 2:
# l2 = 7200, l3 = 4320
# p2 = 0.60, q2 = 0.40

# ------------------------------------------------------------
# tpx_tab
# ------------------------------------------------------------

test_that("tpx_tab matches Chapter 6 formulas under all assumptions", {

  expect_equal(
    tpx_tab(tbl, x = 2, t = 0.25, assumption = "udd"),
    1 - 0.25 * 0.40
  )

  expect_equal(
    tpx_tab(tbl, x = 2, t = 0.25, assumption = "cf"),
    0.60^0.25
  )

  expect_equal(
    tpx_tab(tbl, x = 2, t = 0.25, assumption = "balducci"),
    (1 - 0.40) / (1 - (1 - 0.25) * 0.40)
  )
})

# ------------------------------------------------------------
# tqx_tab
# ------------------------------------------------------------

test_that("tqx_tab equals 1 - tpx_tab under all assumptions", {

  udd_p <- tpx_tab(tbl, x = 2, t = 0.40, assumption = "udd")
  cf_p  <- tpx_tab(tbl, x = 2, t = 0.40, assumption = "cf")
  bal_p <- tpx_tab(tbl, x = 2, t = 0.40, assumption = "balducci")

  expect_equal(
    tqx_tab(tbl, x = 2, t = 0.40, assumption = "udd"),
    1 - udd_p
  )

  expect_equal(
    tqx_tab(tbl, x = 2, t = 0.40, assumption = "cf"),
    1 - cf_p
  )

  expect_equal(
    tqx_tab(tbl, x = 2, t = 0.40, assumption = "balducci"),
    1 - bal_p
  )
})

# ------------------------------------------------------------
# mux_tab
# ------------------------------------------------------------

test_that("mux_tab matches Chapter 6 formulas under all assumptions", {

  qx <- 0.40
  px <- 0.60
  t  <- 0.25

  expect_equal(
    mux_tab(tbl, x = 2, t = t, assumption = "udd"),
    qx / (1 - t * qx)
  )

  expect_equal(
    mux_tab(tbl, x = 2, t = t, assumption = "cf"),
    -log(px)
  )

  expect_equal(
    mux_tab(tbl, x = 2, t = t, assumption = "balducci"),
    qx / (1 - (1 - t) * qx)
  )
})

# ------------------------------------------------------------
# fx_tab
# ------------------------------------------------------------

test_that("fx_tab equals tpx_tab * mux_tab", {

  t <- 0.25

  expect_equal(
    fx_tab(tbl, x = 2, t = t, assumption = "udd"),
    tpx_tab(tbl, x = 2, t = t, assumption = "udd") *
      mux_tab(tbl, x = 2, t = t, assumption = "udd")
  )

  expect_equal(
    fx_tab(tbl, x = 2, t = t, assumption = "cf"),
    tpx_tab(tbl, x = 2, t = t, assumption = "cf") *
      mux_tab(tbl, x = 2, t = t, assumption = "cf")
  )

  expect_equal(
    fx_tab(tbl, x = 2, t = t, assumption = "balducci"),
    tpx_tab(tbl, x = 2, t = t, assumption = "balducci") *
      mux_tab(tbl, x = 2, t = t, assumption = "balducci")
  )
})

# ------------------------------------------------------------
# boundary values
# ------------------------------------------------------------

test_that("fractional probabilities have correct boundary values", {

  expect_equal(tpx_tab(tbl, x = 2, t = 0, assumption = "udd"), 1)
  expect_equal(tpx_tab(tbl, x = 2, t = 0, assumption = "cf"), 1)
  expect_equal(tpx_tab(tbl, x = 2, t = 0, assumption = "balducci"), 1)

  expect_equal(tqx_tab(tbl, x = 2, t = 0, assumption = "udd"), 0)
  expect_equal(tqx_tab(tbl, x = 2, t = 0, assumption = "cf"), 0)
  expect_equal(tqx_tab(tbl, x = 2, t = 0, assumption = "balducci"), 0)

  # At t = 1 these should equal p_x and q_x
  expect_equal(tpx_tab(tbl, x = 2, t = 1, assumption = "udd"), 0.60)
  expect_equal(tpx_tab(tbl, x = 2, t = 1, assumption = "cf"), 0.60)
  expect_equal(tpx_tab(tbl, x = 2, t = 1, assumption = "balducci"), 0.60)

  expect_equal(tqx_tab(tbl, x = 2, t = 1, assumption = "udd"), 0.40)
  expect_equal(tqx_tab(tbl, x = 2, t = 1, assumption = "cf"), 0.40)
  expect_equal(tqx_tab(tbl, x = 2, t = 1, assumption = "balducci"), 0.40)
})

# ------------------------------------------------------------
# nmxq
# ------------------------------------------------------------

test_that("nmxq computes deferred death probabilities correctly", {

  # {}_{2|1} q_1 = {}_2 p_1 * q_3
  # p1 = 0.80, p2 = 0.60 so {}_2 p_1 = 0.48
  # q3 = 0.70
  expect_equal(
    nmxq(tbl, x = 1, n = 2, m = 1),
    0.48 * 0.70
  )
})

# ------------------------------------------------------------
# nkqx
# ------------------------------------------------------------

test_that("nkqx computes curtate death probabilities correctly", {

  # {}_{2|} q_1 = {}_2 p_1 - {}_3 p_1
  # {}_2 p_1 = 0.80 * 0.60 = 0.48
  # {}_3 p_1 = 0.80 * 0.60 * 0.30 = 0.144
  expect_equal(
    nkqx(tbl, x = 1, k = 2),
    0.48 - 0.144
  )
})

# ------------------------------------------------------------
# vectorization
# ------------------------------------------------------------

test_that("fractional functions vectorize correctly", {

  tvals <- c(0.25, 0.50, 0.75)

  expect_equal(
    length(tpx_tab(tbl, x = 2, t = tvals, assumption = "udd")),
    3
  )

  expect_equal(
    length(mux_tab(tbl, x = 2, t = tvals, assumption = "cf")),
    3
  )

  expect_equal(
    length(fx_tab(tbl, x = c(1, 2, 3), t = 0.25, assumption = "balducci")),
    3
  )
})

# ------------------------------------------------------------
# invalid inputs
# ------------------------------------------------------------

test_that("fractional functions reject invalid t", {

  expect_error(
    tpx_tab(tbl, x = 2, t = -0.1, assumption = "udd")
  )

  expect_error(
    tpx_tab(tbl, x = 2, t = 1.2, assumption = "cf")
  )
})

test_that("fractional functions reject invalid assumption", {

  expect_error(
    tpx_tab(tbl, x = 2, t = 0.5, assumption = "wrong")
  )
})
