library(testthat)
library(mqriskR)

# ------------------------------------------------------------
# Test life table
# ------------------------------------------------------------
tbl <- life_table(
  x = 0:4,
  px = c(0.90, 0.80, 0.60, 0.30, 0.00),
  radix = 10000
)

# lx:
# 0: 10000
# 1:  9000
# 2:  7200
# 3:  4320
# 4:  1296
# 5:     0

# ------------------------------------------------------------
# ex_curtate_tab
# ------------------------------------------------------------

test_that("ex_curtate_tab computes curtate expectation correctly", {

  # e_2 = {}_1 p_2 + {}_2 p_2 + {}_3 p_2
  #     = 0.60 + 0.18 + 0
  expect_equal(
    ex_curtate_tab(tbl, 2),
    0.60 + 0.18 + 0
  )

  # e_4 = 0, because l_5 = 0
  expect_equal(
    ex_curtate_tab(tbl, 4),
    0
  )
})

# ------------------------------------------------------------
# ex_temp_curtate_tab
# ------------------------------------------------------------

test_that("ex_temp_curtate_tab computes temporary curtate expectation correctly", {

  # e_{1:3} = {}_1 p_1 + {}_2 p_1 + {}_3 p_1
  #         = 0.80 + 0.48 + 0.144
  expect_equal(
    ex_temp_curtate_tab(tbl, x = 1, n = 3),
    0.80 + 0.48 + 0.144
  )

  # zero temporary period
  expect_equal(
    ex_temp_curtate_tab(tbl, x = 2, n = 0),
    0
  )
})

# ------------------------------------------------------------
# ex_temp_complete_tab under UDD
# ------------------------------------------------------------

test_that("ex_temp_complete_tab under UDD matches e_x + 1/2 within one year", {

  # For one year under UDD:
  # \int_0^1 {}_t p_x dt = 1 - q_x/2
  # For x = 2, q2 = 0.40 => 1 - 0.20 = 0.80
  expect_equal(
    ex_temp_complete_tab(tbl, x = 2, n = 1, assumption = "udd"),
    1 - 0.40 / 2
  )
})



test_that("ex_complete_tab under UDD matches year-by-year decomposition", {

  # For x = 2:
  # year (2,3]: 1 - q2/2 = 1 - 0.4/2 = 0.80
  # year (3,4]: p2 * (1 - q3/2) = 0.6 * 0.65 = 0.39
  # year (4,5]: p2*p3 * (1 - q4/2) = 0.18 * 0.50 = 0.09
  # total = 1.28
  expect_equal(
    ex_complete_tab(tbl, x = 2, assumption = "udd"),
    0.80 + 0.39 + 0.09
  )
})


# ------------------------------------------------------------
# ex_temp_complete_tab under constant force
# ------------------------------------------------------------

test_that("ex_temp_complete_tab under constant force matches closed form over one year", {

  # For one year under CF:
  # integral_0^1 p_x^t dt = (1 - p_x)/(-log p_x)
  px2 <- 0.60
  expected <- (1 - px2) / (-log(px2))

  expect_equal(
    ex_temp_complete_tab(tbl, x = 2, n = 1, assumption = "cf"),
    expected,
    tolerance = 1e-8
  )
})

# ------------------------------------------------------------
# ex_temp_complete_tab under Balducci
# ------------------------------------------------------------

test_that("ex_temp_complete_tab under Balducci is positive and less than 1 over one year", {

  val <- ex_temp_complete_tab(tbl, x = 2, n = 1, assumption = "balducci")

  expect_true(val > 0)
  expect_true(val < 1)
})

# ------------------------------------------------------------
# temporary complete expectation with fractional n
# ------------------------------------------------------------

test_that("ex_temp_complete_tab handles fractional n", {

  val <- ex_temp_complete_tab(tbl, x = 2, n = 1.5, assumption = "udd")

  # should be bigger than one-year value, less than full complete value
  val1 <- ex_temp_complete_tab(tbl, x = 2, n = 1, assumption = "udd")
  valfull <- ex_complete_tab(tbl, x = 2, assumption = "udd")

  expect_true(val > val1)
  expect_true(val < valfull)
})

# ------------------------------------------------------------
# vectorization
# ------------------------------------------------------------

test_that("expectation functions vectorize correctly", {

  vals <- ex_curtate_tab(tbl, c(1, 2, 3))
  expect_equal(length(vals), 3)

  vals2 <- ex_temp_curtate_tab(tbl, x = c(1, 2), n = c(2, 3))
  expect_equal(length(vals2), 2)

  vals3 <- ex_complete_tab(tbl, x = c(1, 2), assumption = "udd")
  expect_equal(length(vals3), 2)

  vals4 <- ex_temp_complete_tab(tbl, x = c(1, 2), n = c(1.5, 2), assumption = "cf")
  expect_equal(length(vals4), 2)
})

# ------------------------------------------------------------
# invalid inputs
# ------------------------------------------------------------

test_that("temporary curtate expectation rejects non-integer n", {

  expect_error(
    ex_temp_curtate_tab(tbl, x = 2, n = 1.5)
  )
})

test_that("temporary complete expectation rejects invalid n", {

  expect_error(
    ex_temp_complete_tab(tbl, x = 2, n = -1, assumption = "udd")
  )
})

test_that("complete expectation rejects invalid assumption", {

  expect_error(
    ex_complete_tab(tbl, x = 2, assumption = "wrong")
  )
})
