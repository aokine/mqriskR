library(testthat)
library(mqriskR)

# -------------------------------------------------------------------
# Build a small deterministic life table used in all tests
# -------------------------------------------------------------------

tbl <- life_table(
  x = 0:4,
  px = c(0.90, 0.80, 0.60, 0.30, 0.00),
  radix = 10000
)

# Expected lx values
lx_expected <- c(
  10000,
  9000,
  7200,
  4320,
  1296,
  0
)

# -------------------------------------------------------------------
# Constructor tests
# -------------------------------------------------------------------

test_that("life_table constructor builds correct lx values", {

  expect_s3_class(tbl, "life_table")

  expect_equal(tbl$lx, lx_expected)

})

# -------------------------------------------------------------------
# lx accessor
# -------------------------------------------------------------------

test_that("lx() returns correct l_x values", {

  expect_equal(lx(tbl, 0), 10000)
  expect_equal(lx(tbl, 2), 7200)
  expect_equal(lx(tbl, 5), 0)

})

# -------------------------------------------------------------------
# dx calculation
# -------------------------------------------------------------------

test_that("dx() computes deaths correctly", {

  expect_equal(dx(tbl, 0), 1000)
  expect_equal(dx(tbl, 2), 2880)
  expect_equal(dx(tbl, 4), 1296)

})

# -------------------------------------------------------------------
# ndx calculation
# -------------------------------------------------------------------

test_that("ndx() computes n-year deaths correctly", {

  # deaths between age 1 and age 4
  expect_equal(ndx(tbl, 1, 3), 7704)

})

# -------------------------------------------------------------------
# qx calculation
# -------------------------------------------------------------------

test_that("qx_tab() returns correct mortality rates", {

  expect_equal(qx_tab(tbl, 0), 0.10)
  expect_equal(qx_tab(tbl, 1), 0.20)
  expect_equal(qx_tab(tbl, 4), 1)

})

# -------------------------------------------------------------------
# npx calculation
# -------------------------------------------------------------------

test_that("npx() computes survival probabilities", {

  # 3-year survival from age 1
  expect_equal(npx(tbl, 1, 3), 1296 / 9000)

  # survival to terminal age
  expect_equal(npx(tbl, 2, 3), 0)

})

# -------------------------------------------------------------------
# nqx calculation
# -------------------------------------------------------------------

test_that("nqx() computes death probabilities", {

  expect_equal(nqx(tbl, 1, 3), 1 - (1296 / 9000))

  # guaranteed death before age 5
  expect_equal(nqx(tbl, 2, 3), 1)

})

# -------------------------------------------------------------------
# Internal identities
# -------------------------------------------------------------------

test_that("life table identities hold", {

  # qx = dx / lx
  expect_equal(qx_tab(tbl, 2), dx(tbl,2)/lx(tbl,2))

  # npx = lx+n / lx
  expect_equal(npx(tbl,1,3), lx(tbl,4)/lx(tbl,1))

  # nqx = 1 - npx
  expect_equal(nqx(tbl,1,3), 1 - npx(tbl,1,3))

})


testthat::test_that("closed life tables return zero lx beyond terminal age", {
  tbl <- life_table(x = 0:3, lx = c(100000, 80000, 50000, 0))

  testthat::expect_equal(lx(tbl, 4), 0)
  testthat::expect_equal(lx(tbl, 10), 0)
  testthat::expect_equal(npx(tbl, x = 2, n = 1), 0)
  testthat::expect_equal(npx(tbl, x = 2, n = 5), 0)
})

testthat::test_that("open life tables return NA beyond final age", {
  tbl <- life_table(x = 0:2, lx = c(100000, 80000, 50000))

  testthat::expect_true(is.na(lx(tbl, 3)))
  testthat::expect_true(is.na(npx(tbl, x = 1, n = 5)))
})

testthat::test_that("life table probabilities behave correctly near terminal age", {
  tbl <- life_table(x = 0:3, lx = c(100000, 80000, 50000, 0))

  testthat::expect_equal(dx(tbl, 2), 50000)
  testthat::expect_equal(qx_tab(tbl, 2), 1)
  testthat::expect_equal(npx(tbl, 2, 1), 0)
  testthat::expect_equal(nqx(tbl, 2, 1), 1)
})

testthat::test_that("life table functions vectorize over x and n", {
  tbl <- life_table(x = 0:4, lx = c(100000, 90000, 80000, 60000, 0))

  testthat::expect_equal(lx(tbl, c(0, 1, 2)), c(100000, 90000, 80000))

  testthat::expect_equal(
    npx(tbl, x = c(0, 1), n = c(1, 2)),
    c(90000 / 100000, 60000 / 90000)
  )

  testthat::expect_equal(
    nqx(tbl, x = c(0, 1), n = c(1, 2)),
    1 - c(90000 / 100000, 60000 / 90000)
  )
})
