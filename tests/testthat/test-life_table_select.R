library(testthat)
library(mqriskR)

# ------------------------------------------------------------
# Small select-and-ultimate style table for testing
# ------------------------------------------------------------
# 2-year select period for simplicity
#
# Row interpretation:
# x_sel = age at selection
# duration = years since selection
# attained_age = x_sel + duration
# lx = survivor count at that select age/duration
#
# For example:
# l_[20]     = 100000
# l_[20]+1   =  99000
# l_22       =  97000   (ultimate at attained age 22)

sel_tbl <- select_life_table(
  x_sel = c(20,20,20,
            21,21,21,
            22,22,22),
  duration = c(0,1,2,
               0,1,2,
               0,1,2),
  attained_age = c(20,21,22,
                   21,22,23,
                   22,23,24),
  lx = c(100000, 99000, 97000,
         99500, 97500, 95000,
         98000, 95500, 92000)
)

# ------------------------------------------------------------
# constructor
# ------------------------------------------------------------

test_that("select_life_table constructor creates valid object", {
  expect_s3_class(sel_tbl, "select_life_table")
  expect_equal(nrow(sel_tbl), 9)
  expect_true(all(c("x_sel", "duration", "attained_age", "lx") %in% names(sel_tbl)))
})

# ------------------------------------------------------------
# lx_select
# ------------------------------------------------------------

test_that("lx_select returns correct select-table survivor values", {
  expect_equal(lx_select(sel_tbl, 20, 0), 100000)
  expect_equal(lx_select(sel_tbl, 20, 1),  99000)
  expect_equal(lx_select(sel_tbl, 20, 2),  97000)

  expect_equal(lx_select(sel_tbl, 21, 0),  99500)
  expect_equal(lx_select(sel_tbl, 21, 2),  95000)
})

test_that("lx_select returns NA for missing combinations", {
  expect_true(is.na(lx_select(sel_tbl, 20, 3)))
})

# ------------------------------------------------------------
# npx_select
# ------------------------------------------------------------

test_that("npx_select computes select survival probabilities correctly", {
  # {}_1 p_[20] = l_[20]+1 / l_[20] = 99000 / 100000
  expect_equal(
    npx_select(sel_tbl, x_sel = 20, t = 0, n = 1),
    99000 / 100000
  )

  # {}_2 p_[20] = l_[20]+2 / l_[20] = 97000 / 100000
  expect_equal(
    npx_select(sel_tbl, x_sel = 20, t = 0, n = 2),
    97000 / 100000
  )

  # {}_1 p_[20]+1 = l_[20]+2 / l_[20]+1 = 97000 / 99000
  expect_equal(
    npx_select(sel_tbl, x_sel = 20, t = 1, n = 1),
    97000 / 99000
  )
})

test_that("npx_select returns NA when future select value is unavailable", {
  expect_true(is.na(npx_select(sel_tbl, x_sel = 20, t = 1, n = 2)))
})

# ------------------------------------------------------------
# nqx_select
# ------------------------------------------------------------

test_that("nqx_select computes select death probabilities correctly", {
  pval <- 99000 / 100000

  expect_equal(
    nqx_select(sel_tbl, x_sel = 20, t = 0, n = 1),
    1 - pval
  )
})

# ------------------------------------------------------------
# nmxq_select
# ------------------------------------------------------------

test_that("nmxq_select computes deferred select death probabilities correctly", {
  # {}_{1|1} q_[20]
  # = {}_1 p_[20] * {}_1 q_[20]+1
  p1 <- 99000 / 100000
  q1_after <- 1 - (97000 / 99000)

  expect_equal(
    nmxq_select(sel_tbl, x_sel = 20, t = 0, n = 1, m = 1),
    p1 * q1_after
  )
})

# ------------------------------------------------------------
# vectorization
# ------------------------------------------------------------

test_that("select-table functions vectorize correctly", {
  vals1 <- lx_select(sel_tbl, x_sel = c(20, 21, 22), t = 0)
  expect_equal(length(vals1), 3)

  vals2 <- npx_select(sel_tbl, x_sel = c(20, 21), t = 0, n = 1)
  expect_equal(length(vals2), 2)

  vals3 <- nqx_select(sel_tbl, x_sel = c(20, 21), t = 0, n = c(1, 2))
  expect_equal(length(vals3), 2)

  vals4 <- nmxq_select(sel_tbl, x_sel = c(20, 21), t = 0, n = 1, m = 1)
  expect_equal(length(vals4), 2)
})

# ------------------------------------------------------------
# invalid constructor inputs
# ------------------------------------------------------------

test_that("select_life_table rejects inconsistent attained ages", {
  expect_error(
    select_life_table(
      x_sel = c(20, 20),
      duration = c(0, 1),
      attained_age = c(20, 25),  # wrong
      lx = c(100000, 99000)
    )
  )
})

test_that("select_life_table rejects unequal input lengths", {
  expect_error(
    select_life_table(
      x_sel = c(20, 20),
      duration = c(0, 1),
      attained_age = c(20, 21),
      lx = c(100000)
    )
  )
})

# ------------------------------------------------------------
# invalid probability inputs
# ------------------------------------------------------------

test_that("select probability functions reject invalid n and m", {
  expect_error(
    npx_select(sel_tbl, x_sel = 20, t = 0, n = -1)
  )

  expect_error(
    npx_select(sel_tbl, x_sel = 20, t = 0, n = 1.5)
  )

  expect_error(
    nmxq_select(sel_tbl, x_sel = 20, t = 0, n = 1, m = -1)
  )

  expect_error(
    nmxq_select(sel_tbl, x_sel = 20, t = 0, n = 1, m = 1.5)
  )
})
