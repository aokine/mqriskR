test_that("double_force_i computes i' = (1+i)^2 - 1", {
  expect_equal(double_force_i(0.05), (1.05)^2 - 1)
  expect_equal(double_force_i(c(0.02, 0.05)),
               c((1.02)^2 - 1, (1.05)^2 - 1))
})

test_that("double_force_i rejects invalid interest rates", {
  expect_error(double_force_i(-1))
  expect_error(double_force_i(c(0.05, NA)))
})

test_that("double_force_delta computes 2 * delta", {
  expect_equal(double_force_delta(0.08), 0.16)
  expect_equal(double_force_delta(c(0.01, 0.03)), c(0.02, 0.06))
})

test_that("double_force_delta rejects non-finite inputs", {
  expect_error(double_force_delta(Inf))
  expect_error(double_force_delta(NA))
})

test_that("udd_continuous_multiplier computes i / delta", {
  i <- 0.05
  delta <- interest_convert(i = i)$delta
  expect_equal(udd_continuous_multiplier(i), i / delta)
})

test_that("udd_continuous_multiplier handles zero interest correctly", {
  expect_equal(udd_continuous_multiplier(0), 1)
})

test_that("udd_continuous_multiplier rejects invalid interest rates", {
  expect_error(udd_continuous_multiplier(-1))
  expect_error(udd_continuous_multiplier(NA))
})

test_that("udd_mthly_multiplier computes i / i^(m)", {
  i <- 0.06
  m <- 12
  im <- interest_convert(i = i, m = m)$im
  expect_equal(udd_mthly_multiplier(i, m), i / im)
})

test_that("udd_mthly_multiplier handles zero interest correctly", {
  expect_equal(udd_mthly_multiplier(0, 12), 1)
})

test_that("udd_mthly_multiplier rejects invalid inputs", {
  expect_error(udd_mthly_multiplier(0.05, 0))
  expect_error(udd_mthly_multiplier(0.05, 2.5))
  expect_error(udd_mthly_multiplier(-1, 12))
})

test_that("Abarx_udd computes (i/delta) * Ax", {
  Ax <- 0.30
  i <- 0.05
  expected <- udd_continuous_multiplier(i) * Ax
  expect_equal(Abarx_udd(Ax, i), expected)
})

test_that("Abarx_udd vectorizes", {
  Ax <- c(0.20, 0.40)
  i <- c(0.03, 0.05)
  expected <- udd_continuous_multiplier(i) * Ax
  expect_equal(Abarx_udd(Ax, i), expected)
})

test_that("Abarx_udd rejects invalid Ax", {
  expect_error(Abarx_udd(NA, 0.05))
})

test_that("Abarxn1_udd computes (i/delta) * Axn1", {
  Axn1 <- 0.18
  i <- 0.04
  expected <- udd_continuous_multiplier(i) * Axn1
  expect_equal(Abarxn1_udd(Axn1, i), expected)
})

test_that("nAbarx_udd computes (i/delta) * nAx", {
  nAx <- 0.11
  i <- 0.05
  expected <- udd_continuous_multiplier(i) * nAx
  expect_equal(nAbarx_udd(nAx, i), expected)
})

test_that("Abarxn_udd computes continuous endowment approximation", {
  Axn1 <- 0.25
  nEx <- 0.40
  i <- 0.05
  expected <- udd_continuous_multiplier(i) * Axn1 + nEx
  expect_equal(Abarxn_udd(Axn1, nEx, i), expected)
})

test_that("Abarxn_udd rejects invalid inputs", {
  expect_error(Abarxn_udd(NA, 0.4, 0.05))
  expect_error(Abarxn_udd(0.2, NA, 0.05))
})

test_that("Ax_m_udd computes (i / i^(m)) * Ax", {
  Ax <- 0.33
  i <- 0.06
  m <- 4
  expected <- udd_mthly_multiplier(i, m) * Ax
  expect_equal(Ax_m_udd(Ax, i, m), expected)
})

test_that("Axn1_m_udd computes (i / i^(m)) * Axn1", {
  Axn1 <- 0.17
  i <- 0.06
  m <- 12
  expected <- udd_mthly_multiplier(i, m) * Axn1
  expect_equal(Axn1_m_udd(Axn1, i, m), expected)
})

test_that("nAx_m_udd computes (i / i^(m)) * nAx", {
  nAx <- 0.09
  i <- 0.04
  m <- 2
  expected <- udd_mthly_multiplier(i, m) * nAx
  expect_equal(nAx_m_udd(nAx, i, m), expected)
})

test_that("Axn_m_udd computes m-thly endowment approximation", {
  Axn1 <- 0.20
  nEx <- 0.50
  i <- 0.05
  m <- 12
  expected <- udd_mthly_multiplier(i, m) * Axn1 + nEx
  expect_equal(Axn_m_udd(Axn1, nEx, i, m), expected)
})

test_that("Axn_m_udd rejects invalid inputs", {
  expect_error(Axn_m_udd(NA, 0.5, 0.05, 12))
  expect_error(Axn_m_udd(0.2, NA, 0.05, 12))
})

test_that("continuous and m-thly multipliers are close for large m", {
  i <- 0.05
  expect_equal(
    udd_mthly_multiplier(i, 10000),
    udd_continuous_multiplier(i),
    tolerance = 1e-4
  )
})
