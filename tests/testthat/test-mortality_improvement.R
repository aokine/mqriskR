test_that("qx_proj returns base rates when projection year equals base year", {
  qx <- c(0.01, 0.02, 0.03)
  AA <- c(0.01, 0.02, 0.03)

  expect_equal(
    qx_proj(qx, AA, base_year = 2000, proj_year = 2000),
    qx
  )
})

test_that("px_proj equals 1 - qx_proj", {
  qx <- c(0.01, 0.02, 0.03)
  AA <- c(0.01, 0.02, 0.03)

  qproj <- qx_proj(qx, AA, 2000, 2005)
  pproj <- px_proj(qx, AA, 2000, 2005)

  expect_equal(pproj, 1 - qproj)
})

test_that("mortality improvement lowers projected q when AA is positive", {
  qx <- c(0.01, 0.02, 0.03)
  AA <- c(0.01, 0.01, 0.01)

  qproj <- qx_proj(qx, AA, 2000, 2010)

  expect_true(all(qproj < qx))
})

test_that("tpx_improved equals 1 when n = 0", {
  expect_equal(
    tpx_improved(
      x0 = 60, n = 0,
      qx_base_vec = c(0.01, 0.02),
      AAx_vec = c(0.01, 0.01),
      base_year = 2000,
      issue_year = 2005
    ),
    1
  )
})

test_that("tpx_improved matches manual product", {
  qx_base <- c(0.10, 0.20, 0.30)
  AA <- c(0.00, 0.00, 0.00)

  expected <- (1 - 0.10) * (1 - 0.20) * (1 - 0.30)

  actual <- tpx_improved(
    x0 = 50, n = 3,
    qx_base_vec = qx_base,
    AAx_vec = AA,
    base_year = 2000,
    issue_year = 2000
  )

  expect_equal(actual, expected)
})

test_that("tpx_improved increases when mortality improvement is introduced", {
  qx_base <- c(0.02, 0.03, 0.04, 0.05)
  AA0 <- c(0, 0, 0, 0)
  AA1 <- c(0.01, 0.01, 0.01, 0.01)

  no_improve <- tpx_improved(
    x0 = 60, n = 4,
    qx_base_vec = qx_base,
    AAx_vec = AA0,
    base_year = 2000,
    issue_year = 2000
  )

  improve <- tpx_improved(
    x0 = 60, n = 4,
    qx_base_vec = qx_base,
    AAx_vec = AA1,
    base_year = 2000,
    issue_year = 2000
  )

  expect_true(improve >= no_improve)
})

test_that("axn_improved matches textbook-style direct summation", {
  x0 <- 60
  n <- 2
  i <- 0.06
  qx_base <- c(0.0479, 0.0539)
  AA <- c(0.01, 0.01)

  p1 <- 1 - qx_proj(qx_base[1], AA[1], 2000, 2020)
  p2 <- 1 - qx_proj(qx_base[2], AA[2], 2000, 2021)

  expected <- discount(i, 1) * p1 +
    discount(i, 2) * p1 * p2

  actual <- axn_improved(
    x0 = x0, n = n, i = i,
    qx_base_vec = qx_base,
    AAx_vec = AA,
    base_year = 2000,
    issue_year = 2020
  )

  expect_equal(actual, expected)
})

test_that("deferred temporary annuity under improvement reduces to shifted sum", {
  x0 <- 60
  u <- 2
  n <- 2
  i <- 0.05
  qx_base <- c(0.01, 0.02, 0.03, 0.04)
  AA <- c(0.00, 0.00, 0.00, 0.00)

  expected <- sum(sapply((u + 1):(u + n), function(t) {
    discount(i, t) * tpx_improved(
      x0 = x0, n = t,
      qx_base_vec = qx_base,
      AAx_vec = AA,
      base_year = 2000,
      issue_year = 2000
    )
  }))

  actual <- naxn_improved(
    x0 = x0, u = u, n = n, i = i,
    qx_base_vec = qx_base,
    AAx_vec = AA,
    base_year = 2000,
    issue_year = 2000
  )

  expect_equal(actual, expected)
})

test_that("whole life annuity under improvement is at least as large as temporary annuity with same horizon", {
  x0 <- 60
  i <- 0.05
  qx_base <- c(0.01, 0.02, 0.03, 0.04, 0.05)
  AA <- c(0.01, 0.01, 0.01, 0.01, 0.01)

  whole <- ax_improved(
    x0 = x0, i = i,
    qx_base_vec = qx_base,
    AAx_vec = AA,
    base_year = 2000,
    issue_year = 2010
  )

  temp <- axn_improved(
    x0 = x0, n = 5, i = i,
    qx_base_vec = qx_base,
    AAx_vec = AA,
    base_year = 2000,
    issue_year = 2010
  )

  expect_true(whole >= temp)
})

test_that("functions reject invalid inputs", {
  expect_error(qx_proj(0.01, 0.01, 2000, 2000.5), "proj_year")
  expect_error(qx_proj(-0.01, 0.01, 2000, 2005), "qx_base")
  expect_error(px_proj(0.01, 1.2, 2000, 2005), "AAx")

  expect_error(
    tpx_improved(
      x0 = 60, n = 3,
      qx_base_vec = c(0.01, 0.02),
      AAx_vec = c(0.01, 0.01),
      base_year = 2000,
      issue_year = 2005
    ),
    "length at least n"
  )
})

test_that("qx_proj vectorizes over qx and AAx", {
  out <- qx_proj(
    qx_base = c(0.01, 0.02, 0.03),
    AAx = 0.01,
    base_year = 2000,
    proj_year = 2010
  )

  expect_length(out, 3)
  expect_true(all(is.finite(out)))
})

test_that("projected annuity values are finite in standard cases", {
  val1 <- axn_improved(
    x0 = 60, n = 3, i = 0.05,
    qx_base_vec = c(0.02, 0.03, 0.04),
    AAx_vec = c(0.01, 0.01, 0.01),
    base_year = 2000,
    issue_year = 2020
  )

  val2 <- naxn_improved(
    x0 = 60, u = 2, n = 3, i = 0.05,
    qx_base_vec = c(0.02, 0.03, 0.04, 0.05, 0.06),
    AAx_vec = c(0.01, 0.01, 0.01, 0.01, 0.01),
    base_year = 2000,
    issue_year = 2020
  )

  expect_true(is.finite(val1))
  expect_true(is.finite(val2))
})
