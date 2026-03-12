test_that("whole life m-thly due and immediate differ by 1/m", {
  i <- 0.05
  m <- 12

  expect_equal(
    adotx_m(40, m, i, model = "uniform", omega = 100),
    ax_m(40, m, i, model = "uniform", omega = 100) + 1 / m,
    tolerance = 1e-10
  )

  expect_equal(
    adotx_m(40, m, i, model = "exponential", lambda = 0.02),
    ax_m(40, m, i, model = "exponential", lambda = 0.02) + 1 / m,
    tolerance = 1e-10
  )
})

test_that("temporary m-thly due and immediate relation matches textbook formula", {
  i <- 0.05
  m <- 4
  x <- 40
  n <- 10

  expect_equal(
    adotxn_m(x, n, m, i, model = "uniform", omega = 100),
    axn_m(x, n, m, i, model = "uniform", omega = 100) + (1 / m) * (1 - nEx(x, n, i, model = "uniform", omega = 100)),
    tolerance = 1e-10
  )
})

test_that("deferred immediate identity holds", {
  i <- 0.05
  m <- 12
  x <- 40
  n <- 10

  expect_equal(
    nax_m(x, n, m, i, model = "uniform", omega = 100),
    nEx(x, n, i, model = "uniform", omega = 100) *
      ax_m(x + n, m, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    nax_m(x, n, m, i, model = "exponential", lambda = 0.02),
    nEx(x, n, i, model = "exponential", lambda = 0.02) *
      ax_m(x + n, m, i, model = "exponential", lambda = 0.02),
    tolerance = 1e-10
  )
})

test_that("deferred due identity holds", {
  i <- 0.05
  m <- 12
  x <- 40
  n <- 10

  expect_equal(
    nadotx_m(x, n, m, i, model = "uniform", omega = 100),
    nEx(x, n, i, model = "uniform", omega = 100) *
      adotx_m(x + n, m, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    nadotx_m(x, n, m, i, model = "exponential", lambda = 0.02),
    nEx(x, n, i, model = "exponential", lambda = 0.02) *
      adotx_m(x + n, m, i, model = "exponential", lambda = 0.02),
    tolerance = 1e-10
  )
})

test_that("whole life decomposition into temporary plus deferred holds", {
  i <- 0.05
  m <- 4
  x <- 40
  n <- 10

  expect_equal(
    ax_m(x, m, i, model = "uniform", omega = 100),
    axn_m(x, n, m, i, model = "uniform", omega = 100) +
      nax_m(x, n, m, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    adotx_m(x, m, i, model = "uniform", omega = 100),
    adotxn_m(x, n, m, i, model = "uniform", omega = 100) +
      nadotx_m(x, n, m, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("actuarial accumulated values equal APV divided by nEx", {
  i <- 0.05
  m <- 12
  x <- 40
  n <- 10

  ex <- nEx(x, n, i, model = "uniform", omega = 100)

  expect_equal(
    sxn_m(x, n, m, i, model = "uniform", omega = 100),
    axn_m(x, n, m, i, model = "uniform", omega = 100) / ex,
    tolerance = 1e-10
  )

  expect_equal(
    sdotxn_m(x, n, m, i, model = "uniform", omega = 100),
    adotxn_m(x, n, m, i, model = "uniform", omega = 100) / ex,
    tolerance = 1e-10
  )
})

test_that("m-thly annuities vectorize over x", {
  i <- 0.05
  m <- 12
  x <- c(40, 45, 50)

  out1 <- ax_m(x, m, i, model = "uniform", omega = 100)
  out2 <- adotx_m(x, m, i, model = "uniform", omega = 100)

  expect_length(out1, 3)
  expect_length(out2, 3)
  expect_true(all(is.finite(out1)))
  expect_true(all(is.finite(out2)))
})

test_that("m-thly temporary and deferred annuities vectorize over x and n", {
  i <- 0.05
  m <- 4
  x <- c(40, 45)
  n <- c(10, 5)

  out1 <- axn_m(x, n, m, i, model = "uniform", omega = 100)
  out2 <- adotxn_m(x, n, m, i, model = "uniform", omega = 100)
  out3 <- nax_m(x, n, m, i, model = "uniform", omega = 100)
  out4 <- nadotx_m(x, n, m, i, model = "uniform", omega = 100)

  expect_length(out1, 2)
  expect_length(out2, 2)
  expect_length(out3, 2)
  expect_length(out4, 2)

  expect_true(all(is.finite(out1)))
  expect_true(all(is.finite(out2)))
  expect_true(all(is.finite(out3)))
  expect_true(all(is.finite(out4)))
})

test_that("m-thly annuity values lie between annual and continuous values in standard cases", {
  i <- 0.05
  x <- 40

  a_annual <- ax(x, i, model = "uniform", omega = 100)
  a_cont <- abarx(x, i, model = "uniform", omega = 100)
  a_m <- ax_m(x, 12, i, model = "uniform", omega = 100)

  expect_true(a_annual < a_m)
  expect_true(a_m < a_cont)

  ad_annual <- adotx(x, i, model = "uniform", omega = 100)
  ad_m <- adotx_m(x, 12, i, model = "uniform", omega = 100)

  expect_true(a_m < ad_m)
  expect_true(ad_m < ad_annual)
})

test_that("m-thly temporary annuity is not larger than whole life annuity", {
  i <- 0.05
  m <- 12
  x <- 40
  n <- 10

  expect_true(
    axn_m(x, n, m, i, model = "uniform", omega = 100) <=
      ax_m(x, m, i, model = "uniform", omega = 100)
  )

  expect_true(
    adotxn_m(x, n, m, i, model = "uniform", omega = 100) <=
      adotx_m(x, m, i, model = "uniform", omega = 100)
  )
})

test_that("deferred annuity is zero when nEx is zero under finite support", {
  i <- 0.05
  m <- 12

  expect_equal(
    nax_m(99, 2, m, i, model = "uniform", omega = 100),
    0
  )

  expect_equal(
    nadotx_m(99, 2, m, i, model = "uniform", omega = 100),
    0
  )
})

test_that("temporary annuity with n = 0 is zero", {
  i <- 0.05
  m <- 12

  expect_equal(axn_m(40, 0, m, i, model = "uniform", omega = 100), 0)
  expect_equal(adotxn_m(40, 0, m, i, model = "uniform", omega = 100), 0)
})

test_that("functions reject invalid m", {
  expect_error(
    ax_m(40, m = 0, i = 0.05, model = "uniform", omega = 100),
    "m must be a positive integer"
  )

  expect_error(
    ax_m(40, m = 2.5, i = 0.05, model = "uniform", omega = 100),
    "m must be a positive integer"
  )
})

test_that("functions reject incompatible x and n lengths", {
  expect_error(
    axn_m(c(40, 45), c(10, 5, 3), m = 12, i = 0.05, model = "uniform", omega = 100),
    "compatible lengths"
  )
})

test_that("functions reject non-integer n * m", {
  expect_error(
    axn_m(40, n = 10.1, m = 4, i = 0.05, model = "uniform", omega = 100),
    "n \\* m must be an integer"
  )
})
