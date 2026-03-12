test_that("ax equals sum of annual pure endowments", {
  library(mqriskR)

  i <- 0.05
  x <- 40

  expected <- sum(sapply(1:(100 - x), function(t) {
    nEx(x, t, i, model = "uniform", omega = 100)
  }))

  expect_equal(
    ax(x, i, model = "uniform", omega = 100),
    expected,
    tolerance = 1e-10
  )
})

test_that("adotx equals 1 + ax", {
  library(mqriskR)

  i <- 0.05
  x <- 40

  expect_equal(
    adotx(x, i, model = "uniform", omega = 100),
    1 + ax(x, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("abarx under exponential model matches closed form", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  lambda <- 0.02
  delta <- log(1 + i)

  expected <- 1 / (delta + lambda)

  expect_equal(
    abarx(x, i, model = "exponential", lambda = lambda),
    expected,
    tolerance = 1e-10
  )
})

test_that("axn equals sum of first n annual pure endowments", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expected <- sum(sapply(1:n, function(t) {
    nEx(x, t, i, model = "uniform", omega = 100)
  }))

  expect_equal(
    axn(x, n, i, model = "uniform", omega = 100),
    expected,
    tolerance = 1e-10
  )
})

test_that("adotxn equals sum from t = 0 to n - 1 of pure endowments", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expected <- sum(sapply(0:(n - 1), function(t) {
    nEx(x, t, i, model = "uniform", omega = 100)
  }))

  expect_equal(
    adotxn(x, n, i, model = "uniform", omega = 100),
    expected,
    tolerance = 1e-10
  )
})

test_that("adotxn equals axn + 1 - nEx", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    adotxn(x, n, i, model = "uniform", omega = 100),
    axn(x, n, i, model = "uniform", omega = 100) +
      1 -
      nEx(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("abarxn equals integral form under exponential model", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02
  delta <- log(1 + i)

  expected <- (1 - exp(-(delta + lambda) * n)) / (delta + lambda)

  expect_equal(
    abarxn(x, n, i, model = "exponential", lambda = lambda),
    expected,
    tolerance = 1e-10
  )
})

test_that("nax equals nEx times ax at age x+n", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    nax(x, n, i, model = "uniform", omega = 100),
    nEx(x, n, i, model = "uniform", omega = 100) *
      ax(x + n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("nadotx equals nEx times adotx at age x+n", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    nadotx(x, n, i, model = "uniform", omega = 100),
    nEx(x, n, i, model = "uniform", omega = 100) *
      adotx(x + n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("nabarx equals nEx times abarx at age x+n", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_equal(
    nabarx(x, n, i, model = "exponential", lambda = lambda),
    nEx(x, n, i, model = "exponential", lambda = lambda) *
      abarx(x + n, i, model = "exponential", lambda = lambda),
    tolerance = 1e-10
  )
})

test_that("whole life annuity decomposes into temporary plus deferred", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    ax(x, i, model = "uniform", omega = 100),
    axn(x, n, i, model = "uniform", omega = 100) +
      nax(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("whole life annuity-due decomposes into temporary plus deferred", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    adotx(x, i, model = "uniform", omega = 100),
    adotxn(x, n, i, model = "uniform", omega = 100) +
      nadotx(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("continuous whole life annuity decomposes into temporary plus deferred", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_equal(
    abarx(x, i, model = "exponential", lambda = lambda),
    abarxn(x, n, i, model = "exponential", lambda = lambda) +
      nabarx(x, n, i, model = "exponential", lambda = lambda),
    tolerance = 1e-10
  )
})

test_that("sxn equals axn divided by nEx", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    sxn(x, n, i, model = "uniform", omega = 100),
    axn(x, n, i, model = "uniform", omega = 100) /
      nEx(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("sdotxn equals adotxn divided by nEx", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    sdotxn(x, n, i, model = "uniform", omega = 100),
    adotxn(x, n, i, model = "uniform", omega = 100) /
      nEx(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("sbarxn equals abarxn divided by nEx", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_equal(
    sbarxn(x, n, i, model = "exponential", lambda = lambda),
    abarxn(x, n, i, model = "exponential", lambda = lambda) /
      nEx(x, n, i, model = "exponential", lambda = lambda),
    tolerance = 1e-10
  )
})

test_that("annuity functions vectorize correctly over x", {
  library(mqriskR)

  i <- 0.05
  x <- c(40, 45, 50)
  n <- 10

  expect_length(ax(x, i, model = "uniform", omega = 100), 3)
  expect_length(adotx(x, i, model = "uniform", omega = 100), 3)
  expect_length(abarx(x, i, model = "exponential", lambda = 0.02), 3)
  expect_length(axn(x, n, i, model = "uniform", omega = 100), 3)
  expect_length(adotxn(x, n, i, model = "uniform", omega = 100), 3)
  expect_length(abarxn(x, n, i, model = "exponential", lambda = 0.02), 3)
  expect_length(nax(x, n, i, model = "uniform", omega = 100), 3)
  expect_length(nadotx(x, n, i, model = "uniform", omega = 100), 3)
  expect_length(nabarx(x, n, i, model = "exponential", lambda = 0.02), 3)
})

test_that("annuity functions vectorize correctly over x and n", {
  library(mqriskR)

  i <- 0.05
  x <- c(40, 45)
  n <- c(10, 5)

  expect_length(axn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(adotxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(abarxn(x, n, i, model = "exponential", lambda = 0.02), 2)
  expect_length(nax(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(nadotx(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(nabarx(x, n, i, model = "exponential", lambda = 0.02), 2)
  expect_length(sxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(sdotxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(sbarxn(x, n, i, model = "exponential", lambda = 0.02), 2)
})

test_that("temporary annuities are not larger than corresponding whole life annuities", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_true(axn(x, n, i, model = "uniform", omega = 100) <=
                ax(x, i, model = "uniform", omega = 100))

  expect_true(adotxn(x, n, i, model = "uniform", omega = 100) <=
                adotx(x, i, model = "uniform", omega = 100))

  expect_true(abarxn(x, n, i, model = "exponential", lambda = lambda) <=
                abarx(x, i, model = "exponential", lambda = lambda))
})

test_that("due annuities are at least as large as immediate annuities", {
  library(mqriskR)

  i <- 0.05
  x <- 40
  n <- 10

  expect_true(adotx(x, i, model = "uniform", omega = 100) >=
                ax(x, i, model = "uniform", omega = 100))

  expect_true(adotxn(x, n, i, model = "uniform", omega = 100) >=
                axn(x, n, i, model = "uniform", omega = 100))
})

test_that("functions reject incompatible x and n lengths", {
  library(mqriskR)

  expect_error(
    axn(c(40, 45, 50), c(10, 5), 0.05, model = "uniform", omega = 100)
  )

  expect_error(
    nax(c(40, 45, 50), c(10, 5), 0.05, model = "uniform", omega = 100)
  )
})

test_that("annual annuity functions reject non-integer n where appropriate", {
  library(mqriskR)

  expect_error(
    axn(40, 10.5, 0.05, model = "uniform", omega = 100),
    "integer"
  )

  expect_error(
    adotxn(40, 10.5, 0.05, model = "uniform", omega = 100),
    "integer"
  )

  expect_error(
    nax(40, 10.5, 0.05, model = "uniform", omega = 100),
    "integer"
  )

  expect_error(
    nadotx(40, 10.5, 0.05, model = "uniform", omega = 100),
    "integer"
  )
})
