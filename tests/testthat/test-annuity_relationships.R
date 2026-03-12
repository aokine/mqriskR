test_that("annuity_identity_ax matches ax", {
  i <- 0.05
  x <- 40

  expect_equal(
    annuity_identity_ax(x, i, model = "uniform", omega = 100),
    ax(x, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_adotx matches adotx", {
  i <- 0.05
  x <- 40

  expect_equal(
    annuity_identity_adotx(x, i, model = "uniform", omega = 100),
    adotx(x, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_abarx matches abarx", {
  i <- 0.05
  x <- 40
  lambda <- 0.02

  expect_equal(
    annuity_identity_abarx(x, i, model = "exponential", lambda = lambda),
    abarx(x, i, model = "exponential", lambda = lambda),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_axn matches axn", {
  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_axn(x, n, i, model = "uniform", omega = 100),
    axn(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_adotxn matches adotxn", {
  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_adotxn(x, n, i, model = "uniform", omega = 100),
    adotxn(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_abarxn matches abarxn", {
  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_equal(
    annuity_identity_abarxn(x, n, i, model = "exponential", lambda = lambda),
    abarxn(x, n, i, model = "exponential", lambda = lambda),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_nax matches nax", {
  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_nax(x, n, i, model = "uniform", omega = 100),
    nax(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_nadotx matches nadotx", {
  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_nadotx(x, n, i, model = "uniform", omega = 100),
    nadotx(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("annuity_identity_nabarx matches nabarx", {
  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_equal(
    annuity_identity_nabarx(x, n, i, model = "exponential", lambda = lambda),
    nabarx(x, n, i, model = "exponential", lambda = lambda),
    tolerance = 1e-10
  )
})

test_that("whole life immediate identity agrees with textbook formula", {
  i <- 0.05
  d <- i / (1 + i)
  v <- 1 / (1 + i)
  x <- 40

  expect_equal(
    annuity_identity_ax(x, i, model = "uniform", omega = 100),
    (v - Ax(x, i, model = "uniform", omega = 100)) / d,
    tolerance = 1e-10
  )
})

test_that("whole life due identity agrees with textbook formula", {
  i <- 0.05
  d <- i / (1 + i)
  x <- 40

  expect_equal(
    annuity_identity_adotx(x, i, model = "uniform", omega = 100),
    (1 - Ax(x, i, model = "uniform", omega = 100)) / d,
    tolerance = 1e-10
  )
})

test_that("continuous whole life identity agrees with textbook formula", {
  i <- 0.05
  delta <- log(1 + i)
  x <- 40
  lambda <- 0.02

  expect_equal(
    annuity_identity_abarx(x, i, model = "exponential", lambda = lambda),
    (1 - Abarx(x, i, model = "exponential", lambda = lambda)) / delta,
    tolerance = 1e-10
  )
})

test_that("temporary immediate identity agrees with textbook formula", {
  i <- 0.05
  d <- i / (1 + i)
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_axn(x, n, i, model = "uniform", omega = 100),
    (1 - Axn(x, n, i, model = "uniform", omega = 100)) / d - 1 +
      nEx(x, n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("temporary due identity agrees with textbook formula", {
  i <- 0.05
  d <- i / (1 + i)
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_adotxn(x, n, i, model = "uniform", omega = 100),
    (1 - Axn(x, n, i, model = "uniform", omega = 100)) / d,
    tolerance = 1e-10
  )
})

test_that("continuous temporary identity agrees with textbook formula", {
  i <- 0.05
  delta <- log(1 + i)
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_equal(
    annuity_identity_abarxn(x, n, i, model = "exponential", lambda = lambda),
    (1 - Abarxn(x, n, i, model = "exponential", lambda = lambda)) / delta,
    tolerance = 1e-10
  )
})

test_that("deferred immediate identity agrees with pure endowment relation", {
  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_nax(x, n, i, model = "uniform", omega = 100),
    nEx(x, n, i, model = "uniform", omega = 100) *
      ax(x + n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("deferred due identity agrees with pure endowment relation", {
  i <- 0.05
  x <- 40
  n <- 10

  expect_equal(
    annuity_identity_nadotx(x, n, i, model = "uniform", omega = 100),
    nEx(x, n, i, model = "uniform", omega = 100) *
      adotx(x + n, i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("deferred continuous identity agrees with pure endowment relation", {
  i <- 0.05
  x <- 40
  n <- 10
  lambda <- 0.02

  expect_equal(
    annuity_identity_nabarx(x, n, i, model = "exponential", lambda = lambda),
    nEx(x, n, i, model = "exponential", lambda = lambda) *
      abarx(x + n, i, model = "exponential", lambda = lambda),
    tolerance = 1e-10
  )
})

test_that("relationship functions vectorize over x", {
  i <- 0.05
  x <- c(40, 45, 50)

  expect_length(
    annuity_identity_ax(x, i, model = "uniform", omega = 100),
    3
  )

  expect_length(
    annuity_identity_adotx(x, i, model = "uniform", omega = 100),
    3
  )

  expect_length(
    annuity_identity_abarx(x, i, model = "exponential", lambda = 0.02),
    3
  )
})

test_that("relationship functions vectorize over x and n", {
  i <- 0.05
  x <- c(40, 45)
  n <- c(10, 5)

  expect_length(
    annuity_identity_axn(x, n, i, model = "uniform", omega = 100),
    2
  )

  expect_length(
    annuity_identity_adotxn(x, n, i, model = "uniform", omega = 100),
    2
  )

  expect_length(
    annuity_identity_abarxn(x, n, i, model = "exponential", lambda = 0.02),
    2
  )

  expect_length(
    annuity_identity_nax(x, n, i, model = "uniform", omega = 100),
    2
  )

  expect_length(
    annuity_identity_nadotx(x, n, i, model = "uniform", omega = 100),
    2
  )

  expect_length(
    annuity_identity_nabarx(x, n, i, model = "exponential", lambda = 0.02),
    2
  )
})

test_that("relationship functions reject incompatible x and n lengths", {
  expect_error(
    annuity_identity_axn(c(40, 45, 50), c(10, 5), 0.05, model = "uniform", omega = 100)
  )

  expect_error(
    annuity_identity_nax(c(40, 45, 50), c(10, 5), 0.05, model = "uniform", omega = 100)
  )
})

test_that("relationship functions return finite values in standard cases", {
  i <- 0.05

  expect_true(is.finite(
    annuity_identity_ax(40, i, model = "uniform", omega = 100)
  ))

  expect_true(is.finite(
    annuity_identity_adotx(40, i, model = "uniform", omega = 100)
  ))

  expect_true(is.finite(
    annuity_identity_abarx(40, i, model = "exponential", lambda = 0.02)
  ))

  expect_true(is.finite(
    annuity_identity_axn(40, 10, i, model = "uniform", omega = 100)
  ))

  expect_true(is.finite(
    annuity_identity_adotxn(40, 10, i, model = "uniform", omega = 100)
  ))

  expect_true(is.finite(
    annuity_identity_abarxn(40, 10, i, model = "exponential", lambda = 0.02)
  ))
})
