test_that("Abarxn decomposes into continuous term insurance plus pure endowment", {
  i <- 0.05

  expect_equal(
    Abarxn(40, 10, i, model = "uniform", omega = 100),
    Abarxn1(40, 10, i, model = "uniform", omega = 100) +
      discount(i, 10) * tpx(10, x = 40, model = "uniform", omega = 100)
  )
})

test_that("nAbarx equals v^n * npx * Abarx at age x+n", {
  i <- 0.05

  expect_equal(
    nAbarx(40, 10, i, model = "uniform", omega = 100),
    discount(i, 10) *
      tpx(10, x = 40, model = "uniform", omega = 100) *
      Abarx(50, i, model = "uniform", omega = 100)
  )
})

test_that("A2barx evaluates Abarx at doubled force", {
  i <- 0.05

  expect_equal(
    A2barx(40, i, model = "uniform", omega = 100),
    Abarx(40, double_force_i(i), model = "uniform", omega = 100)
  )
})

test_that("A2barxn1 evaluates Abarxn1 at doubled force", {
  i <- 0.05

  expect_equal(
    A2barxn1(40, 10, i, model = "uniform", omega = 100),
    Abarxn1(40, 10, double_force_i(i), model = "uniform", omega = 100)
  )
})

test_that("A2nAbarx evaluates nAbarx at doubled force", {
  i <- 0.05

  expect_equal(
    A2nAbarx(40, 10, i, model = "uniform", omega = 100),
    nAbarx(40, 10, double_force_i(i), model = "uniform", omega = 100)
  )
})

test_that("A2barxn decomposes into second moments of term and pure endowment", {
  i <- 0.05
  i2 <- double_force_i(i)

  expect_equal(
    A2barxn(40, 10, i, model = "uniform", omega = 100),
    A2barxn1(40, 10, i, model = "uniform", omega = 100) +
      discount(i2, 10) * tpx(10, x = 40, model = "uniform", omega = 100)
  )
})

test_that("var_Abarx equals A2barx - Abarx^2", {
  i <- 0.05

  expect_equal(
    var_Abarx(40, i, model = "uniform", omega = 100),
    A2barx(40, i, model = "uniform", omega = 100) -
      Abarx(40, i, model = "uniform", omega = 100)^2
  )
})

test_that("var_Abarxn1 equals A2barxn1 - Abarxn1^2", {
  i <- 0.05

  expect_equal(
    var_Abarxn1(40, 10, i, model = "uniform", omega = 100),
    A2barxn1(40, 10, i, model = "uniform", omega = 100) -
      Abarxn1(40, 10, i, model = "uniform", omega = 100)^2
  )
})

test_that("var_nAbarx equals A2nAbarx - nAbarx^2", {
  i <- 0.05

  expect_equal(
    var_nAbarx(40, 10, i, model = "uniform", omega = 100),
    A2nAbarx(40, 10, i, model = "uniform", omega = 100) -
      nAbarx(40, 10, i, model = "uniform", omega = 100)^2
  )
})

test_that("var_Abarxn equals A2barxn - Abarxn^2", {
  i <- 0.05

  expect_equal(
    var_Abarxn(40, 10, i, model = "uniform", omega = 100),
    A2barxn(40, 10, i, model = "uniform", omega = 100) -
      Abarxn(40, 10, i, model = "uniform", omega = 100)^2
  )
})

test_that("continuous insurance functions vectorize correctly", {
  i <- 0.05

  expect_length(Abarxn1(c(40, 45), c(10, 5), i, model = "uniform", omega = 100), 2)
  expect_length(nAbarx(c(40, 45), c(10, 5), i, model = "uniform", omega = 100), 2)
  expect_length(Abarxn(c(40, 45), c(10, 5), i, model = "uniform", omega = 100), 2)
})

test_that("continuous insurance functions reject invalid inputs", {
  i <- 0.05

  expect_error(Abarx(-1, i, model = "uniform", omega = 100))
  expect_error(Abarxn1(40, -1, i, model = "uniform", omega = 100))
  expect_error(nAbarx(40, -1, i, model = "uniform", omega = 100))
  expect_error(Abarx(40, -1, model = "uniform", omega = 100))
})

test_that("Abarx matches closed form under exponential survival", {
  i <- 0.05
  lambda <- 0.25
  delta <- interest_convert(i = i)$delta

  expected <- lambda / (delta + lambda)

  expect_equal(
    Abarx(40, i, model = "exponential", lambda = lambda),
    expected,
    tolerance = 1e-8
  )
})

test_that("Abarxn matches closed form under exponential survival", {
  i <- 0.05
  lambda <- 0.25
  delta <- interest_convert(i = i)$delta
  n <- 10

  expected <- (lambda + delta * exp(-(delta + lambda) * n)) / (delta + lambda)

  expect_equal(
    Abarxn(40, n, i, model = "exponential", lambda = lambda),
    expected,
    tolerance = 1e-8
  )
})

test_that("Abarxn1 equals Abarx minus nAbarx", {
  i <- 0.05

  expect_equal(
    Abarxn1(40, 10, i, model = "uniform", omega = 100),
    Abarx(40, i, model = "uniform", omega = 100) -
      nAbarx(40, 10, i, model = "uniform", omega = 100)
  )
})

test_that("continuous whole life APV under uniform survival is between 0 and 1", {
  i <- 0.05
  val <- Abarx(40, i, model = "uniform", omega = 100)

  expect_true(is.numeric(val))
  expect_true(val > 0)
  expect_true(val < 1)
})

test_that("continuous term insurance APV is not greater than continuous whole life APV", {
  i <- 0.05

  term_val <- Abarxn1(40, 10, i, model = "uniform", omega = 100)
  whole_val <- Abarx(40, i, model = "uniform", omega = 100)

  expect_true(term_val <= whole_val)
})

test_that("continuous endowment APV is at least the pure endowment component", {
  i <- 0.05
  pure_val <- discount(i, 10) * tpx(10, x = 40, model = "uniform", omega = 100)
  endow_val <- Abarxn(40, 10, i, model = "uniform", omega = 100)

  expect_true(endow_val >= pure_val)
})
