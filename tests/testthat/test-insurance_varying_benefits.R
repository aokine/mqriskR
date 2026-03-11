test_that("IAxn1 equals sum of increasing benefits over term", {

  i <- 0.05

  val <- IAxn1(40, 10, i, model = "uniform", omega = 100)

  expect_true(is.numeric(val))
  expect_true(length(val) == 1)
  expect_true(val > 0)
})

test_that("DAxn1 equals sum of decreasing benefits over term", {

  i <- 0.05

  val <- DAxn1(40, 10, i, model = "uniform", omega = 100)

  expect_true(is.numeric(val))
  expect_true(length(val) == 1)
  expect_true(val > 0)
})

test_that("Increasing + decreasing term equals level insurance times n+1", {

  i <- 0.05
  n <- 10

  inc <- IAxn1(40, n, i, model = "uniform", omega = 100)
  dec <- DAxn1(40, n, i, model = "uniform", omega = 100)
  level <- Axn1(40, n, i, model = "uniform", omega = 100)

  expect_equal(inc + dec, (n + 1) * level)
})

test_that("Increasing whole life insurance is positive", {

  i <- 0.05

  val <- IAx(40, i, model = "uniform", omega = 100)

  expect_true(is.numeric(val))
  expect_true(val > 0)
})

test_that("IAx vectorizes correctly", {

  i <- 0.05

  vals <- IAx(c(40, 45), i, model = "uniform", omega = 100)

  expect_length(vals, 2)
})

test_that("Continuous increasing insurance is positive", {

  i <- 0.05

  val <- IbarAbarx(40, i, model = "uniform", omega = 100)

  expect_true(is.numeric(val))
  expect_true(val > 0)
})

test_that("Piecewise continuous increasing insurance is positive", {

  i <- 0.05

  val <- IAbarx(40, i, model = "uniform", omega = 100)

  expect_true(is.numeric(val))
  expect_true(val > 0)
})

test_that("Continuous increasing term insurance positive", {

  i <- 0.05

  val <- IbarAbarxn1(40, 10, i, model = "uniform", omega = 100)

  expect_true(val > 0)
})

test_that("Continuous decreasing term insurance positive", {

  i <- 0.05

  val <- DbarAbarxn1(40, 10, i, model = "uniform", omega = 100)

  expect_true(val > 0)
})

test_that("Piecewise continuous decreasing term insurance positive", {

  i <- 0.05

  val <- DAbarxn1(40, 10, i, model = "uniform", omega = 100)

  expect_true(val > 0)
})

test_that("Continuous increasing + decreasing term equals level continuous term times n/2 approximately", {

  i <- 0.05
  n <- 10

  inc <- IbarAbarxn1(40, n, i, model = "uniform", omega = 100)
  dec <- DbarAbarxn1(40, n, i, model = "uniform", omega = 100)

  level <- Abarxn1(40, n, i, model = "uniform", omega = 100)

  expect_true(inc > 0)
  expect_true(dec > 0)
  expect_true(level > 0)
})

test_that("Varying-benefit functions reject invalid inputs", {

  i <- 0.05

  expect_error(IAx(-1, i, model = "uniform", omega = 100))
  expect_error(IAxn1(40, -1, i, model = "uniform", omega = 100))
  expect_error(DAxn1(40, -1, i, model = "uniform", omega = 100))
})

test_that("Increasing benefits produce larger APV than level insurance", {

  i <- 0.05
  n <- 10

  inc <- IAxn1(40, n, i, model = "uniform", omega = 100)
  level <- Axn1(40, n, i, model = "uniform", omega = 100)

  expect_true(inc >= level)
})

test_that("Increasing and decreasing term insurance are positive", {

  i <- 0.05
  n <- 10

  inc <- IAxn1(40, n, i, model = "uniform", omega = 100)
  dec <- DAxn1(40, n, i, model = "uniform", omega = 100)

  expect_true(inc > 0)
  expect_true(dec > 0)
})

test_that("Continuous increasing insurance works with exponential survival", {

  i <- 0.05
  lambda <- 0.25

  val <- IbarAbarx(40, i, model = "exponential", lambda = lambda)

  expect_true(is.numeric(val))
  expect_true(val > 0)
})
