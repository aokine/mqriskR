test_that("Iax equals manual sum under uniform model", {
  x <- 40
  i <- 0.05
  v <- 1 / (1 + i)
  omega <- 100

  times <- 1:(omega - x)
  expected <- sum(times * (v^times) * (1 - times / (omega - x)))

  expect_equal(
    Iax(x, i, model = "uniform", omega = omega),
    expected,
    tolerance = 1e-10
  )
})

test_that("Iaxn equals manual finite sum", {
  x <- 40
  n <- 10
  i <- 0.05
  v <- 1 / (1 + i)
  omega <- 100

  times <- 1:n
  expected <- sum(times * (v^times) * (1 - times / (omega - x)))

  expect_equal(
    Iaxn(x, n, i, model = "uniform", omega = omega),
    expected,
    tolerance = 1e-10
  )
})

test_that("Daxn equals manual finite sum", {
  x <- 40
  n <- 10
  i <- 0.05
  v <- 1 / (1 + i)
  omega <- 100

  times <- 1:n
  expected <- sum((n + 1 - times) * (v^times) * (1 - times / (omega - x)))

  expect_equal(
    Daxn(x, n, i, model = "uniform", omega = omega),
    expected,
    tolerance = 1e-10
  )
})

test_that("Iadotx equals manual sum under uniform model", {
  x <- 40
  i <- 0.05
  v <- 1 / (1 + i)
  omega <- 100

  times <- 0:(omega - x)
  expected <- sum((times + 1) * (v^times) * (1 - times / (omega - x)))

  expect_equal(
    Iadotx(x, i, model = "uniform", omega = omega),
    expected,
    tolerance = 1e-10
  )
})

test_that("Iadotxn equals manual finite sum", {
  x <- 40
  n <- 10
  i <- 0.05
  v <- 1 / (1 + i)
  omega <- 100

  times <- 0:(n - 1)
  expected <- sum((times + 1) * (v^times) * (1 - times / (omega - x)))

  expect_equal(
    Iadotxn(x, n, i, model = "uniform", omega = omega),
    expected,
    tolerance = 1e-10
  )
})

test_that("Dadotxn equals manual finite sum", {
  x <- 40
  n <- 10
  i <- 0.05
  v <- 1 / (1 + i)
  omega <- 100

  times <- 0:(n - 1)
  expected <- sum((n - times) * (v^times) * (1 - times / (omega - x)))

  expect_equal(
    Dadotxn(x, n, i, model = "uniform", omega = omega),
    expected,
    tolerance = 1e-10
  )
})

test_that("temporary increasing and decreasing immediate annuities add to (n+1) times level annuity", {
  x <- 40
  n <- 10
  i <- 0.05

  inc <- Iaxn(x, n, i, model = "uniform", omega = 100)
  dec <- Daxn(x, n, i, model = "uniform", omega = 100)
  lvl <- axn(x, n, i, model = "uniform", omega = 100)

  expect_equal(inc + dec, (n + 1) * lvl, tolerance = 1e-10)
})

test_that("temporary increasing and decreasing due annuities add to (n+1) times level due annuity", {
  x <- 40
  n <- 10
  i <- 0.05

  inc <- Iadotxn(x, n, i, model = "uniform", omega = 100)
  dec <- Dadotxn(x, n, i, model = "uniform", omega = 100)
  lvl <- adotxn(x, n, i, model = "uniform", omega = 100)

  expect_equal(inc + dec, (n + 1) * lvl, tolerance = 1e-10)
})

test_that("continuous temporary increasing and decreasing annuities add to n times level continuous annuity", {
  x <- 40
  n <- 10
  i <- 0.05

  inc <- Iabarxn(x, n, i, model = "uniform", omega = 100)
  dec <- Dabarxn(x, n, i, model = "uniform", omega = 100)
  lvl <- abarxn(x, n, i, model = "uniform", omega = 100)

  expect_equal(inc + dec, n * lvl, tolerance = 1e-8)
})

test_that("Iabarx matches numerical integral under exponential model", {
  x <- 40
  i <- 0.05
  delta <- log(1 + i)
  lambda <- 0.02

  expected <- 1 / (delta + lambda)^2

  expect_equal(
    Iabarx(x, i, model = "exponential", lambda = lambda),
    expected,
    tolerance = 1e-8
  )
})

test_that("Iabarxn and Dabarxn are nonnegative", {
  x <- 40
  n <- 10
  i <- 0.05

  expect_true(Iabarxn(x, n, i, model = "uniform", omega = 100) >= 0)
  expect_true(Dabarxn(x, n, i, model = "uniform", omega = 100) >= 0)
})

test_that("average of increasing and decreasing temporary annuities matches level annuity identity", {
  x <- 40
  n <- 10
  i <- 0.05

  inc_immed <- Iaxn(x, n, i, model = "uniform", omega = 100)
  dec_immed <- Daxn(x, n, i, model = "uniform", omega = 100)
  lvl_immed <- axn(x, n, i, model = "uniform", omega = 100)

  inc_due <- Iadotxn(x, n, i, model = "uniform", omega = 100)
  dec_due <- Dadotxn(x, n, i, model = "uniform", omega = 100)
  lvl_due <- adotxn(x, n, i, model = "uniform", omega = 100)

  expect_equal(inc_immed + dec_immed, (n + 1) * lvl_immed, tolerance = 1e-10)
  expect_equal(inc_due + dec_due, (n + 1) * lvl_due, tolerance = 1e-10)
})

test_that("functions vectorize over x and n", {
  x <- c(40, 45)
  n <- c(10, 5)
  i <- 0.05

  expect_length(Iax(x, i, model = "uniform", omega = 100), 2)
  expect_length(Iaxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(Daxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(Iadotx(x, i, model = "uniform", omega = 100), 2)
  expect_length(Iadotxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(Dadotxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(Iabarx(x, i, model = "uniform", omega = 100), 2)
  expect_length(Iabarxn(x, n, i, model = "uniform", omega = 100), 2)
  expect_length(Dabarxn(x, n, i, model = "uniform", omega = 100), 2)
})

test_that("functions reject incompatible x and n lengths", {
  expect_error(
    Iaxn(c(40, 45), c(10, 5, 3), 0.05, model = "uniform", omega = 100),
    "same length|length 1"
  )
  expect_error(
    Daxn(c(40, 45), c(10, 5, 3), 0.05, model = "uniform", omega = 100),
    "same length|length 1"
  )
  expect_error(
    Iadotxn(c(40, 45), c(10, 5, 3), 0.05, model = "uniform", omega = 100),
    "same length|length 1"
  )
  expect_error(
    Dadotxn(c(40, 45), c(10, 5, 3), 0.05, model = "uniform", omega = 100),
    "same length|length 1"
  )
  expect_error(
    Iabarxn(c(40, 45), c(10, 5, 3), 0.05, model = "uniform", omega = 100),
    "same length|length 1"
  )
  expect_error(
    Dabarxn(c(40, 45), c(10, 5, 3), 0.05, model = "uniform", omega = 100),
    "same length|length 1"
  )
})

test_that("functions reject invalid interest", {
  expect_error(Iax(40, -1, model = "uniform", omega = 100))
  expect_error(Iadotx(40, -1, model = "uniform", omega = 100))
  expect_error(Iabarx(40, -1, model = "uniform", omega = 100))
})

test_that("annual temporary functions reject non-integer n", {
  expect_error(Iaxn(40, 10.5, 0.05, model = "uniform", omega = 100), "integer")
  expect_error(Daxn(40, 10.5, 0.05, model = "uniform", omega = 100), "integer")
  expect_error(Iadotxn(40, 10.5, 0.05, model = "uniform", omega = 100), "integer")
  expect_error(Dadotxn(40, 10.5, 0.05, model = "uniform", omega = 100), "integer")
})

test_that("varying annuity values are finite in standard cases", {
  vals <- c(
    Iax(40, 0.05, model = "uniform", omega = 100),
    Iaxn(40, 10, 0.05, model = "uniform", omega = 100),
    Daxn(40, 10, 0.05, model = "uniform", omega = 100),
    Iadotx(40, 0.05, model = "uniform", omega = 100),
    Iadotxn(40, 10, 0.05, model = "uniform", omega = 100),
    Dadotxn(40, 10, 0.05, model = "uniform", omega = 100),
    Iabarx(40, 0.05, model = "uniform", omega = 100),
    Iabarxn(40, 10, 0.05, model = "uniform", omega = 100),
    Dabarxn(40, 10, 0.05, model = "uniform", omega = 100)
  )

  expect_true(all(is.finite(vals)))
})
