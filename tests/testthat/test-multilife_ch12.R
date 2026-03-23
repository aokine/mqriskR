test_that("joint-life survival equals product of single-life survivals", {
  val <- tpxy(40, 50, t = 10, model = "uniform", omega = 100)
  rhs <- tpx(40, t = 10, model = "uniform", omega = 100) *
    tpx(50, t = 10, model = "uniform", omega = 100)

  expect_equal(val, rhs, tolerance = 1e-10)
})

test_that("joint-life and last-survivor survival identity holds", {
  lhs <- tpxybar(40, 50, t = 10, model = "uniform", omega = 100)
  rhs <- tpx(40, t = 10, model = "uniform", omega = 100) +
    tpx(50, t = 10, model = "uniform", omega = 100) -
    tpxy(40, 50, t = 10, model = "uniform", omega = 100)

  expect_equal(lhs, rhs, tolerance = 1e-10)
})

test_that("joint-life and last-survivor failure probabilities complement survival", {
  expect_equal(
    tqxy(40, 50, t = 10, model = "uniform", omega = 100),
    1 - tpxy(40, 50, t = 10, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    tqxybar(40, 50, t = 10, model = "uniform", omega = 100),
    1 - tpxybar(40, 50, t = 10, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("joint-life whole life insurance and annuity satisfy identity", {
  x <- 40
  y <- 50
  i <- 0.05
  d <- i / (1 + i)

  expect_equal(
    Axy(x, y, i = i, model = "uniform", omega = 100),
    1 - d * adotxy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("last-survivor whole life insurance identity holds", {
  x <- 40
  y <- 50
  i <- 0.05

  expect_equal(
    Axybar(x, y, i = i, model = "uniform", omega = 100),
    Ax(x, i = i, model = "uniform", omega = 100) +
      Ax(y, i = i, model = "uniform", omega = 100) -
      Axy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("last-survivor annuity-due identity holds", {
  x <- 40
  y <- 50
  i <- 0.05

  expect_equal(
    adotxybar(x, y, i = i, model = "uniform", omega = 100),
    adotx(x, i = i, model = "uniform", omega = 100) +
      adotx(y, i = i, model = "uniform", omega = 100) -
      adotxy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("joint-life temporary endowment equals term plus pure endowment", {
  x <- 40
  y <- 50
  n <- 10
  i <- 0.05

  expect_equal(
    Axyn(x, y, n = n, i = i, model = "uniform", omega = 100),
    Axyn1(x, y, n = n, i = i, model = "uniform", omega = 100) +
      nExy(x, y, n = n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("last-survivor temporary endowment equals term plus pure endowment", {
  x <- 40
  y <- 50
  n <- 10
  i <- 0.05

  expect_equal(
    Axybarn(x, y, n = n, i = i, model = "uniform", omega = 100),
    Axybarn1(x, y, n = n, i = i, model = "uniform", omega = 100) +
      nExybar(x, y, n = n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("reversionary annuity identities hold", {
  x <- 40
  y <- 50
  i <- 0.05

  expect_equal(
    ax_y(x, y, i = i, model = "uniform", omega = 100),
    ax(y, i = i, model = "uniform", omega = 100) -
      axy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )

  expect_equal(
    ay_x(x, y, i = i, model = "uniform", omega = 100),
    ax(x, i = i, model = "uniform", omega = 100) -
      axy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("joint-life survival functions vectorize", {
  out <- tpxy(
    x = c(40, 45),
    y = c(50, 55),
    t = c(5, 10),
    model = "uniform",
    omega = 100
  )

  expect_length(out, 2)
})

test_that("joint-life annuity and insurance functions vectorize", {
  out1 <- Axy(c(40, 45), c(50, 55), i = c(0.05, 0.06), model = "uniform", omega = 100)
  out2 <- adotxy(c(40, 45), c(50, 55), i = c(0.05, 0.06), model = "uniform", omega = 100)

  expect_length(out1, 2)
  expect_length(out2, 2)
})

test_that("multilife functions reject invalid inputs", {
  expect_error(
    tpxy(40, 50, t = -1, model = "uniform", omega = 100),
    "nonnegative"
  )

  expect_error(
    nExy(40, 50, n = -2, i = 0.05, model = "uniform", omega = 100),
    "nonnegative integer-like"
  )

  expect_error(
    Axy(40, 50, i = -1.2, model = "uniform", omega = 100),
    "greater than -1"
  )
})
