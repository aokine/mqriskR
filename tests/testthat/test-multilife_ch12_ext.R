test_that("contingent probabilities add to joint-life failure probability", {
  val <- tqxy1(40, 50, n = 10, model = "uniform", omega = 100) +
    tqyx1(40, 50, n = 10, model = "uniform", omega = 100)

  expect_equal(
    val,
    tqxy(40, 50, t = 10, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("second-failure contingent probabilities add to last-survivor failure probability", {
  val <- tqxy2(40, 50, n = 10, model = "uniform", omega = 100) +
    tqyx2(40, 50, n = 10, model = "uniform", omega = 100)

  expect_equal(
    val,
    tqxybar(40, 50, t = 10, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("continuous joint-life insurance-annuity identity holds", {
  x <- 40
  y <- 50
  i <- 0.05
  delta <- log(1 + i)

  expect_equal(
    Abarxy(x, y, i = i, model = "uniform", omega = 100),
    1 - delta * abarxy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("continuous last-survivor identity holds", {
  x <- 40
  y <- 50
  i <- 0.05

  expect_equal(
    Abarxybar(x, y, i = i, model = "uniform", omega = 100),
    Abarx(x, i = i, model = "uniform", omega = 100) +
      Abarx(y, i = i, model = "uniform", omega = 100) -
      Abarxy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})

test_that("continuous reversionary annuity identities hold", {
  x <- 40
  y <- 50
  i <- 0.05

  expect_equal(
    abarx_y(x, y, i = i, model = "uniform", omega = 100),
    abarx(y, i = i, model = "uniform", omega = 100) -
      abarxy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )

  expect_equal(
    abary_x(x, y, i = i, model = "uniform", omega = 100),
    abarx(x, i = i, model = "uniform", omega = 100) -
      abarxy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-8
  )
})


test_that("continuous contingent insurance identities hold", {
  x <- 40
  y <- 50
  i <- 0.05

  expect_equal(
    Abarxy1(x, y, i = i, model = "uniform", omega = 100) +
      Abaryx1(x, y, i = i, model = "uniform", omega = 100),
    Abarxy(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-7
  )

  expect_equal(
    Abarxy2(x, y, i = i, model = "uniform", omega = 100) +
      Abaryx2(x, y, i = i, model = "uniform", omega = 100),
    Abarxybar(x, y, i = i, model = "uniform", omega = 100),
    tolerance = 1e-6
  )
})

test_that("multilife extension functions vectorize", {
  out1 <- abarxy(c(40, 45), c(50, 55), i = c(0.05, 0.06), model = "uniform", omega = 100)
  out2 <- Abarxy(c(40, 45), c(50, 55), i = c(0.05, 0.06), model = "uniform", omega = 100)
  out3 <- tqxy1(c(40, 45), c(50, 55), n = c(5, 10), model = "uniform", omega = 100)

  expect_length(out1, 2)
  expect_length(out2, 2)
  expect_length(out3, 2)
})

test_that("multilife extension functions reject invalid inputs", {
  expect_error(
    tqxy1(40, 50, n = -1, model = "uniform", omega = 100),
    "nonnegative"
  )

  expect_error(
    abarxy(40, 50, i = -1.2, model = "uniform", omega = 100),
    "greater than -1"
  )
})
