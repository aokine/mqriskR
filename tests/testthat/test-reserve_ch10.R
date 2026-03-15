test_that("whole life reserve is zero at issue", {
  expect_equal(
    tVx(40, t = 0, i = 0.05, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )
})

test_that("whole life reserve matches alternative Chapter 10 formulas", {
  x <- 40
  t <- 10
  i <- 0.05

  v1 <- tVx(x, t = t, i = i, model = "uniform", omega = 100)

  v2 <- 1 - adotx(x + t, i = i, model = "uniform", omega = 100) /
    adotx(x, i = i, model = "uniform", omega = 100)

  v3 <- (Px(x + t, i = i, model = "uniform", omega = 100) -
           Px(x, i = i, model = "uniform", omega = 100)) *
    adotx(x + t, i = i, model = "uniform", omega = 100)

  v4 <- (Ax(x + t, i = i, model = "uniform", omega = 100) -
           Ax(x, i = i, model = "uniform", omega = 100)) /
    (1 - Ax(x, i = i, model = "uniform", omega = 100))

  expect_equal(v1, v2, tolerance = 1e-10)
  expect_equal(v1, v3, tolerance = 1e-10)
  expect_equal(v1, v4, tolerance = 1e-10)
})

test_that("term reserve is zero at expiry", {
  expect_equal(
    tVxn1(40, n = 20, t = 20, i = 0.05, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )
})

test_that("pure endowment reserve is one at maturity", {
  expect_equal(
    tVnEx(40, n = 20, t = 20, i = 0.05, model = "uniform", omega = 100),
    1,
    tolerance = 1e-10
  )
})

test_that("endowment reserve is one at maturity", {
  expect_equal(
    tVxn(40, n = 20, t = 20, i = 0.05, model = "uniform", omega = 100),
    1,
    tolerance = 1e-10
  )
})

test_that("endowment reserve equals term plus pure endowment reserve", {
  val_endow <- tVxn(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
  val_term  <- tVxn1(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
  val_pe    <- tVnEx(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)

  expect_equal(val_endow, val_term + val_pe, tolerance = 1e-10)
})

test_that("limited pay whole life switches to benefit-only reserve after premium period", {
  x <- 40
  h <- 10
  i <- 0.05

  expect_equal(
    htVx(x, h = h, t = 12, i = i, model = "uniform", omega = 100),
    Ax(x + 12, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("continuous premium reserve matches Chapter 10 identity", {
  x <- 40
  t <- 10
  i <- 0.05

  val <- tVbarAbarx(x, t = t, i = i, model = "uniform", omega = 100)
  rhs <- 1 - abarx(x + t, i = i, model = "uniform", omega = 100) /
    abarx(x, i = i, model = "uniform", omega = 100)

  expect_equal(val, rhs, tolerance = 1e-10)
})

test_that("m-thly reserve is zero at issue", {
  expect_equal(
    tVx_m(40, t = 0, m = 12, i = 0.05, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )
})

test_that("gain decomposition holds for discrete insurance", {
  Vt <- 0.085697
  Vt1 <- 0.096115
  P <- 0.008013
  i_actual <- 0.065
  i_assumed <- 0.06
  q_actual <- 0.00300
  q_assumed <- 0.00356

  gt <- GT_disc(Vt, Vt1, P, i_actual = i_actual, q_actual = q_actual, B = 1)
  gm <- GM_disc(Vt, Vt1, P, i_assumed = i_assumed, q_actual = q_actual, B = 1)
  gi <- GI_disc(Vt, Vt1, P, i_actual = i_actual, q_assumed = q_assumed, B = 1)

  expect_lt(abs((gm + gi) - gt), 5e-7)
})

test_that("reserve functions vectorize over x, n, and t", {
  out1 <- tVx(c(40, 45), t = c(5, 10), i = 0.05, model = "uniform", omega = 100)
  out2 <- tVxn1(c(40, 45), n = 20, t = c(5, 10), i = 0.05, model = "uniform", omega = 100)
  out3 <- tVxn(c(40, 45), n = c(15, 20), t = c(5, 10), i = 0.05, model = "uniform", omega = 100)

  expect_length(out1, 2)
  expect_length(out2, 2)
  expect_length(out3, 2)
})

test_that("reserve functions reject invalid durations", {
  expect_error(
    tVxn1(40, n = 10, t = 11, i = 0.05, model = "uniform", omega = 100),
    "t must satisfy t <= n"
  )

  expect_error(
    tVxn(40, n = 10, t = -1, i = 0.05, model = "uniform", omega = 100),
    "nonnegative integer-like"
  )
})


test_that("ELtx equals tVx when premium is the equivalence-principle premium", {
  x <- 40
  t <- 10
  i <- 0.05

  premium <- Px(x, i = i, model = "uniform", omega = 100)

  expect_equal(
    ELtx(x, t = t, i = i, P = premium, model = "uniform", omega = 100),
    tVx(x, t = t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("varLtx matches Chapter 10 variance formula", {
  x <- 40
  t <- 10
  i <- 0.05

  premium <- Px(x, i = i, model = "uniform", omega = 100)
  d <- i / (1 + i)

  rhs <- (1 + premium / d)^2 *
    (A2x(x + t, i = i, model = "uniform", omega = 100) -
       Ax(x + t, i = i, model = "uniform", omega = 100)^2)

  expect_equal(
    varLtx(x, t = t, i = i, P = premium, model = "uniform", omega = 100),
    rhs,
    tolerance = 1e-10
  )
})

test_that("ELtx and varLtx vectorize over x, t, and P", {
  i <- 0.05
  x <- c(40, 45)
  t <- c(5, 10)
  P <- c(
    Px(40, i = i, model = "uniform", omega = 100),
    Px(45, i = i, model = "uniform", omega = 100)
  )

  out_mean <- ELtx(x, t = t, i = i, P = P, model = "uniform", omega = 100)
  out_var  <- varLtx(x, t = t, i = i, P = P, model = "uniform", omega = 100)

  expect_length(out_mean, 2)
  expect_length(out_var, 2)
})

test_that("ELtx and varLtx reject invalid inputs", {
  premium <- Px(40, i = 0.05, model = "uniform", omega = 100)

  expect_error(
    ELtx(40, t = -1, i = 0.05, P = premium, model = "uniform", omega = 100),
    "nonnegative integer-like"
  )

  expect_error(
    varLtx(40, t = 5, i = 0.05, P = -0.01, model = "uniform", omega = 100),
    "nonnegative finite"
  )
})
