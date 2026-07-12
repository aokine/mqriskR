test_that("whole life retrospective reserve equals prospective reserve", {
  x <- 40
  t <- 10
  i <- 0.05

  expect_equal(
    tVx_ret(x, t = t, i = i, model = "uniform", omega = 100),
    tVx(x, t = t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("whole life retrospective reserve is zero at issue", {
  expect_equal(
    tVx_ret(40, t = 0, i = 0.05, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )
})

test_that("term retrospective reserve equals prospective reserve", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  expect_equal(
    tVxn1_ret(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    tVxn1(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("term retrospective reserve is zero at expiry", {
  expect_equal(
    tVxn1_ret(40, n = 20, t = 20, i = 0.05, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )
})

test_that("endowment retrospective reserve equals prospective reserve", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  expect_equal(
    tVxn_ret(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    tVxn(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("endowment retrospective reserve is one at maturity", {
  expect_equal(
    tVxn_ret(40, n = 20, t = 20, i = 0.05, model = "uniform", omega = 100),
    1,
    tolerance = 1e-10
  )
})

test_that("deferred insurance reserve matches prospective formula before deferral ends", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  direct <- nAx(x = x + t, n = n - t, i = i, model = "uniform", omega = 100) -
    PnAx(x = x, n = n, i = i, model = "uniform", omega = 100) *
    adotxn(x = x + t, n = n - t, i = i, model = "uniform", omega = 100)

  expect_equal(
    tVnAx(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    direct,
    tolerance = 1e-10
  )
})

test_that("deferred insurance reserve becomes whole life reserve after deferral period", {
  x <- 40
  n <- 20
  t <- 25
  i <- 0.05

  expect_equal(
    tVnAx(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    Ax(x = x + t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("h-pay deferred insurance reserve behaves correctly in all regions", {
  x <- 40
  n <- 20
  h <- 10
  i <- 0.05

  t1 <- 5
  val1 <- htVnAx(x, n = n, h = h, t = t1, i = i, model = "uniform", omega = 100)
  rhs1 <- nAx(x = x + t1, n = n - t1, i = i, model = "uniform", omega = 100) -
    tPnAx(x = x, n = n, t = h, i = i, model = "uniform", omega = 100) *
    adotxn(x = x + t1, n = h - t1, i = i, model = "uniform", omega = 100)

  expect_equal(val1, rhs1, tolerance = 1e-10)

  t2 <- 15
  val2 <- htVnAx(x, n = n, h = h, t = t2, i = i, model = "uniform", omega = 100)
  rhs2 <- nAx(x = x + t2, n = n - t2, i = i, model = "uniform", omega = 100)

  expect_equal(val2, rhs2, tolerance = 1e-10)

  t3 <- 25
  val3 <- htVnAx(x, n = n, h = h, t = t3, i = i, model = "uniform", omega = 100)
  rhs3 <- Ax(x = x + t3, i = i, model = "uniform", omega = 100)

  expect_equal(val3, rhs3, tolerance = 1e-10)
})

test_that("deferred annuity-due premium matches definition", {
  x <- 40
  n <- 20
  i <- 0.05

  expect_equal(
    PnAdotx(x, n = n, i = i, model = "uniform", omega = 100),
    nadotx(x, n = n, i = i, model = "uniform", omega = 100) /
      adotxn(x, n = n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("deferred annuity-due reserve matches prospective formula", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  direct <- nadotx(x = x + t, n = n - t, i = i, model = "uniform", omega = 100) -
    PnAdotx(x = x, n = n, i = i, model = "uniform", omega = 100) *
    adotxn(x = x + t, n = n - t, i = i, model = "uniform", omega = 100)

  expect_equal(
    tVnAdotx(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    direct,
    tolerance = 1e-10
  )
})

test_that("deferred annuity-immediate premium matches definition", {
  x <- 40
  n <- 20
  i <- 0.05

  expect_equal(
    Pnax(x, n = n, i = i, model = "uniform", omega = 100),
    nax(x, n = n, i = i, model = "uniform", omega = 100) /
      adotxn(x, n = n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("deferred annuity-immediate reserve matches prospective formula", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  direct <- nax(x = x + t, n = n - t, i = i, model = "uniform", omega = 100) -
    Pnax(x = x, n = n, i = i, model = "uniform", omega = 100) *
    adotxn(x = x + t, n = n - t, i = i, model = "uniform", omega = 100)

  expect_equal(
    tVnax(x, n = n, t = t, i = i, model = "uniform", omega = 100),
    direct,
    tolerance = 1e-10
  )
})

test_that("continuous total, mortality, and interest gains coincide when actual equals assumed", {
  Vt <- 10
  Vt1 <- 11
  P <- 1
  delta <- 0.05
  p <- 0.99
  benefit <- 0.5
  h <- 1

  gt <- GT_cont(Vt, Vt1, P, delta_actual = delta, p_actual = p, benefit = benefit, h = h)
  gm <- GM_cont(Vt, Vt1, P, delta_assumed = delta, p_actual = p, benefit = benefit, h = h)
  gi <- GI_cont(Vt, Vt1, P, delta_actual = delta, p_assumed = p, benefit = benefit, h = h)

  expect_equal(gt, gm, tolerance = 1e-12)
  expect_equal(gt, gi, tolerance = 1e-12)
})

test_that("Thiele derivative matches direct formula", {
  expect_equal(
    thiele_dVdt(V = 900, P = 25, delta = 0.05, mu = 0.002, benefit = 1000),
    25 + 0.05 * 900 - 0.002 * (1000 - 900),
    tolerance = 1e-12
  )
})

test_that("backward Thiele step matches explicit formula", {
  V_next <- 1000
  P <- 26.96
  delta <- 0.058
  mu <- 0.002
  benefit <- 1000
  h <- 1

  expected <- (V_next - h * P + h * mu * benefit) / (1 + h * (delta + mu))

  expect_equal(
    thiele_backward_step(
      V_next = V_next,
      P = P,
      delta = delta,
      mu = mu,
      benefit = benefit,
      h = h
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("backward Thiele path returns terminal value at final time", {
  times <- seq(19, 20, by = 0.25)

  vals <- thiele_backward_path(
    times = times,
    V_terminal = 1000,
    P = 26.96,
    delta = 0.058,
    mu = 0.002,
    benefit = 1000
  )

  expect_length(vals, length(times))
  expect_equal(tail(vals, 1), 1000, tolerance = 1e-12)
})

test_that("backward Thiele path one-step case matches backward step", {
  times <- c(19, 20)

  path_vals <- thiele_backward_path(
    times = times,
    V_terminal = 1000,
    P = 26.96,
    delta = 0.058,
    mu = 0.002,
    benefit = 1000
  )

  one_step <- thiele_backward_step(
    V_next = 1000,
    P = 26.96,
    delta = 0.058,
    mu = 0.002,
    benefit = 1000,
    h = 1
  )

  expect_equal(path_vals[1], one_step, tolerance = 1e-12)
  expect_equal(path_vals[2], 1000, tolerance = 1e-12)
})

test_that("extended reserve functions vectorize cleanly", {
  out1 <- tVx_ret(c(40, 45), t = c(5, 10), i = 0.05, model = "uniform", omega = 100)
  out2 <- tVnAx(c(40, 45), n = c(20, 25), t = c(10, 15), i = 0.05, model = "uniform", omega = 100)
  out3 <- tVnAdotx(c(40, 45), n = c(20, 25), t = c(5, 10), i = 0.05, model = "uniform", omega = 100)
  out4 <- tVnax(c(40, 45), n = c(20, 25), t = c(5, 10), i = 0.05, model = "uniform", omega = 100)

  expect_length(out1, 2)
  expect_length(out2, 2)
  expect_length(out3, 2)
  expect_length(out4, 2)
})

test_that("extended reserve functions reject invalid inputs", {
  expect_error(
    tVnAdotx(40, n = 10, t = 10, i = 0.05, model = "uniform", omega = 100),
    "t must satisfy t < n"
  )

  expect_error(
    tVnax(40, n = 10, t = 12, i = 0.05, model = "uniform", omega = 100),
    "t must satisfy t < n"
  )

  expect_error(
    htVnAx(40, n = 10, h = 12, t = 5, i = 0.05, model = "uniform", omega = 100),
    "h must satisfy h <= n"
  )

  expect_error(
    thiele_backward_path(c(20, 19), V_terminal = 1000, P = 10, delta = 0.05, mu = 0.002),
    "strictly increasing"
  )
})


testthat::test_that("extended reserve functions work with life tables", {
  tbl <- life_table(
    x = 0:5,
    lx = c(100000, 85000, 65000, 40000, 15000, 0)
  )

  whole_life <- tVx_ret(
    x = 0,
    t = 1,
    i = 0.05,
    tbl = tbl
  )

  deferred <- tVnAx(
    x = 0,
    n = 3,
    t = 1,
    i = 0.05,
    tbl = tbl
  )

  testthat::expect_length(whole_life, 1)
  testthat::expect_length(deferred, 1)
  testthat::expect_true(is.finite(whole_life))
  testthat::expect_true(is.finite(deferred))
})


testthat::test_that("retrospective reserve functions vectorize over interest rates", {
  out <- tVxn_ret(
    x = 40,
    n = 20,
    t = 10,
    i = c(0.03, 0.05),
    model = "uniform",
    omega = 100
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("retrospective reserves retain maturity values", {
  testthat::expect_equal(
    tVxn_ret(
      x = 40,
      n = 20,
      t = 20,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    1
  )

  testthat::expect_equal(
    tVxn1_ret(
      x = 40,
      n = 20,
      t = 20,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    0
  )
})


testthat::test_that("extended reserve functions reject incompatible lengths", {
  testthat::expect_error(
    tVnAx(
      x = c(40, 45),
      n = c(10, 15, 20),
      t = 5,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "length 1 or the common length"
  )

  testthat::expect_error(
    htVnAx(
      x = c(40, 45),
      n = 20,
      h = c(5, 10, 15),
      t = 2,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "length 1 or the common length"
  )
})


testthat::test_that("extended reserve functions require one mortality basis", {
  testthat::expect_error(
    tVx_ret(
      x = 40,
      t = 10,
      i = 0.05
    ),
    "Supply either tbl or model"
  )

  tbl <- life_table(
    x = 0:3,
    lx = c(100000, 70000, 30000, 0)
  )

  testthat::expect_error(
    tVx_ret(
      x = 0,
      t = 1,
      i = 0.05,
      model = "uniform",
      omega = 100,
      tbl = tbl
    ),
    "Supply only one"
  )
})

testthat::test_that("deferred annuity premiums require positive deferral periods", {
  testthat::expect_error(
    PnAdotx(
      x = 40,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "positive integer-like"
  )

  testthat::expect_error(
    Pnax(
      x = 40,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "positive integer-like"
  )
})


testthat::test_that("continuous gain helpers allow finite negative values", {
  out <- GT_cont(
    Vt = -1,
    Vt1 = 0.5,
    P = -0.1,
    delta_actual = -0.01,
    p_actual = 0.99,
    benefit = 0,
    h = 0.25
  )

  expected <- (-1 - 0.1 * 0.25) *
    exp(-0.01 * 0.25) -
    0.99 * 0.5

  testthat::expect_equal(out, expected)
})

testthat::test_that("Thiele derivative allows negative reserves and interest forces", {
  out <- thiele_dVdt(
    V = -10,
    P = -1,
    delta = -0.01,
    mu = 0.02,
    benefit = 100
  )

  expected <- -1 +
    (-0.01) * (-10) -
    0.02 * (100 - (-10))

  testthat::expect_equal(out, expected)
})


testthat::test_that("backward Thiele step matches its defining formula", {
  V_next <- 1000
  P <- 26.96
  delta <- 0.058
  mu <- 0.002
  benefit <- 1000
  h <- 0.25

  expected <- (
    V_next -
      h * P +
      h * mu * benefit
  ) / (
    1 + h * (delta + mu)
  )

  testthat::expect_equal(
    thiele_backward_step(
      V_next = V_next,
      P = P,
      delta = delta,
      mu = mu,
      benefit = benefit,
      h = h
    ),
    expected
  )
})


testthat::test_that("backward Thiele path preserves the terminal reserve", {
  times <- seq(19, 20, by = 0.25)

  out <- thiele_backward_path(
    times = times,
    V_terminal = 1000,
    P = 26.96,
    delta = 0.058,
    mu = 0.002,
    benefit = 1000
  )

  testthat::expect_length(out, length(times))
  testthat::expect_equal(out[length(out)], 1000)
  testthat::expect_true(all(is.finite(out)))
})

testthat::test_that("backward Thiele path checks step-vector lengths", {
  times <- seq(0, 1, by = 0.25)

  testthat::expect_error(
    thiele_backward_path(
      times = times,
      V_terminal = 1000,
      P = c(10, 20),
      delta = 0.05,
      mu = 0.01,
      benefit = 1000
    ),
    "P must have length 1 or 4"
  )

  testthat::expect_error(
    thiele_backward_path(
      times = times,
      V_terminal = c(1000, 1001),
      P = 10,
      delta = 0.05,
      mu = 0.01,
      benefit = 1000
    ),
    "V_terminal must be a finite scalar"
  )
})

