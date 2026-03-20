test_that("alphaF equals one-year term APV", {
  expect_equal(
    alphaF(40, i = 0.05, model = "uniform", omega = 100),
    Axn1(40, n = 1, i = 0.05, model = "uniform", omega = 100),
    tolerance = 1e-12
  )
})

test_that("betaF equals premium at age x+1", {
  expect_equal(
    betaF(40, i = 0.05, model = "uniform", omega = 100),
    Px(41, i = 0.05, model = "uniform", omega = 100),
    tolerance = 1e-12
  )
})

test_that("FPT whole life reserve is zero in years 0 and 1", {
  expect_equal(
    tVFx(40, t = 0, i = 0.05, model = "uniform", omega = 100),
    0,
    tolerance = 1e-12
  )

  expect_equal(
    tVFx(40, t = 1, i = 0.05, model = "uniform", omega = 100),
    0,
    tolerance = 1e-12
  )
})

test_that("FPT whole life reserve shifts to NLP reserve at age x+1", {
  expect_equal(
    tVFx(40, t = 5, i = 0.05, model = "uniform", omega = 100),
    tVx(41, t = 4, i = 0.05, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("fractional whole life reserve gives initial reserve at s = 0", {
  x <- 40
  t <- 10
  i <- 0.05

  expect_equal(
    tsVx(x, t = t, s = 0, i = i, model = "uniform", omega = 100),
    tVx(x, t = t, i = i, model = "uniform", omega = 100) +
      Px(x, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("fractional whole life reserve gives next terminal reserve at s = 1", {
  x <- 40
  t <- 10
  i <- 0.05

  expect_equal(
    tsVx(x, t = t, s = 1, i = i, model = "uniform", omega = 100),
    tVx(x, t = t + 1, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("meanVx matches tsVx at s = 0.5", {
  expect_equal(
    meanVx(40, t = 10, i = 0.05, model = "uniform", omega = 100),
    tsVx(40, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100),
    tolerance = 1e-12
  )
})

test_that("fractional endowment reserve gives correct endpoints", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  expect_equal(
    tsVxn(x, n = n, t = t, s = 0, i = i, model = "uniform", omega = 100),
    tVxn(x, n = n, t = t, i = i, model = "uniform", omega = 100) +
      Pxn(x, n = n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    tsVxn(x, n = n, t = t, s = 1, i = i, model = "uniform", omega = 100),
    tVxn(x, n = n, t = t + 1, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("fractional term reserve gives correct endpoints", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  expect_equal(
    tsVxn1(x, n = n, t = t, s = 0, i = i, model = "uniform", omega = 100),
    tVxn1(x, n = n, t = t, i = i, model = "uniform", omega = 100) +
      Pxn1(x, n = n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    tsVxn1(x, n = n, t = t, s = 1, i = i, model = "uniform", omega = 100),
    tVxn1(x, n = n, t = t + 1, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("gross reserve equals net reserve plus expense reserve", {
  x <- 40
  t <- 10
  i <- 0.05
  G <- 0.04
  benefit <- 1
  r <- 0.10
  e <- 0.002
  s <- 0.02

  expect_equal(
    tVGx(
      x = x, t = t, i = i, G = G, benefit = benefit,
      renewal_premium_pct = r,
      renewal_policy_exp = e,
      settlement_exp = s,
      model = "uniform", omega = 100
    ),
    benefit * tVx(x, t = t, i = i, model = "uniform", omega = 100) +
      tVEx(
        x = x, t = t, i = i, G = G, benefit = benefit,
        renewal_premium_pct = r,
        renewal_policy_exp = e,
        settlement_exp = s,
        model = "uniform", omega = 100
      ),
    tolerance = 1e-10
  )
})

test_that("gross total gain reproduces Example 11.9 value", {
  gt <- GTg_disc(
    VtG = 3950.73,
    Vt1G = 4607.07,
    G = 685.00,
    i_actual = 0.065,
    q_actual = 0.005,
    r_actual = 0.06,
    e_actual = 0,
    s_actual = 100,
    b = 50000
  )

  expect_equal(gt, 58.75, tolerance = 0.02)
})

test_that("ordered gross gain decomposition sums to total gain", {
  out <- decompGg_disc(
    VtG = 3950.73,
    Vt1G = 4607.07,
    G = 685.00,
    i_assumed = 0.06,
    q_assumed = 0.00592,
    r_assumed = 0.05,
    e_assumed = 0,
    s_assumed = 300,
    i_actual = 0.065,
    q_actual = 0.005,
    r_actual = 0.06,
    e_actual = 0,
    s_actual = 100,
    b = 50000,
    order = c("interest", "mortality", "expense")
  )

  expect_equal(
    unname(out["total_gain"]),
    unname(out["check"]),
    tolerance = 1e-10
  )
})

test_that("different decomposition orders still sum to same total gain", {
  out1 <- decompGg_disc(
    VtG = 3950.73,
    Vt1G = 4607.07,
    G = 685.00,
    i_assumed = 0.06,
    q_assumed = 0.00592,
    r_assumed = 0.05,
    e_assumed = 0,
    s_assumed = 300,
    i_actual = 0.065,
    q_actual = 0.005,
    r_actual = 0.06,
    e_actual = 0,
    s_actual = 100,
    b = 50000,
    order = c("interest", "mortality", "expense")
  )

  out2 <- decompGg_disc(
    VtG = 3950.73,
    Vt1G = 4607.07,
    G = 685.00,
    i_assumed = 0.06,
    q_assumed = 0.00592,
    r_assumed = 0.05,
    e_assumed = 0,
    s_assumed = 300,
    i_actual = 0.065,
    q_actual = 0.005,
    r_actual = 0.06,
    e_actual = 0,
    s_actual = 100,
    b = 50000,
    order = c("mortality", "expense", "interest")
  )

  expect_equal(
    unname(out1["total_gain"]),
    unname(out2["total_gain"]),
    tolerance = 1e-12
  )

  expect_equal(
    unname(out2["total_gain"]),
    unname(out2["check"]),
    tolerance = 1e-10
  )
})

test_that("Chapter 11 extended functions reject invalid inputs", {
  expect_error(
    tsVx(40, t = 10, s = 1.2, i = 0.05, model = "uniform", omega = 100),
    "s must contain values in \\[0, 1\\]"
  )

  expect_error(
    tsVxn(40, n = 10, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100),
    "t must satisfy t < n"
  )

  expect_error(
    decompGg_disc(
      VtG = 1, Vt1G = 1, G = 1,
      i_assumed = 0.05, q_assumed = 0.01,
      r_assumed = 0.03, e_assumed = 0, s_assumed = 0,
      i_actual = 0.06, q_actual = 0.02,
      r_actual = 0.04, e_actual = 0, s_actual = 0,
      order = c("interest", "interest", "expense")
    ),
    "exactly once"
  )
})
