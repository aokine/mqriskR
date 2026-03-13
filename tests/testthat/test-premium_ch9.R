test_that("Px matches Chapter 9 whole life premium identity", {
  x <- 40
  i <- 0.05

  expect_equal(
    Px(x, i = i, model = "uniform", omega = 100),
    Ax(x, i = i, model = "uniform", omega = 100) /
      adotx(x, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  d <- i / (1 + i)

  expect_equal(
    Px(x, i = i, model = "uniform", omega = 100),
    1 / adotx(x, i = i, model = "uniform", omega = 100) - d,
    tolerance = 1e-10
  )
})

test_that("Pxn1, PnEx, and Pxn match Chapter 9 temporary formulas", {
  x <- 40
  n <- 20
  i <- 0.05

  p_term <- Pxn1(x, n, i = i, model = "uniform", omega = 100)
  p_pe   <- PnEx(x, n, i = i, model = "uniform", omega = 100)
  p_end  <- Pxn(x, n, i = i, model = "uniform", omega = 100)

  expect_equal(
    p_term,
    Axn1(x, n, i = i, model = "uniform", omega = 100) /
      adotxn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    p_pe,
    nEx(x, n, i = i, model = "uniform", omega = 100) /
      adotxn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    p_end,
    Axn(x, n, i = i, model = "uniform", omega = 100) /
      adotxn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(p_end, p_term + p_pe, tolerance = 1e-10)
})

test_that("limited-payment premiums use temporary premium annuity denominator", {
  x <- 35
  n <- 20
  t <- 10
  i <- 0.05

  expect_equal(
    tPx(x, t, i = i, model = "uniform", omega = 100),
    Ax(x, i = i, model = "uniform", omega = 100) /
      adotxn(x, t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    tPxn1(x, n, t, i = i, model = "uniform", omega = 100),
    Axn1(x, n, i = i, model = "uniform", omega = 100) /
      adotxn(x, t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    tPnEx(x, n, t, i = i, model = "uniform", omega = 100),
    nEx(x, n, i = i, model = "uniform", omega = 100) /
      adotxn(x, t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    tPxn(x, n, t, i = i, model = "uniform", omega = 100),
    Axn(x, n, i = i, model = "uniform", omega = 100) /
      adotxn(x, t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("deferred insurance premium functions match Chapter 9 identities", {
  x <- 40
  n <- 15
  t <- 10
  i <- 0.05

  expect_equal(
    PnAx(x, n, i = i, model = "uniform", omega = 100),
    nAx(x, n, i = i, model = "uniform", omega = 100) /
      adotx(x, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    tPnAx(x, n, t, i = i, model = "uniform", omega = 100),
    nAx(x, n, i = i, model = "uniform", omega = 100) /
      adotxn(x, t, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("continuous payment premium functions match formulas", {
  x <- 40
  n <- 20
  i <- 0.05

  expect_equal(
    Pbarx(x, i = i, model = "uniform", omega = 100),
    Ax(x, i = i, model = "uniform", omega = 100) /
      abarx(x, i = i, model = "uniform", omega = 100),
    tolerance = 1e-9
  )

  expect_equal(
    Pbarxn1(x, n, i = i, model = "uniform", omega = 100),
    Axn1(x, n, i = i, model = "uniform", omega = 100) /
      abarxn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-9
  )

  expect_equal(
    Pbarxn(x, n, i = i, model = "uniform", omega = 100),
    Axn(x, n, i = i, model = "uniform", omega = 100) /
      abarxn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-9
  )
})

test_that("fully continuous premium functions match formulas and identity", {
  x <- 40
  n <- 20
  i <- 0.05
  delta <- log(1 + i)

  expect_equal(
    PbarAbarx(x, i = i, model = "uniform", omega = 100),
    Abarx(x, i = i, model = "uniform", omega = 100) /
      abarx(x, i = i, model = "uniform", omega = 100),
    tolerance = 1e-9
  )

  expect_equal(
    PbarAbarx(x, i = i, model = "uniform", omega = 100),
    1 / abarx(x, i = i, model = "uniform", omega = 100) - delta,
    tolerance = 1e-9
  )

  expect_equal(
    PbarAbarxn1(x, n, i = i, model = "uniform", omega = 100),
    Abarxn1(x, n, i = i, model = "uniform", omega = 100) /
      abarxn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-9
  )

  expect_equal(
    PbarAbarxn(x, n, i = i, model = "uniform", omega = 100),
    Abarxn(x, n, i = i, model = "uniform", omega = 100) /
      abarxn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-9
  )
})

test_that("true fractional premiums match Chapter 9 formulas", {
  x <- 40
  n <- 20
  i <- 0.05
  m <- 12

  expect_equal(
    Px_m(x, m = m, i = i, model = "uniform", omega = 100),
    Ax(x, i = i, model = "uniform", omega = 100) /
      adotx_m(x, m = m, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    Pxn1_m(x, n, m = m, i = i, model = "uniform", omega = 100),
    Axn1(x, n, i = i, model = "uniform", omega = 100) /
      adotxn_m(x, n, m = m, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    Pxn_m(x, n, m = m, i = i, model = "uniform", omega = 100),
    Axn(x, n, i = i, model = "uniform", omega = 100) /
      adotxn_m(x, n, m = m, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    PnAx_m(x, n, m = m, i = i, model = "uniform", omega = 100),
    nAx(x, n, i = i, model = "uniform", omega = 100) /
      adotx_m(x, m = m, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )
})

test_that("expected loss is zero under equivalence-principle premium", {
  x <- 40
  n <- 20
  i <- 0.05

  expect_equal(
    EL0x(x, P = Px(x, i = i, model = "uniform", omega = 100),
         i = i, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )

  expect_equal(
    EL0xn1(x, n, P = Pxn1(x, n, i = i, model = "uniform", omega = 100),
           i = i, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )

  expect_equal(
    EL0xn(x, n, P = Pxn(x, n, i = i, model = "uniform", omega = 100),
          i = i, model = "uniform", omega = 100),
    0,
    tolerance = 1e-10
  )

  expect_equal(
    EL0barAbarx(x, P = PbarAbarx(x, i = i, model = "uniform", omega = 100),
                i = i, model = "uniform", omega = 100),
    0,
    tolerance = 1e-9
  )
})

test_that("loss variances match Chapter 9 formulas", {
  x <- 40
  n <- 20
  i <- 0.05
  d <- i / (1 + i)
  delta <- log(1 + i)

  P_wl <- Px(x, i = i, model = "uniform", omega = 100)
  P_term <- Pxn1(x, n, i = i, model = "uniform", omega = 100)
  P_end <- Pxn(x, n, i = i, model = "uniform", omega = 100)
  P_cont <- PbarAbarx(x, i = i, model = "uniform", omega = 100)

  expect_equal(
    varL0x(x, P = P_wl, i = i, model = "uniform", omega = 100),
    (1 + P_wl / d)^2 * var_Ax(x, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    varL0xn1(x, n, P = P_term, i = i, model = "uniform", omega = 100),
    (1 + P_term / d)^2 * var_Axn1(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    varL0xn(x, n, P = P_end, i = i, model = "uniform", omega = 100),
    (1 + P_end / d)^2 * var_Axn(x, n, i = i, model = "uniform", omega = 100),
    tolerance = 1e-10
  )

  expect_equal(
    varL0barAbarx(x, P = P_cont, i = i, model = "uniform", omega = 100),
    (1 + P_cont / delta)^2 * var_Abarx(x, i = i, model = "uniform", omega = 100),
    tolerance = 1e-9
  )
})

test_that("gross premium whole life matches Example 9.10 formula structure", {
  x <- 40
  i <- 0.05

  Ax_val <- Ax(x, i = i, model = "uniform", omega = 100)
  ax_val <- ax(x, i = i, model = "uniform", omega = 100)
  adotx_val <- adotx(x, i = i, model = "uniform", omega = 100)

  G_manual <- ((1000 + 20) * Ax_val + 10 + 2 * ax_val) /
    (adotx_val - 0.75 - 0.10 * ax_val)

  G_pkg <- Gx(
    x = x,
    i = i,
    benefit = 1000,
    first_premium_pct = 0.75,
    renewal_premium_pct = 0.10,
    first_policy_exp = 10,
    renewal_policy_exp = 2,
    settlement_exp = 20,
    model = "uniform",
    omega = 100
  )

  expect_equal(G_pkg, G_manual, tolerance = 1e-10)
})

test_that("premium functions vectorize over x and n", {
  x <- c(30, 40, 50)
  n <- c(10, 15, 20)
  t <- c(5, 10, 10)
  i <- 0.05

  expect_length(Px(x, i = i, model = "uniform", omega = 100), 3)
  expect_length(Pxn1(x, n, i = i, model = "uniform", omega = 100), 3)
  expect_length(PnEx(x, n, i = i, model = "uniform", omega = 100), 3)
  expect_length(Pxn(x, n, i = i, model = "uniform", omega = 100), 3)
  expect_length(tPxn1(x, n, t, i = i, model = "uniform", omega = 100), 3)
  expect_length(PnAx(x, n, i = i, model = "uniform", omega = 100), 3)
})

test_that("functions reject invalid t-pay structures", {
  expect_error(
    tPxn1(40, n = 10, t = 12, i = 0.05, model = "uniform", omega = 100),
    "t must satisfy t <= n"
  )

  expect_error(
    tPnAx(40, n = 10, t = 12, i = 0.05, model = "uniform", omega = 100),
    "t must satisfy t <= n"
  )
})

test_that("premium functions return finite values in standard cases", {
  i <- 0.05

  vals <- c(
    Px(40, i = i, model = "uniform", omega = 100),
    Pxn1(40, 20, i = i, model = "uniform", omega = 100),
    PnEx(40, 20, i = i, model = "uniform", omega = 100),
    Pxn(40, 20, i = i, model = "uniform", omega = 100),
    PbarAbarx(40, i = i, model = "uniform", omega = 100),
    Px_m(40, m = 12, i = i, model = "uniform", omega = 100)
  )

  expect_true(all(is.finite(vals)))
  expect_true(all(vals > 0))
})
