test_that("alpha beta gamma are finite and sensible", {
  i <- 0.05
  m <- 12

  expect_true(is.finite(.alpha_m(i, m)))
  expect_true(is.finite(.beta_m(i, m)))
  expect_true(is.finite(.gamma_m(i, m)))
  expect_true(.alpha_m(i, m) > 0)
})

test_that("UDD m-thly due approximation equals formula", {
  x <- 40
  i <- 0.05
  m <- 12

  expected <- .alpha_m(i, m) * adotx(x, i, model = "uniform", omega = 100) -
    .beta_m(i, m)

  expect_equal(
    adotx_m_udd(x, m, i, model = "uniform", omega = 100),
    expected
  )
})

test_that("UDD m-thly immediate approximation equals formula", {
  x <- 40
  i <- 0.05
  m <- 12

  expected <- .alpha_m(i, m) * ax(x, i, model = "uniform", omega = 100) +
    .gamma_m(i, m)

  expect_equal(
    ax_m_udd(x, m, i, model = "uniform", omega = 100),
    expected
  )
})

test_that("UDD temporary due approximation equals formula", {
  x <- 40
  n <- 10
  i <- 0.05
  m <- 12

  ex <- nEx(x, n, i, model = "uniform", omega = 100)
  expected <- .alpha_m(i, m) * adotxn(x, n, i, model = "uniform", omega = 100) -
    .beta_m(i, m) * (1 - ex)

  expect_equal(
    adotxn_m_udd(x, n, m, i, model = "uniform", omega = 100),
    expected
  )
})

test_that("UDD temporary immediate approximation equals formula", {
  x <- 40
  n <- 10
  i <- 0.05
  m <- 12

  ex <- nEx(x, n, i, model = "uniform", omega = 100)
  expected <- .alpha_m(i, m) * axn(x, n, i, model = "uniform", omega = 100) +
    .gamma_m(i, m) * (1 - ex)

  expect_equal(
    axn_m_udd(x, n, m, i, model = "uniform", omega = 100),
    expected
  )
})

test_that("UDD deferred due approximation equals formula", {
  x <- 40
  n <- 10
  i <- 0.05
  m <- 12

  ex <- nEx(x, n, i, model = "uniform", omega = 100)
  expected <- .alpha_m(i, m) * nadotx(x, n, i, model = "uniform", omega = 100) -
    .beta_m(i, m) * ex

  expect_equal(
    nadotx_m_udd(x, n, m, i, model = "uniform", omega = 100),
    expected
  )
})

test_that("UDD deferred immediate approximation equals formula", {
  x <- 40
  n <- 10
  i <- 0.05
  m <- 12

  ex <- nEx(x, n, i, model = "uniform", omega = 100)
  expected <- .alpha_m(i, m) * nax(x, n, i, model = "uniform", omega = 100) +
    .gamma_m(i, m) * ex

  expect_equal(
    nax_m_udd(x, n, m, i, model = "uniform", omega = 100),
    expected
  )
})

test_that("UDD accumulated value approximations equal formulas", {
  x <- 40
  n <- 10
  i <- 0.05
  m <- 12
  ex <- nEx(x, n, i, model = "uniform", omega = 100)

  expect_equal(
    sdotxn_m_udd(x, n, m, i, model = "uniform", omega = 100),
    .alpha_m(i, m) * sdotxn(x, n, i, model = "uniform", omega = 100) -
      .beta_m(i, m) * (1 / ex - 1)
  )

  expect_equal(
    sxn_m_udd(x, n, m, i, model = "uniform", omega = 100),
    .alpha_m(i, m) * sxn(x, n, i, model = "uniform", omega = 100) +
      .gamma_m(i, m) * (1 / ex - 1)
  )
})

test_that("continuous UDD whole life approximation equals formula", {
  x <- 40
  i <- 0.05
  delta <- interest_convert(i = i)$delta
  d <- interest_convert(i = i)$d

  expect_equal(
    abarx_udd(x, i, model = "uniform", omega = 100),
    ((i * d) / delta^2) * adotx(x, i, model = "uniform", omega = 100) -
      (i - delta) / delta^2
  )
})

test_that("continuous UDD temporary approximation equals identity with existing package UDD insurance function", {
  x <- 40
  n <- 10
  i <- 0.05
  delta <- interest_convert(i = i)$delta

  expect_equal(
    abarxn_udd(x, n, i, model = "uniform", omega = 100),
    (1 - Abarxn_udd(x, n, i)) / delta
  )
})

test_that("continuous UDD deferred approximation equals deferred identity", {
  x <- 40
  n <- 10
  i <- 0.05

  expect_equal(
    nabarx_udd(x, n, i, model = "uniform", omega = 100),
    nEx(x, n, i, model = "uniform", omega = 100) *
      abarx_udd(x + n, i, model = "uniform", omega = 100)
  )
})

test_that("Woolhouse 2-term whole life approximations equal formulas", {
  x <- 40
  i <- 0.05
  m <- 12

  expect_equal(
    ax_m_woolhouse2(x, m, i, model = "uniform", omega = 100),
    ax(x, i, model = "uniform", omega = 100) + (m - 1) / (2 * m)
  )

  expect_equal(
    adotx_m_woolhouse2(x, m, i, model = "uniform", omega = 100),
    adotx(x, i, model = "uniform", omega = 100) - (m - 1) / (2 * m)
  )
})

test_that("Woolhouse 2-term temporary and deferred approximations equal formulas", {
  x <- 40
  n <- 10
  i <- 0.05
  m <- 12
  ex <- nEx(x, n, i, model = "uniform", omega = 100)

  expect_equal(
    nax_m_woolhouse2(x, n, m, i, model = "uniform", omega = 100),
    nax(x, n, i, model = "uniform", omega = 100) + (m - 1) / (2 * m) * ex
  )

  expect_equal(
    nadotx_m_woolhouse2(x, n, m, i, model = "uniform", omega = 100),
    nadotx(x, n, i, model = "uniform", omega = 100) - (m - 1) / (2 * m) * ex
  )

  expect_equal(
    axn_m_woolhouse2(x, n, m, i, model = "uniform", omega = 100),
    axn(x, n, i, model = "uniform", omega = 100) + (m - 1) / (2 * m) * (1 - ex)
  )

  expect_equal(
    adotxn_m_woolhouse2(x, n, m, i, model = "uniform", omega = 100),
    adotxn(x, n, i, model = "uniform", omega = 100) - (m - 1) / (2 * m) * (1 - ex)
  )
})

test_that("Woolhouse 2-term accumulated-value approximations equal formulas", {
  x <- 40
  n <- 10
  i <- 0.05
  m <- 12
  ex <- nEx(x, n, i, model = "uniform", omega = 100)

  expect_equal(
    sxn_m_woolhouse2(x, n, m, i, model = "uniform", omega = 100),
    sxn(x, n, i, model = "uniform", omega = 100) + (m - 1) / (2 * m) * (1 / ex - 1)
  )

  expect_equal(
    sdotxn_m_woolhouse2(x, n, m, i, model = "uniform", omega = 100),
    sdotxn(x, n, i, model = "uniform", omega = 100) - (m - 1) / (2 * m) * (1 / ex - 1)
  )
})

test_that("Woolhouse 2-term continuous whole life equals ddot ax minus one-half", {
  x <- 40
  i <- 0.05

  expect_equal(
    abarx_woolhouse2(x, i, model = "uniform", omega = 100),
    adotx(x, i, model = "uniform", omega = 100) - 0.5
  )
})

test_that("Woolhouse 3-term whole life approximations equal formulas", {
  x <- 40
  i <- 0.05
  m <- 12
  delta <- interest_convert(i = i)$delta
  mu <- hazard0(x, model = "uniform", omega = 100)

  expect_equal(
    ax_m_woolhouse3(x, m, i, model = "uniform", omega = 100),
    ax(x, i, model = "uniform", omega = 100) +
      (m - 1) / (2 * m) -
      ((m^2 - 1) / (12 * m^2)) * (mu + delta)
  )

  expect_equal(
    adotx_m_woolhouse3(x, m, i, model = "uniform", omega = 100),
    adotx(x, i, model = "uniform", omega = 100) -
      (m - 1) / (2 * m) -
      ((m^2 - 1) / (12 * m^2)) * (mu + delta)
  )
})

test_that("Woolhouse 3-term continuous whole life equals formula", {
  x <- 40
  i <- 0.05
  delta <- interest_convert(i = i)$delta
  mu <- hazard0(x, model = "uniform", omega = 100)

  expect_equal(
    abarx_woolhouse3(x, i, model = "uniform", omega = 100),
    adotx(x, i, model = "uniform", omega = 100) - 0.5 - (mu + delta) / 12
  )
})

test_that("approximation functions vectorize over x and n", {
  i <- 0.05
  m <- 12
  x <- c(40, 45)
  n <- c(10, 5)

  expect_length(ax_m_udd(x, m, i, model = "uniform", omega = 100), 2)
  expect_length(adotx_m_woolhouse2(x, m, i, model = "uniform", omega = 100), 2)
  expect_length(axn_m_udd(x, n, m, i, model = "uniform", omega = 100), 2)
  expect_length(nax_m_woolhouse3(x, n, m, i, model = "uniform", omega = 100), 2)
  expect_length(abarxn_udd(x, n, i, model = "uniform", omega = 100), 2)
})

test_that("approximation functions reject bad inputs", {
  expect_error(ax_m_udd(40, m = 0, i = 0.05, model = "uniform", omega = 100))
  expect_error(adotx_m_woolhouse2(40, m = 12, i = -1, model = "uniform", omega = 100))
  expect_error(
    axn_m_udd(c(40, 45), c(10, 5, 3), m = 12, i = 0.05, model = "uniform", omega = 100),
    "same length|length 1"
  )
})

test_that("approximation values are finite in standard cases", {
  i <- 0.05
  m <- 12

  vals <- c(
    ax_m_udd(40, m, i, model = "uniform", omega = 100),
    adotx_m_udd(40, m, i, model = "uniform", omega = 100),
    abarx_udd(40, i, model = "uniform", omega = 100),
    ax_m_woolhouse2(40, m, i, model = "uniform", omega = 100),
    ax_m_woolhouse3(40, m, i, model = "uniform", omega = 100)
  )

  expect_true(all(is.finite(vals)))
})

test_that("continuous approximation lies between annual immediate and annual due in standard case", {
  i <- 0.05

  a_immed <- ax(40, i, model = "uniform", omega = 100)
  a_due <- adotx(40, i, model = "uniform", omega = 100)
  a_cont <- abarx_woolhouse2(40, i, model = "uniform", omega = 100)

  expect_true(a_immed < a_cont)
  expect_true(a_cont < a_due)
})
