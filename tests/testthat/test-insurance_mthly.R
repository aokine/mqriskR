test_that("Axn_m decomposes into m-thly term insurance plus pure endowment", {
  i <- 0.05
  m <- 4

  expect_equal(
    Axn_m(40, 10, i, m, model = "uniform", omega = 100),
    Axn1_m(40, 10, i, m, model = "uniform", omega = 100) +
      discount(i, 10) * tpx(10, x = 40, model = "uniform", omega = 100)
  )
})

test_that("nAx_m equals v^n * npx * Ax_m at age x+n", {
  i <- 0.05
  m <- 4

  expect_equal(
    nAx_m(40, 10, i, m, model = "uniform", omega = 100),
    discount(i, 10) *
      tpx(10, x = 40, model = "uniform", omega = 100) *
      Ax_m(50, i, m, model = "uniform", omega = 100)
  )
})

test_that("A2x_m evaluates Ax_m at doubled force", {
  i <- 0.05
  m <- 4

  expect_equal(
    A2x_m(40, i, m, model = "uniform", omega = 100),
    Ax_m(40, double_force_i(i), m, model = "uniform", omega = 100)
  )
})

test_that("A2xn1_m evaluates Axn1_m at doubled force", {
  i <- 0.05
  m <- 4

  expect_equal(
    A2xn1_m(40, 10, i, m, model = "uniform", omega = 100),
    Axn1_m(40, 10, double_force_i(i), m, model = "uniform", omega = 100)
  )
})

test_that("A2nAx_m evaluates nAx_m at doubled force", {
  i <- 0.05
  m <- 4

  expect_equal(
    A2nAx_m(40, 10, i, m, model = "uniform", omega = 100),
    nAx_m(40, 10, double_force_i(i), m, model = "uniform", omega = 100)
  )
})

test_that("A2xn_m decomposes into second moments of m-thly term and pure endowment", {
  i <- 0.05
  i2 <- double_force_i(i)
  m <- 4

  expect_equal(
    A2xn_m(40, 10, i, m, model = "uniform", omega = 100),
    A2xn1_m(40, 10, i, m, model = "uniform", omega = 100) +
      discount(i2, 10) * tpx(10, x = 40, model = "uniform", omega = 100)
  )
})

test_that("var_Ax_m equals A2x_m - Ax_m^2", {
  i <- 0.05
  m <- 4

  expect_equal(
    var_Ax_m(40, i, m, model = "uniform", omega = 100),
    A2x_m(40, i, m, model = "uniform", omega = 100) -
      Ax_m(40, i, m, model = "uniform", omega = 100)^2
  )
})

test_that("var_Axn1_m equals A2xn1_m - Axn1_m^2", {
  i <- 0.05
  m <- 4

  expect_equal(
    var_Axn1_m(40, 10, i, m, model = "uniform", omega = 100),
    A2xn1_m(40, 10, i, m, model = "uniform", omega = 100) -
      Axn1_m(40, 10, i, m, model = "uniform", omega = 100)^2
  )
})

test_that("var_nAx_m equals A2nAx_m - nAx_m^2", {
  i <- 0.05
  m <- 4

  expect_equal(
    var_nAx_m(40, 10, i, m, model = "uniform", omega = 100),
    A2nAx_m(40, 10, i, m, model = "uniform", omega = 100) -
      nAx_m(40, 10, i, m, model = "uniform", omega = 100)^2
  )
})

test_that("var_Axn_m equals A2xn_m - Axn_m^2", {
  i <- 0.05
  m <- 4

  expect_equal(
    var_Axn_m(40, 10, i, m, model = "uniform", omega = 100),
    A2xn_m(40, 10, i, m, model = "uniform", omega = 100) -
      Axn_m(40, 10, i, m, model = "uniform", omega = 100)^2
  )
})

test_that("m-thly insurance functions vectorize correctly", {
  i <- 0.05
  m <- 4

  expect_length(Axn1_m(c(40, 45), c(10, 5), i, m, model = "uniform", omega = 100), 2)
  expect_length(nAx_m(c(40, 45), c(10, 5), i, m, model = "uniform", omega = 100), 2)
  expect_length(Axn_m(c(40, 45), c(10, 5), i, m, model = "uniform", omega = 100), 2)
})

test_that("m-thly insurance functions reject invalid inputs", {
  i <- 0.05
  m <- 4

  expect_error(Ax_m(-1, i, m, model = "uniform", omega = 100))
  expect_error(Axn1_m(40, -1, i, m, model = "uniform", omega = 100))
  expect_error(nAx_m(40, -1, i, m, model = "uniform", omega = 100))
  expect_error(Ax_m(40, -1, m, model = "uniform", omega = 100))
  expect_error(Ax_m(40, i, 0, model = "uniform", omega = 100))
  expect_error(Ax_m(40, i, 2.5, model = "uniform", omega = 100))
})

test_that("m-thly whole life APV under uniform survival is between 0 and 1", {
  i <- 0.05
  m <- 4
  val <- Ax_m(40, i, m, model = "uniform", omega = 100)

  expect_true(is.numeric(val))
  expect_true(val > 0)
  expect_true(val < 1)
})

test_that("m-thly term insurance APV is not greater than m-thly whole life APV", {
  i <- 0.05
  m <- 4

  term_val <- Axn1_m(40, 10, i, m, model = "uniform", omega = 100)
  whole_val <- Ax_m(40, i, m, model = "uniform", omega = 100)

  expect_true(term_val <= whole_val)
})

test_that("m-thly endowment APV is at least the pure endowment component", {
  i <- 0.05
  m <- 4
  pure_val <- discount(i, 10) * tpx(10, x = 40, model = "uniform", omega = 100)
  endow_val <- Axn_m(40, 10, i, m, model = "uniform", omega = 100)

  expect_true(endow_val >= pure_val)
})

test_that("Axn1_m equals Ax_m minus nAx_m", {
  i <- 0.05
  m <- 4

  expect_equal(
    Axn1_m(40, 10, i, m, model = "uniform", omega = 100),
    Ax_m(40, i, m, model = "uniform", omega = 100) -
      nAx_m(40, 10, i, m, model = "uniform", omega = 100)
  )
})

test_that("m-thly APV lies between annual and continuous values under uniform survival", {
  i <- 0.05
  x <- 40

  annual_val <- Ax(x, i, model = "uniform", omega = 100)
  mthly_val  <- Ax_m(x, i, 4, model = "uniform", omega = 100)
  cont_val   <- Abarx(x, i, model = "uniform", omega = 100)

  expect_true(mthly_val >= annual_val)
  expect_true(mthly_val <= cont_val)
})

test_that("larger m gives value closer to continuous insurance", {
  i <- 0.05
  x <- 40

  val_m2  <- Ax_m(x, i, 2, model = "uniform", omega = 100)
  val_m12 <- Ax_m(x, i, 12, model = "uniform", omega = 100)
  val_cont <- Abarx(x, i, model = "uniform", omega = 100)

  expect_true(abs(val_m12 - val_cont) <= abs(val_m2 - val_cont))
})

test_that("m-thly whole life APV works for exponential model", {
  i <- 0.05
  m <- 12
  lambda <- 0.25

  val <- Ax_m(40, i, m, model = "exponential", lambda = lambda)

  expect_true(is.numeric(val))
  expect_true(length(val) == 1)
  expect_true(val > 0)
  expect_true(val < 1)
})
