test_that("nEx computes v^n * npx from life table", {
  tbl <- life_table(
    x = 0:4,
    px = c(0.90, 0.80, 0.60, 0.30, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    nEx(1, 3, i, tbl = tbl),
    discount(i, 3) * npx(tbl, 1, 3)
  )
})

test_that("Axn decomposes into term insurance plus pure endowment", {
  tbl <- life_table(
    x = 0:4,
    px = c(0.90, 0.80, 0.60, 0.30, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    Axn(1, 3, i, tbl = tbl),
    Axn1(1, 3, i, tbl = tbl) + nEx(1, 3, i, tbl = tbl)
  )
})

test_that("Ax decomposes into term plus deferred insurance", {
  tbl <- life_table(
    x = 0:4,
    px = c(0.90, 0.80, 0.60, 0.30, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    Ax(1, i, tbl = tbl),
    Axn1(1, 2, i, tbl = tbl) + nAx(1, 2, i, tbl = tbl)
  )
})

test_that("nAx equals nEx times Ax at age x+n", {
  tbl <- life_table(
    x = 0:4,
    px = c(0.90, 0.80, 0.60, 0.30, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    nAx(1, 2, i, tbl = tbl),
    nEx(1, 2, i, tbl = tbl) * Ax(3, i, tbl = tbl)
  )
})

test_that("Ax satisfies recursion Ax = v qx + v px Ax+1", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  v <- discount(i, 1)

  lhs <- Ax(1, i, tbl = tbl)
  rhs <- v * nqx(tbl, 1, 1) + v * npx(tbl, 1, 1) * Ax(2, i, tbl = tbl)

  expect_equal(lhs, rhs)
})

test_that("A2x evaluates Ax at doubled force", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    A2x(1, i, tbl = tbl),
    Ax(1, double_force_i(i), tbl = tbl)
  )
})

test_that("A2xn1 evaluates Axn1 at doubled force", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    A2xn1(1, 3, i, tbl = tbl),
    Axn1(1, 3, double_force_i(i), tbl = tbl)
  )
})

test_that("A2nEx evaluates nEx at doubled force", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    A2nEx(1, 3, i, tbl = tbl),
    nEx(1, 3, double_force_i(i), tbl = tbl)
  )
})

test_that("A2nAx evaluates nAx at doubled force", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    A2nAx(1, 2, i, tbl = tbl),
    nAx(1, 2, double_force_i(i), tbl = tbl)
  )
})

test_that("A2xn decomposes into second moments of term and pure endowment", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    A2xn(1, 3, i, tbl = tbl),
    A2xn1(1, 3, i, tbl = tbl) + A2nEx(1, 3, i, tbl = tbl)
  )
})

test_that("var_Ax equals A2x - Ax^2", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    var_Ax(1, i, tbl = tbl),
    A2x(1, i, tbl = tbl) - Ax(1, i, tbl = tbl)^2
  )
})

test_that("var_Axn1 equals A2xn1 - Axn1^2", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    var_Axn1(1, 3, i, tbl = tbl),
    A2xn1(1, 3, i, tbl = tbl) - Axn1(1, 3, i, tbl = tbl)^2
  )
})

test_that("var_nEx equals A2nEx - nEx^2", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    var_nEx(1, 3, i, tbl = tbl),
    A2nEx(1, 3, i, tbl = tbl) - nEx(1, 3, i, tbl = tbl)^2
  )
})

test_that("var_nAx equals A2nAx - nAx^2", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    var_nAx(1, 2, i, tbl = tbl),
    A2nAx(1, 2, i, tbl = tbl) - nAx(1, 2, i, tbl = tbl)^2
  )
})

test_that("var_Axn equals A2xn - Axn^2", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    var_Axn(1, 3, i, tbl = tbl),
    A2xn(1, 3, i, tbl = tbl) - Axn(1, 3, i, tbl = tbl)^2
  )
})

test_that("cov_term_deferred equals -Axn1 * nAx", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    cov_term_deferred(1, 2, i, tbl = tbl),
    -Axn1(1, 2, i, tbl = tbl) * nAx(1, 2, i, tbl = tbl)
  )
})

test_that("cov_term_endow equals -Axn1 * nEx", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05
  expect_equal(
    cov_term_endow(1, 3, i, tbl = tbl),
    -Axn1(1, 3, i, tbl = tbl) * nEx(1, 3, i, tbl = tbl)
  )
})

test_that("discrete insurance functions work with parametric survival models", {
  i <- 0.05

  ax_val <- Ax(40, i, model = "uniform", omega = 100)
  term_val <- Axn1(40, 10, i, model = "uniform", omega = 100)
  endow_val <- Axn(40, 10, i, model = "uniform", omega = 100)
  pure_val <- nEx(40, 10, i, model = "uniform", omega = 100)

  expect_true(is.numeric(ax_val))
  expect_true(is.numeric(term_val))
  expect_true(is.numeric(endow_val))
  expect_true(is.numeric(pure_val))

  expect_equal(endow_val, term_val + pure_val)
})

test_that("x and n vectorize correctly", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  i <- 0.05

  expect_length(Axn1(c(1, 2), c(2, 3), i, tbl = tbl), 2)
  expect_length(nEx(c(1, 2), c(2, 3), i, tbl = tbl), 2)
  expect_length(Axn(c(1, 2), c(2, 3), i, tbl = tbl), 2)
})

test_that("invalid inputs are rejected", {
  tbl <- life_table(
    x = 0:5,
    px = c(0.90, 0.80, 0.60, 0.30, 0.20, 0.00),
    radix = 10000
  )

  expect_error(Ax(-1, 0.05, tbl = tbl))
  expect_error(Axn1(1, -2, 0.05, tbl = tbl))
  expect_error(nEx(1, 2.5, 0.05, tbl = tbl))
  expect_error(Ax(1, -1, tbl = tbl))
})
