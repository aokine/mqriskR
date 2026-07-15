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
    nExy(
      40, 50,
      n = -2,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "n must contain nonnegative finite values"
  )

  expect_error(
    Axy(40, 50, i = -1.2, model = "uniform", omega = 100),
    "greater than -1"
  )
})



testthat::test_that("temporary joint-life annuity-immediate satisfies the correct identity", {
  x <- 40
  y <- 50
  n <- 10
  i <- 0.05

  expected <- adotxyn(
    x = x,
    y = y,
    n = n,
    i = i,
    model = "uniform",
    omega = 100
  ) - 1 + nExy(
    x = x,
    y = y,
    n = n,
    i = i,
    model = "uniform",
    omega = 100
  )

  testthat::expect_equal(
    axyn(
      x = x,
      y = y,
      n = n,
      i = i,
      model = "uniform",
      omega = 100
    ),
    expected,
    tolerance = 1e-12
  )
})


testthat::test_that("zero-term multi-life functions have correct values", {
  testthat::expect_equal(
    adotxyn(
      40, 50,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    0
  )

  testthat::expect_equal(
    axyn(
      40, 50,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    0
  )

  testthat::expect_equal(
    nExy(
      40, 50,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    1
  )

  testthat::expect_equal(
    nExybar(
      40, 50,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    1
  )

  testthat::expect_equal(
    Axyn1(
      40, 50,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    0
  )

  testthat::expect_equal(
    Axyn(
      40, 50,
      n = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    1
  )
})


testthat::test_that("multi-life functions vectorize over all main arguments", {
  out <- nExy(
    x = c(40, 45),
    y = c(50, 55),
    n = c(10, 5),
    i = c(0.03, 0.05),
    model = "uniform",
    omega = 100
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))

  annuity <- adotxyn(
    x = c(40, 45),
    y = 50,
    n = c(10, 5),
    i = c(0.03, 0.05),
    model = "uniform",
    omega = 100
  )

  testthat::expect_length(annuity, 2)
  testthat::expect_true(all(is.finite(annuity)))
})


testthat::test_that("multi-life functions reject incompatible argument lengths", {
  testthat::expect_error(
    nExy(
      x = c(40, 45),
      y = c(50, 55, 60),
      n = 10,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "length 1 or the common length"
  )

  testthat::expect_error(
    adotxyn(
      x = c(40, 45),
      y = 50,
      n = c(5, 10, 15),
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "length 1 or the common length"
  )
})


testthat::test_that("multi-life functions reject conflicting mortality bases", {
  tbl <- life_table(
    x = 0:4,
    lx = c(100000, 80000, 50000, 20000, 0)
  )

  testthat::expect_error(
    tpxy(
      x = 0,
      y = 1,
      t = 1,
      tbl = tbl,
      model = "uniform",
      omega = 100
    ),
    "only one of tbl or model"
  )

  testthat::expect_error(
    nExy(
      x = 0,
      y = 1,
      n = 1,
      i = 0.05
    ),
    "either tbl or model"
  )
})


testthat::test_that("multi-life functions work with finite life tables", {
  tbl <- life_table(
    x = 0:4,
    lx = c(100000, 80000, 50000, 20000, 0)
  )

  testthat::expect_equal(
    tpxy(
      x = 0,
      y = 1,
      t = 1,
      tbl = tbl
    ),
    0.8 * 0.625,
    tolerance = 1e-12
  )

  values <- c(
    adotxy(0, 1, i = 0.05, tbl = tbl),
    axy(0, 1, i = 0.05, tbl = tbl),
    Axy(0, 1, i = 0.05, tbl = tbl),
    adotxybar(0, 1, i = 0.05, tbl = tbl),
    axybar(0, 1, i = 0.05, tbl = tbl),
    Axybar(0, 1, i = 0.05, tbl = tbl)
  )

  testthat::expect_true(all(is.finite(values)))
})


testthat::test_that("joint and last-survivor probabilities satisfy identities", {
  px_joint <- tpxy(
    40, 50,
    t = 10,
    model = "uniform",
    omega = 100
  )

  px_last <- tpxybar(
    40, 50,
    t = 10,
    model = "uniform",
    omega = 100
  )

  px <- tpx(
    x = 40,
    t = 10,
    model = "uniform",
    omega = 100
  )

  py <- tpx(
    x = 50,
    t = 10,
    model = "uniform",
    omega = 100
  )

  testthat::expect_equal(
    px_joint,
    px * py,
    tolerance = 1e-12
  )

  testthat::expect_equal(
    px_last,
    px + py - px * py,
    tolerance = 1e-12
  )

  testthat::expect_equal(
    tqxy(40, 50, 10, model = "uniform", omega = 100),
    1 - px_joint,
    tolerance = 1e-12
  )

  testthat::expect_equal(
    tqxybar(40, 50, 10, model = "uniform", omega = 100),
    1 - px_last,
    tolerance = 1e-12
  )
})


testthat::test_that("multi-life insurance and annuity identities hold", {
  x <- 40
  y <- 50
  n <- 10
  i <- 0.05
  d <- i / (1 + i)

  testthat::expect_equal(
    Axyn1(
      x, y, n, i,
      model = "uniform",
      omega = 100
    ),
    1 -
      d * adotxyn(
        x, y, n, i,
        model = "uniform",
        omega = 100
      ) -
      nExy(
        x, y, n, i,
        model = "uniform",
        omega = 100
      ),
    tolerance = 1e-10
  )

  testthat::expect_equal(
    Axy(
      x, y, i,
      model = "uniform",
      omega = 100
    ),
    1 -
      d * adotxy(
        x, y, i,
        model = "uniform",
        omega = 100
      ),
    tolerance = 1e-10
  )

  testthat::expect_equal(
    Axybar(
      x, y, i,
      model = "uniform",
      omega = 100
    ),
    Ax(
      x, i,
      model = "uniform",
      omega = 100
    ) +
      Ax(
        y, i,
        model = "uniform",
        omega = 100
      ) -
      Axy(
        x, y, i,
        model = "uniform",
        omega = 100
      ),
    tolerance = 1e-10
  )
})
