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


testthat::test_that(
  "contingent probabilities add to joint-life failure probability",
  {
    x <- 40
    y <- 50
    n <- 10

    first_x <- tqxy1(
      x, y, n = n,
      model = "uniform",
      omega = 100
    )

    first_y <- tqyx1(
      x, y, n = n,
      model = "uniform",
      omega = 100
    )

    testthat::expect_equal(
      first_x + first_y,
      tqxy(
        x, y, t = n,
        model = "uniform",
        omega = 100
      ),
      tolerance = 1e-8
    )
  }
)


testthat::test_that(
  "second-death probabilities satisfy their defining identities",
  {
    x <- 40
    y <- 50
    n <- 10

    qx <- 1 - tpx(
      t = n,
      x = x,
      model = "uniform",
      omega = 100
    )

    qy <- 1 - tpx(
      t = n,
      x = y,
      model = "uniform",
      omega = 100
    )

    testthat::expect_equal(
      tqxy1(
        x, y, n,
        model = "uniform",
        omega = 100
      ) +
        tqxy2(
          x, y, n,
          model = "uniform",
          omega = 100
        ),
      qx,
      tolerance = 1e-8
    )

    testthat::expect_equal(
      tqyx1(
        x, y, n,
        model = "uniform",
        omega = 100
      ) +
        tqyx2(
          x, y, n,
          model = "uniform",
          omega = 100
        ),
      qy,
      tolerance = 1e-8
    )
  }
)


testthat::test_that(
  "continuous joint-life insurance-annuity identity holds",
  {
    x <- 40
    y <- 50
    i <- 0.05
    delta <- log1p(i)

    testthat::expect_equal(
      Abarxy(
        x, y, i,
        model = "uniform",
        omega = 100
      ),
      1 - delta * abarxy(
        x, y, i,
        model = "uniform",
        omega = 100
      ),
      tolerance = 1e-8
    )
  }
)


testthat::test_that(
  "continuous contingent insurances add to joint-life insurance",
  {
    x <- 40
    y <- 50
    i <- 0.05

    testthat::expect_equal(
      Abarxy1(
        x, y, i,
        model = "uniform",
        omega = 100
      ) +
        Abaryx1(
          x, y, i,
          model = "uniform",
          omega = 100
        ),
      Abarxy(
        x, y, i,
        model = "uniform",
        omega = 100
      ),
      tolerance = 1e-8
    )
  }
)


testthat::test_that(
  "continuous multi-life functions vectorize consistently",
  {
    out <- abarxy(
      x = c(40, 45),
      y = c(50, 55),
      i = c(0.04, 0.05),
      model = "uniform",
      omega = 100
    )

    testthat::expect_length(out, 2)
    testthat::expect_true(all(is.finite(out)))
    testthat::expect_true(all(out >= 0))
  }
)


testthat::test_that(
  "continuous multi-life functions handle limiting ages",
  {
    testthat::expect_equal(
      abarxy(
        x = 100,
        y = 50,
        i = 0.05,
        model = "uniform",
        omega = 100
      ),
      0
    )

    testthat::expect_equal(
      Abarxy1(
        x = 100,
        y = 50,
        i = 0.05,
        model = "uniform",
        omega = 100
      ),
      0
    )
  }
)


testthat::test_that(
  "continuous multi-life functions reject incompatible vector lengths",
  {
    testthat::expect_error(
      abarxy(
        x = c(40, 45),
        y = c(50, 55, 60),
        i = 0.05,
        model = "uniform",
        omega = 100
      ),
      "length 1 or the common length"
    )
  }
)


testthat::test_that(
  "continuous multi-life functions reject conflicting mortality bases",
  {
    tbl <- life_table(
      x = 0:3,
      lx = c(100000, 80000, 50000, 0)
    )

    testthat::expect_error(
      abarxy(
        x = 0,
        y = 1,
        i = 0.05,
        tbl = tbl,
        model = "uniform",
        omega = 100
      ),
      "only one of tbl or model"
    )
  }
)


testthat::test_that(
  "continuous multi-life functions explain life-table limitation",
  {
    tbl <- life_table(
      x = 0:3,
      lx = c(100000, 80000, 50000, 0)
    )

    testthat::expect_error(
      abarxy(
        x = 0,
        y = 1,
        i = 0.05,
        tbl = tbl
      ),
      "require a parametric model"
    )

    testthat::expect_error(
      tqxy1(
        x = 0,
        y = 1,
        n = 1,
        tbl = tbl
      ),
      "require a parametric model"
    )
  }
)



