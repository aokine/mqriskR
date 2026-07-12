testthat::test_that("alphaF equals one-year term insurance APV", {
  testthat::expect_equal(
    alphaF(
      x = 40,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    Axn1(
      x = 40,
      n = 1,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    tolerance = 1e-12
  )
})


testthat::test_that("betaF equals whole life premium at age x plus one", {
  testthat::expect_equal(
    betaF(
      x = 40,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    Px(
      x = 41,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    tolerance = 1e-12
  )
})


testthat::test_that("FPT whole life reserve is zero at durations zero and one", {
  testthat::expect_equal(
    tVFx(
      x = 40,
      t = 0,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    0,
    tolerance = 1e-12
  )

  testthat::expect_equal(
    tVFx(
      x = 40,
      t = 1,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    0,
    tolerance = 1e-12
  )
})


testthat::test_that("FPT whole life reserve shifts to NLP reserve at age x plus one", {
  testthat::expect_equal(
    tVFx(
      x = 40,
      t = 5,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    tVx(
      x = 41,
      t = 4,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    tolerance = 1e-10
  )
})


testthat::test_that("FPT functions support life tables", {
  tbl <- life_table(
    x = 0:4,
    lx = c(100000, 80000, 50000, 20000, 0)
  )

  testthat::expect_true(
    is.finite(
      alphaF(
        x = 0,
        i = 0.05,
        tbl = tbl
      )
    )
  )

  testthat::expect_true(
    is.finite(
      betaF(
        x = 0,
        i = 0.05,
        tbl = tbl
      )
    )
  )

  testthat::expect_equal(
    tVFx(
      x = 0,
      t = 1,
      i = 0.05,
      tbl = tbl
    ),
    0
  )

  testthat::expect_true(
    is.finite(
      tVFx(
        x = 0,
        t = 2,
        i = 0.05,
        tbl = tbl
      )
    )
  )
})


testthat::test_that("fractional whole life reserve gives correct endpoints", {
  x <- 40
  t <- 10
  i <- 0.05

  reserve_t <- tVx(
    x = x,
    t = t,
    i = i,
    model = "uniform",
    omega = 100
  )

  reserve_t1 <- tVx(
    x = x,
    t = t + 1,
    i = i,
    model = "uniform",
    omega = 100
  )

  premium <- Px(
    x = x,
    i = i,
    model = "uniform",
    omega = 100
  )

  testthat::expect_equal(
    tsVx(
      x = x,
      t = t,
      s = 0,
      i = i,
      model = "uniform",
      omega = 100
    ),
    reserve_t + premium,
    tolerance = 1e-10
  )

  testthat::expect_equal(
    tsVx(
      x = x,
      t = t,
      s = 1,
      i = i,
      model = "uniform",
      omega = 100
    ),
    reserve_t1,
    tolerance = 1e-10
  )
})


testthat::test_that("mean reserve equals fractional reserve at one-half", {
  testthat::expect_equal(
    meanVx(
      x = 40,
      t = 10,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    tsVx(
      x = 40,
      t = 10,
      s = 0.5,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    tolerance = 1e-12
  )
})


testthat::test_that("fractional endowment reserve gives correct endpoints", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  testthat::expect_equal(
    tsVxn(
      x = x,
      n = n,
      t = t,
      s = 0,
      i = i,
      model = "uniform",
      omega = 100
    ),
    tVxn(
      x = x,
      n = n,
      t = t,
      i = i,
      model = "uniform",
      omega = 100
    ) +
      Pxn(
        x = x,
        n = n,
        i = i,
        model = "uniform",
        omega = 100
      ),
    tolerance = 1e-10
  )

  testthat::expect_equal(
    tsVxn(
      x = x,
      n = n,
      t = t,
      s = 1,
      i = i,
      model = "uniform",
      omega = 100
    ),
    tVxn(
      x = x,
      n = n,
      t = t + 1,
      i = i,
      model = "uniform",
      omega = 100
    ),
    tolerance = 1e-10
  )
})


testthat::test_that("fractional term reserve gives correct endpoints", {
  x <- 40
  n <- 20
  t <- 10
  i <- 0.05

  testthat::expect_equal(
    tsVxn1(
      x = x,
      n = n,
      t = t,
      s = 0,
      i = i,
      model = "uniform",
      omega = 100
    ),
    tVxn1(
      x = x,
      n = n,
      t = t,
      i = i,
      model = "uniform",
      omega = 100
    ) +
      Pxn1(
        x = x,
        n = n,
        i = i,
        model = "uniform",
        omega = 100
      ),
    tolerance = 1e-10
  )

  testthat::expect_equal(
    tsVxn1(
      x = x,
      n = n,
      t = t,
      s = 1,
      i = i,
      model = "uniform",
      omega = 100
    ),
    tVxn1(
      x = x,
      n = n,
      t = t + 1,
      i = i,
      model = "uniform",
      omega = 100
    ),
    tolerance = 1e-10
  )
})


testthat::test_that("fractional term and endowment reserves vectorize", {
  endowment <- tsVxn(
    x = c(40, 45),
    n = 20,
    t = c(5, 10),
    s = c(0.25, 0.75),
    i = 0.05,
    model = "uniform",
    omega = 100
  )

  term <- tsVxn1(
    x = c(40, 45),
    n = 20,
    t = c(5, 10),
    s = c(0.25, 0.75),
    i = 0.05,
    model = "uniform",
    omega = 100
  )

  testthat::expect_length(endowment, 2)
  testthat::expect_length(term, 2)
  testthat::expect_true(all(is.finite(endowment)))
  testthat::expect_true(all(is.finite(term)))
})


testthat::test_that("gross reserve equals net reserve plus expense reserve", {
  x <- 40
  t <- 10
  i <- 0.05
  G <- 0.04
  benefit <- 1
  renewal_premium_pct <- 0.10
  renewal_policy_exp <- 0.002
  settlement_exp <- 0.02

  gross_reserve <- tVGx(
    x = x,
    t = t,
    i = i,
    G = G,
    benefit = benefit,
    renewal_premium_pct = renewal_premium_pct,
    renewal_policy_exp = renewal_policy_exp,
    settlement_exp = settlement_exp,
    model = "uniform",
    omega = 100
  )

  expense_reserve <- tVEx(
    x = x,
    t = t,
    i = i,
    G = G,
    benefit = benefit,
    renewal_premium_pct = renewal_premium_pct,
    renewal_policy_exp = renewal_policy_exp,
    settlement_exp = settlement_exp,
    model = "uniform",
    omega = 100
  )

  net_reserve <- tVx(
    x = x,
    t = t,
    i = i,
    model = "uniform",
    omega = 100
  )

  testthat::expect_equal(
    gross_reserve,
    benefit * net_reserve + expense_reserve,
    tolerance = 1e-10
  )
})


testthat::test_that("gross reserve functions support life tables", {
  tbl <- life_table(
    x = 0:4,
    lx = c(100000, 80000, 50000, 20000, 0)
  )

  gross_reserve <- tVGx(
    x = 0,
    t = 1,
    i = 0.05,
    G = 0.20,
    tbl = tbl
  )

  expense_reserve <- tVEx(
    x = 0,
    t = 1,
    i = 0.05,
    G = 0.20,
    tbl = tbl
  )

  testthat::expect_true(is.finite(gross_reserve))
  testthat::expect_true(is.finite(expense_reserve))
})


testthat::test_that("gross total gain reproduces Example 11.9 value", {
  total_gain <- GTg_disc(
    VtG = 3950.73,
    Vt1G = 4607.07,
    G = 685,
    i_actual = 0.065,
    q_actual = 0.005,
    r_actual = 0.06,
    e_actual = 0,
    s_actual = 100,
    b = 50000
  )

  testthat::expect_equal(
    total_gain,
    58.75,
    tolerance = 0.02
  )
})


testthat::test_that("gross gain function vectorizes consistently", {
  out <- GTg_disc(
    VtG = c(0.10, 0.15),
    Vt1G = c(0.12, 0.17),
    G = 0.02,
    i_actual = c(0.04, 0.05),
    q_actual = 0.01,
    r_actual = 0.03,
    e_actual = 0,
    s_actual = 0.01,
    b = 1
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("ordered gross gain decomposition sums to total gain", {
  out <- decompGg_disc(
    VtG = 3950.73,
    Vt1G = 4607.07,
    G = 685,
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

  testthat::expect_type(out, "double")

  testthat::expect_named(
    out,
    c(
      "total_gain",
      "interest",
      "mortality",
      "expense",
      "check"
    )
  )

  testthat::expect_equal(
    unname(out["total_gain"]),
    unname(out["check"]),
    tolerance = 1e-10
  )
})


testthat::test_that("different decomposition orders preserve total gain", {
  out1 <- decompGg_disc(
    VtG = 3950.73,
    Vt1G = 4607.07,
    G = 685,
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
    G = 685,
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

  testthat::expect_equal(
    unname(out1["total_gain"]),
    unname(out2["total_gain"]),
    tolerance = 1e-12
  )

  testthat::expect_equal(
    unname(out1["total_gain"]),
    unname(out1["check"]),
    tolerance = 1e-10
  )

  testthat::expect_equal(
    unname(out2["total_gain"]),
    unname(out2["check"]),
    tolerance = 1e-10
  )
})


testthat::test_that("Chapter 11 functions reject invalid inputs", {
  testthat::expect_error(
    tsVx(
      x = 40,
      t = 10,
      s = 1.2,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "values in \\[0, 1\\]"
  )

  testthat::expect_error(
    tsVxn(
      x = 40,
      n = 10,
      t = 10,
      s = 0.5,
      i = 0.05,
      model = "uniform",
      omega = 100
    ),
    "t must satisfy t < n"
  )

  testthat::expect_error(
    decompGg_disc(
      VtG = 1,
      Vt1G = 1,
      G = 1,
      i_assumed = 0.05,
      q_assumed = 0.01,
      r_assumed = 0.03,
      e_assumed = 0,
      s_assumed = 0,
      i_actual = 0.06,
      q_actual = 0.02,
      r_actual = 0.04,
      e_actual = 0,
      s_actual = 0,
      order = c("interest", "interest", "expense")
    ),
    "exactly once"
  )
})


testthat::test_that("Chapter 11 reserve functions reject conflicting mortality bases", {
  tbl <- life_table(
    x = 0:3,
    lx = c(100000, 80000, 50000, 0)
  )

  testthat::expect_error(
    tsVx(
      x = 0,
      t = 1,
      s = 0.5,
      i = 0.05,
      tbl = tbl,
      model = "uniform",
      omega = 100
    ),
    "either tbl or model, not both"
  )
})
