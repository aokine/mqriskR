test_that("Pr_vector_disc reproduces Example 17.1", {
  V <- c(0, 5.66, 6.17, 0)
  qx <- c(.00142, .00153, .00166)

  Pr <- Pr_vector_disc(
    V = V,
    G = 95,
    i = 0.06,
    r = 0.05,
    e = 10,
    q1 = qx,
    b1 = 50000,
    pre_contract_expense = 15
  )

  expect_equal(unname(Pr["Pr0"]), -15.00, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr1"]), 8.42, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr2"]), 8.40, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr3"]), 8.61, tolerance = 1e-2)
})

test_that("Pi_signature reproduces Example 17.2", {
  Pr <- c(-15.00, 8.42, 8.40, 8.61)
  Pi <- Pi_signature(Pr, p_tau = c(.99858, .99847, .99834))

  expect_equal(unname(Pi["Pi0"]), -15.00, tolerance = 1e-2)
  expect_equal(unname(Pi["Pi1"]), 8.42, tolerance = 1e-2)
  expect_equal(unname(Pi["Pi2"]), 8.39, tolerance = 1e-2)
  expect_equal(unname(Pi["Pi3"]), 8.58, tolerance = 1e-2)
})

test_that("NPV_profit reproduces Example 17.3", {
  Pi <- c(-15.00, 8.42, 8.39, 8.58)
  val <- NPV_profit(Pi, r = 0.10)

  expect_equal(val, 6.03, tolerance = 1e-2)
})

test_that("IRR_profit reproduces Example 17.4", {
  Pi <- c(-15.00, 8.42, 8.39, 8.58)
  irr <- IRR_profit(Pi, interval = c(0, 1))

  expect_equal(100 * irr, 31.635, tolerance = 1e-2)
})

test_that("APV_gross_premiums and profit_margin reproduce Example 17.5", {
  apv_gp <- APV_gross_premiums(
    G = rep(95, 3),
    r = 0.10,
    p_tau = c(.99858, .99847, .99834)
  )

  pm <- profit_margin(NPV = 6.03, APV_GP = apv_gp)

  expect_equal(apv_gp, 259.52, tolerance = 1e-2)
  expect_equal(round(pm, 5), 0.02324)
})

test_that("NPV_partial and discounted_payback_period reproduce Example 17.6", {
  Pi <- c(-15.00, 8.42, 8.39, 8.58)
  partial <- NPV_partial(Pi, r = 0.10)
  dpp <- discounted_payback_period(Pi, r = 0.10)

  expect_equal(unname(partial["NPV(0)"]), -15.00, tolerance = 1e-2)
  expect_equal(unname(partial["NPV(1)"]), -7.35, tolerance = 1e-2)
  expect_lt(abs(unname(partial["NPV(2)"]) - (-0.42)), 0.01)
  expect_equal(unname(partial["NPV(3)"]), 6.03, tolerance = 1e-2)
  expect_equal(dpp, 3L)
})

test_that("Pr_vector_disc reproduces Example 17.7", {
  V <- c(0, 2500, 4000, 5000, 4000, 0)
  qx <- c(.015, .017, .019, .021, .024)

  Pr <- Pr_vector_disc(
    V = V,
    G = 19250,
    i = 0.06,
    r = 0,
    e = 240,
    q1 = qx,
    b1 = 1000000,
    pre_contract_expense = 5000
  )

  expect_equal(unname(Pr["Pr0"]), -5000.00, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr1"]), 2688.10, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr2"]), 1868.60, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr3"]), 485.60, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr4"]), 534.60, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr5"]), 390.60, tolerance = 1e-2)
})

test_that("Pi_signature, NPV_profit, and profit_margin reproduce Example 17.8", {
  Pr <- c(-5000.00, 2688.10, 1868.60, 485.60, 534.60, 390.60)
  p_tau <- c(.98500, .98300, .98100, .97900, .97600)

  Pi <- Pi_signature(Pr, p_tau = p_tau)
  npv <- NPV_profit(Pi, r = 0.10)
  apv_gp <- APV_gross_premiums(
    G = rep(19250, 5),
    r = 0.10,
    p_tau = p_tau
  )
  pm <- profit_margin(npv, apv_gp)

  expect_equal(unname(Pi["Pi2"]), 1840.57, tolerance = 1e-2)
  expect_equal(unname(Pi["Pi3"]), 470.19, tolerance = 1e-2)
  expect_equal(unname(Pi["Pi4"]), 507.80, tolerance = 1e-2)
  expect_equal(unname(Pi["Pi5"]), 363.22, tolerance = 1e-2)

  expect_equal(npv, -109.52, tolerance = 1e-2)
  expect_equal(apv_gp, 77855.49, tolerance = 1e-2)
  expect_equal(round(pm, 5), -0.00141)
})

test_that("V_zeroized reproduces Chapter 17 zeroized reserves", {
  Vz <- V_zeroized(
    qx = c(.015, .017, .019, .021, .024),
    i = 0.06,
    G = 19279,
    benefit = 1000000,
    r = 0,
    e = 240,
    V_terminal = 0,
    floor_zero = TRUE
  )

  expect_equal(unname(Vz["V0"]), 0.00, tolerance = 1e-2)
  expect_equal(unname(Vz["V1"]), 0.00, tolerance = 1e-2)
  expect_equal(unname(Vz["V2"]), 2679.54, tolerance = 1e-2)
  expect_equal(unname(Vz["V3"]), 4099.54, tolerance = 1e-2)
  expect_equal(unname(Vz["V4"]), 3602.51, tolerance = 1e-2)
  expect_equal(unname(Vz["V5"]), 0.00, tolerance = 1e-2)
})

test_that("Pr_vector_disc with zeroized reserves reproduces zeroized profit pattern", {
  Vz <- c(0.00, 0.00, 2679.54, 4099.54, 3602.51, 0.00)
  qx <- c(.015, .017, .019, .021, .024)

  Pr <- Pr_vector_disc(
    V = Vz,
    G = 19279,
    i = 0.06,
    r = 0,
    e = 240,
    q1 = qx,
    b1 = 1000000,
    pre_contract_expense = 5000
  )

  Pi <- Pi_signature(Pr, p_tau = c(.98500, .98300, .98100, .97900, .97600))
  npv <- NPV_profit(Pi, r = 0.10)

  expect_equal(unname(Pr["Pr1"]), 5181.34, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr2"]), 547.35, tolerance = 1e-1)
  expect_equal(unname(Pr["Pr3"]), 0.00, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr4"]), 0.00, tolerance = 1e-2)
  expect_equal(unname(Pr["Pr5"]), 0.00, tolerance = 1e-2)

  expect_equal(unname(Pi["Pi2"]), 539.14, tolerance = 1e-1)
  expect_equal(npv, 155.88, tolerance = 1e-1)
})



testthat::test_that("profit vector follows the discrete profit formula", {
  V <- c(0, 5.66, 6.17, 0)
  q1 <- c(0.00142, 0.00153, 0.00166)

  out <- Pr_vector_disc(
    V = V,
    G = 95,
    i = 0.06,
    r = 0.05,
    e = 10,
    q1 = q1,
    b1 = 50000,
    pre_contract_expense = 15
  )

  expected_year_1 <- (
    V[1] + 95 * (1 - 0.05) - 10
  ) * 1.06 - (
    50000 * q1[1] +
      V[2] * (1 - q1[1])
  )

  testthat::expect_equal(unname(out[1]), -15)
  testthat::expect_equal(unname(out[2]), expected_year_1, tolerance = 1e-12)
  testthat::expect_named(out, paste0("Pr", 0:3))
})


testthat::test_that("profit vector supports two decrements", {
  out <- Pr_vector_disc(
    V = c(0, 10, 0),
    G = 100,
    i = 0.05,
    r = 0.10,
    e = 5,
    q1 = c(0.01, 0.02),
    q2 = c(0.05, 0.04),
    b1 = 1000,
    b2 = 50,
    s1 = 10,
    s2 = 2
  )

  testthat::expect_length(out, 3)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("profit vector accepts negative interest rates above minus one", {
  out <- Pr_vector_disc(
    V = c(0, 0),
    G = 100,
    i = -0.01,
    q1 = 0.01,
    b1 = 1000
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("profit vector rejects excessive total decrement probability", {
  testthat::expect_error(
    Pr_vector_disc(
      V = c(0, 0),
      G = 100,
      i = 0.05,
      q1 = 0.70,
      q2 = 0.40,
      b1 = 1000
    ),
    "q1 \\+ q2 must not exceed 1"
  )
})


testthat::test_that("profit vector rejects invalid expense percentages", {
  testthat::expect_error(
    Pr_vector_disc(
      V = c(0, 0),
      G = 100,
      i = 0.05,
      r = 1.20,
      q1 = 0.01,
      b1 = 1000
    ),
    "values in \\[0, 1\\]"
  )
})


testthat::test_that("pre-contract expense must be scalar", {
  testthat::expect_error(
    Pr_vector_disc(
      V = c(0, 0),
      G = 100,
      i = 0.05,
      q1 = 0.01,
      b1 = 1000,
      pre_contract_expense = c(10, 20)
    ),
    "must be a numeric scalar"
  )
})


testthat::test_that("profit signature uses start-of-year persistency", {
  Pr <- c(-15, 8, 9, 10)
  p_tau <- c(0.90, 0.80, 0.70)

  out <- Pi_signature(Pr, p_tau)

  expected <- c(
    -15,
    8,
    9 * 0.90,
    10 * 0.90 * 0.80
  )

  testthat::expect_equal(unname(out), expected)
  testthat::expect_named(out, paste0("Pi", 0:3))
})


testthat::test_that("profit signature handles one-year contracts", {
  testthat::expect_equal(
    Pi_signature(
      Pr = c(-10, 25),
      p_tau = numeric(0)
    ),
    c(Pi0 = -10, Pi1 = 25)
  )

  testthat::expect_equal(
    Pi_signature(
      Pr = c(-10, 25),
      p_tau = 0.95
    ),
    c(Pi0 = -10, Pi1 = 25)
  )
})


testthat::test_that("NPV calculation matches direct discounting", {
  Pi <- c(-15, 8.42, 8.39, 8.58)
  r <- 0.10

  expected <- sum(Pi / (1 + r)^(0:3))

  testthat::expect_equal(
    NPV_profit(Pi, r),
    expected,
    tolerance = 1e-12
  )
})


testthat::test_that("NPV vectorizes over discount rates", {
  Pi <- c(-15, 8.42, 8.39, 8.58)
  rates <- c(0.08, 0.10, 0.12)

  out <- NPV_profit(Pi, rates)

  expected <- vapply(
    rates,
    function(rate) sum(Pi / (1 + rate)^(0:3)),
    numeric(1)
  )

  testthat::expect_equal(out, expected, tolerance = 1e-12)
})


testthat::test_that("NPV accepts negative rates above minus one", {
  Pi <- c(-100, 110)

  testthat::expect_equal(
    NPV_profit(Pi, r = -0.10),
    -100 + 110 / 0.90
  )

  testthat::expect_error(
    NPV_profit(Pi, r = -1),
    "greater than -1"
  )
})


testthat::test_that("partial NPV preserves vector output for scalar rate", {
  Pi <- c(-15, 8.42, 8.39, 8.58)

  out <- NPV_partial(Pi, 0.10)

  expected <- cumsum(Pi / 1.10^(0:3))

  testthat::expect_type(out, "double")
  testthat::expect_equal(unname(out), expected)
  testthat::expect_named(
    out,
    paste0("NPV(", 0:3, ")")
  )
})


testthat::test_that("partial NPV returns a matrix for vectorized rates", {
  Pi <- c(-15, 8.42, 8.39, 8.58)

  out <- NPV_partial(
    Pi,
    r = c(0.08, 0.10)
  )

  testthat::expect_true(is.matrix(out))
  testthat::expect_equal(dim(out), c(2, 4))
  testthat::expect_equal(
    colnames(out),
    paste0("NPV(", 0:3, ")")
  )
})


testthat::test_that("discounted payback period uses partial NPV", {
  Pi <- c(-15, 8.42, 8.39, 8.58)

  partial <- NPV_partial(Pi, r = 0.10)
  expected <- which(partial >= 0)[1] - 1L

  testthat::expect_equal(
    discounted_payback_period(Pi, r = 0.10),
    as.integer(expected)
  )
})


testthat::test_that("discounted payback period vectorizes over rates", {
  Pi <- c(-15, 8.42, 8.39, 8.58)

  out <- discounted_payback_period(
    Pi,
    r = c(0.05, 0.10, 0.20)
  )

  testthat::expect_length(out, 3)
  testthat::expect_true(all(is.na(out) | out >= 0))
})


testthat::test_that("IRR produces approximately zero NPV", {
  Pi <- c(-15, 8.42, 8.39, 8.58)

  irr <- IRR_profit(Pi)

  testthat::expect_equal(
    NPV_profit(Pi, irr),
    0,
    tolerance = 1e-8
  )
})


testthat::test_that("IRR can search over a negative interval", {
  Pi <- c(-100, 90)

  irr <- IRR_profit(
    Pi,
    interval = c(-0.50, 0)
  )

  testthat::expect_equal(irr, -0.10, tolerance = 1e-8)
})


testthat::test_that("IRR rejects an interval without a sign change", {
  testthat::expect_error(
    IRR_profit(
      Pi = c(10, 10),
      interval = c(0, 1)
    ),
    "does not change sign"
  )
})


testthat::test_that("gross premium APV uses start-of-year persistency", {
  G <- rep(95, 3)
  p_tau <- c(0.90, 0.80, 0.70)
  r <- 0.10

  expected <- sum(
    G *
      c(1, 0.90, 0.90 * 0.80) /
      1.10^(0:2)
  )

  testthat::expect_equal(
    APV_gross_premiums(G, r, p_tau),
    expected,
    tolerance = 1e-12
  )
})


testthat::test_that("gross premium APV handles one premium", {
  testthat::expect_equal(
    APV_gross_premiums(
      G = 100,
      r = 0.10,
      p_tau = numeric(0)
    ),
    100
  )
})


testthat::test_that("gross premium APV vectorizes over rates", {
  out <- APV_gross_premiums(
    G = rep(95, 3),
    r = c(0.08, 0.10),
    p_tau = c(0.90, 0.80)
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("profit margin vectorizes with scalar recycling", {
  out <- profit_margin(
    NPV = c(5, 10),
    APV_GP = 100
  )

  testthat::expect_equal(out, c(0.05, 0.10))
})


testthat::test_that("profit margin requires positive premium APV", {
  testthat::expect_error(
    profit_margin(
      NPV = 5,
      APV_GP = 0
    ),
    "positive values"
  )
})


testthat::test_that("zeroized reserves satisfy the backward recursion", {
  qx <- c(0.01, 0.02)
  i <- 0.05
  G <- 20
  benefit <- 1000
  e <- 2

  out <- V_zeroized(
    qx = qx,
    i = i,
    G = G,
    benefit = benefit,
    e = e,
    floor_zero = FALSE
  )

  expected_v1 <- (
    benefit * qx[2] +
      0 * (1 - qx[2])
  ) / 1.05 - (G - e)

  expected_v0 <- (
    benefit * qx[1] +
      expected_v1 * (1 - qx[1])
  ) / 1.05 - (G - e)

  testthat::expect_equal(
    unname(out),
    c(expected_v0, expected_v1, 0),
    tolerance = 1e-12
  )
})


testthat::test_that("zeroized reserve floor removes negative reserves", {
  out <- V_zeroized(
    qx = c(0.001, 0.001),
    i = 0.05,
    G = 100,
    benefit = 1000,
    floor_zero = TRUE
  )

  testthat::expect_true(all(out >= 0))
})


testthat::test_that("zeroized reserves support vectorized yearly assumptions", {
  out <- V_zeroized(
    qx = c(0.01, 0.02, 0.03),
    i = c(0.03, 0.04, 0.05),
    G = c(20, 21, 22),
    benefit = c(1000, 900, 800),
    r = c(0.05, 0.06, 0.07),
    e = c(2, 2, 3)
  )

  testthat::expect_length(out, 4)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("zeroized reserves validate logical and percentage inputs", {
  testthat::expect_error(
    V_zeroized(
      qx = 0.01,
      i = 0.05,
      G = 20,
      benefit = 1000,
      r = 1.2
    ),
    "values in \\[0, 1\\]"
  )

  testthat::expect_error(
    V_zeroized(
      qx = 0.01,
      i = 0.05,
      G = 20,
      benefit = 1000,
      floor_zero = NA
    ),
    "must be TRUE or FALSE"
  )
})


testthat::test_that("profit-analysis functions reject incompatible lengths", {
  testthat::expect_error(
    Pr_vector_disc(
      V = c(0, 1, 2, 0),
      G = c(100, 100),
      i = 0.05,
      q1 = c(0.01, 0.02, 0.03),
      b1 = 1000
    ),
    "length 1 or 3"
  )

  testthat::expect_error(
    profit_margin(
      NPV = c(1, 2),
      APV_GP = c(100, 200, 300)
    ),
    "length 1 or the common length"
  )
})
