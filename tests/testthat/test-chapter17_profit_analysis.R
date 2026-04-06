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
