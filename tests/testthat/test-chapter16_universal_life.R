test_that("AV_path_ul_typeB reproduces Example 16.1", {
  qx <- c(.00076, .00081, .00085, .00090, .00095)
  r <- c(.75, .10, .10, .10, .10)
  e <- c(100, 20, 20, 20, 20)
  G <- rep(5000, 5)

  out <- AV_path_ul_typeB(
    G = G,
    r = r,
    e = e,
    qx = qx,
    ic = 0.03,
    iq = 0.03,
    B = 100000
  )

  expect_equal(out$AV[2], 1108.50, tolerance = 1e-2)
  expect_equal(out$AV[3], 5675.16, tolerance = 1e-2)
  expect_equal(out$AV[4], 10374.82, tolerance = 1e-2)
  expect_equal(out$AV[5], 15210.46, tolerance = 1e-2)
  expect_equal(out$AV[6], 20186.18, tolerance = 1e-2)
})

test_that("AV_path_ul_typeA reproduces Example 16.2", {
  qx <- c(.00076, .00081, .00085, .00090, .00095)
  r <- c(.75, .10, .10, .10, .10)
  e <- c(100, 20, 20, 20, 20)
  G <- rep(5000, 5)

  out <- AV_path_ul_typeA(
    G = G,
    r = r,
    e = e,
    qx = qx,
    ic = 0.03,
    iq = 0.03,
    B = 100000
  )

  expect_equal(out$AV[2], 1109.34, tolerance = 1e-2)
  expect_equal(out$AV[3], 5680.62, tolerance = 1e-2)
  expect_equal(out$AV[4], 10389.27, tolerance = 1e-2)
  expect_equal(out$AV[5], 15239.06, tolerance = 1e-2)
  expect_equal(out$AV[6], 20234.86, tolerance = 1e-2)
})

test_that("coi_ul_typeB matches Equation 16.5a", {
  val <- coi_ul_typeB(B = 100000, qx = 0.00076, iq = 0.03)
  expect_equal(val, 73.79, tolerance = 1e-2)
})

test_that("iP_eiul and i_credit_eiul reproduce Example 16.4", {
  idx <- c(1000, 1050, 1200, 1100, 950, 1060, 1150)

  raw <- iP_eiul(idx)
  cred <- i_credit_eiul(
    i_raw = raw,
    part = 1.10,
    floor = 0.01,
    cap = 0.10
  )

  expect_equal(round(100 * raw, 2), c(5.00, 14.29, -8.33, -13.64, 11.58, 8.49))
  expect_equal(round(100 * cred, 2), c(5.50, 10.00, 1.00, 1.00, 10.00, 9.34))
})

test_that("iMA_eiul reproduces the monthly-average growth in Example 16.5", {
  idx <- c(1000, 1020, 1100, 1150, 1080, 1040, 960, 1030, 1000, 1070, 1150, 1200, 1150)

  raw <- iMA_eiul(idx)

  expect_true(abs(raw - 0.0791666667) < 1e-10)
})

test_that("pxtau_ul and tpxtau_ul reproduce Example 16.6", {
  qd <- c(.001, .002, .003, .004, .005)
  qw <- c(.02, .02, .03, .04, .05)

  p <- pxtau_ul(qd = qd, qw = qw, year_end_withdrawal = TRUE)
  tp <- tpxtau_ul(qd = qd, qw = qw, year_end_withdrawal = TRUE)

  expect_equal(round(p, 5), c(.97902, .97804, .96709, .95616, .94525))
  expect_equal(round(tp, 5), c(.97902, .95752, .92601, .88541, .83694))
})

test_that("GMF_rollforward_ul, rt_ul, and Vprefloor_crvm_ul reproduce Example 16.9", {
  gmf10 <- GMF_rollforward_ul(
    GMF_prev = 140.40,
    GMP = 14.49,
    r = 0.04,
    policy_charge = 11.80,
    i = 0.03
  )

  r10 <- rt_ul(AV = 49.182, GMF = gmf10)
  pre <- Vprefloor_crvm_ul(r = r10, pvfb_minus_pvfp = 70)

  expect_equal(gmf10, 146.7857, tolerance = 1e-4)
  expect_equal(r10, 0.33506, tolerance = 1e-5)
  expect_equal(pre, 23.4542, tolerance = 1e-4)
})

test_that("ag38_reserve_ul reproduces Example 16.11", {
  out <- ag38_reserve_ul(
    basic_reserve = 10000,
    deficiency_reserve = 0,
    excess_payment = 60000,
    nsp_required = 100000,
    valuation_nsp = 150000,
    surrender_charge = 5000
  )

  expect_equal(out$prefunding_ratio, 0.60, tolerance = 1e-12)
  expect_equal(out$net_amount_additional, 84000, tolerance = 1e-12)
  expect_equal(out$reduced_deficiency_reserve, 0, tolerance = 1e-12)
  expect_equal(out$step8_reserve, 89000, tolerance = 1e-12)
  expect_equal(out$final_reserve, 89000, tolerance = 1e-12)
})

test_that("ag38_prefunding_ratio is capped at 1", {
  expect_equal(ag38_prefunding_ratio(120000, 100000), 1, tolerance = 1e-12)
})


testthat::test_that("Type B cost of insurance follows the direct formula", {
  testthat::expect_equal(
    coi_ul_typeB(
      B = 100000,
      qx = 0.00076,
      iq = 0.03
    ),
    100000 * 0.00076 / 1.03,
    tolerance = 1e-12
  )
})


testthat::test_that("Type B cost of insurance vectorizes consistently", {
  out <- coi_ul_typeB(
    B = 100000,
    qx = c(0.00076, 0.00081),
    iq = c(0.03, 0.04)
  )

  expected <- 100000 * c(0.00076, 0.00081) /
    (1 + c(0.03, 0.04))

  testthat::expect_equal(out, expected, tolerance = 1e-12)
})


testthat::test_that("Type B account value matches one-period formula", {
  out <- AV_path_ul_typeB(
    G = 5000,
    r = 0.10,
    e = 20,
    qx = 0.001,
    ic = 0.03,
    B = 100000,
    iq = 0.03,
    AV0 = 1000
  )

  coi <- 100000 * 0.001 / 1.03

  expected <- (
    1000 +
      5000 * (1 - 0.10) -
      20 -
      coi
  ) * 1.03

  testthat::expect_equal(out$AV[2], expected, tolerance = 1e-10)
})


testthat::test_that("Type A account value matches explicit formula", {
  out <- AV_path_ul_typeA(
    G = 5000,
    r = 0.10,
    e = 20,
    qx = 0.001,
    ic = 0.03,
    B = 100000,
    iq = 0.03,
    AV0 = 1000
  )

  mortality_factor <- 0.001 * 1.03 / 1.03

  expected <- (
    (1000 + 5000 * 0.90 - 20) * 1.03 -
      100000 * mortality_factor
  ) / (1 - mortality_factor)

  testthat::expect_equal(out$AV[2], expected, tolerance = 1e-10)
})


testthat::test_that("account-value paths recycle period inputs", {
  qx <- c(0.001, 0.002, 0.003)

  type_b <- AV_path_ul_typeB(
    G = 5000,
    r = 0.10,
    e = 20,
    qx = qx,
    ic = 0.03,
    B = c(100000, 95000, 90000)
  )

  type_a <- AV_path_ul_typeA(
    G = 5000,
    r = 0.10,
    e = 20,
    qx = qx,
    ic = 0.03,
    B = c(100000, 95000, 90000)
  )

  testthat::expect_equal(nrow(type_b), 4)
  testthat::expect_equal(nrow(type_a), 4)
  testthat::expect_true(all(is.finite(type_b$AV)))
  testthat::expect_true(all(is.finite(type_a$AV)))
})


testthat::test_that("credited-rate helper works with the default infinite cap", {
  raw <- c(-0.10, 0.05, 0.20)

  out <- i_credit_eiul(raw)

  testthat::expect_equal(
    out,
    c(0, 0.05, 0.20),
    tolerance = 1e-12
  )
})


testthat::test_that("credited rates apply floor cap participation and margin", {
  raw <- c(-0.10, 0.05, 0.20)

  out <- i_credit_eiul(
    raw,
    part = 1.10,
    floor = 0.01,
    cap = 0.10,
    margin = 0.005
  )

  expected <- pmin(
    0.10,
    pmax(0.01, 1.10 * raw - 0.005)
  )

  testthat::expect_equal(out, expected, tolerance = 1e-12)
})


testthat::test_that("point-to-point rates use consecutive index values", {
  index <- c(1000, 1050, 1200, 1100)

  testthat::expect_equal(
    iP_eiul(index),
    index[-1] / index[-length(index)] - 1,
    tolerance = 1e-12
  )
})


testthat::test_that("monthly-average growth uses twelve closing values", {
  index <- c(
    1000, 1020, 1100, 1150, 1080, 1040, 960,
    1030, 1000, 1070, 1150, 1200, 1150
  )

  testthat::expect_equal(
    iMA_eiul(index),
    mean(index[-1]) / index[1] - 1,
    tolerance = 1e-12
  )
})


testthat::test_that("persistency formulas are correct", {
  qd <- c(0.01, 0.02)
  qw <- c(0.03, 0.04)

  year_end <- (1 - qd) * (1 - qw)
  competing <- 1 - qd - qw

  testthat::expect_equal(
    pxtau_ul(qd, qw, year_end_withdrawal = TRUE),
    year_end
  )

  testthat::expect_equal(
    pxtau_ul(qd, qw, year_end_withdrawal = FALSE),
    competing
  )

  testthat::expect_equal(
    tpxtau_ul(qd, qw),
    cumprod(year_end)
  )
})


testthat::test_that("persistency probabilities recycle scalar inputs", {
  out <- pxtau_ul(
    qd = c(0.01, 0.02, 0.03),
    qw = 0.05
  )

  testthat::expect_equal(
    out,
    (1 - c(0.01, 0.02, 0.03)) * 0.95
  )
})


testthat::test_that("guaranteed maturity fund roll-forward vectorizes", {
  out <- GMF_rollforward_ul(
    GMF_prev = c(100, 120),
    GMP = 15,
    r = 0.04,
    policy_charge = c(10, 11),
    i = 0.03
  )

  expected <- (
    c(100, 120) +
      15 * (1 - 0.04) -
      c(10, 11)
  ) * 1.03

  testthat::expect_equal(out, expected, tolerance = 1e-12)
})


testthat::test_that("funding ratio helpers cap results at one", {
  testthat::expect_equal(
    rt_ul(
      AV = c(50, 120),
      GMF = 100
    ),
    c(0.5, 1)
  )

  testthat::expect_equal(
    ag38_prefunding_ratio(
      excess_payment = c(50, 120),
      nsp_required = 100
    ),
    c(0.5, 1)
  )
})


testthat::test_that("pre-floor reserve permits negative APV differences", {
  testthat::expect_equal(
    Vprefloor_crvm_ul(
      r = c(0.5, 0.75),
      pvfb_minus_pvfp = c(100, -20)
    ),
    c(50, -15)
  )
})


testthat::test_that("scalar AG 38 calculation preserves list output", {
  out <- ag38_reserve_ul(
    basic_reserve = 10000,
    deficiency_reserve = 0,
    excess_payment = 60000,
    nsp_required = 100000,
    valuation_nsp = 150000,
    surrender_charge = 5000
  )

  testthat::expect_type(out, "list")

  testthat::expect_named(
    out,
    c(
      "prefunding_ratio",
      "net_amount_additional",
      "reduced_deficiency_reserve",
      "step8_reserve",
      "increased_basic_reserve",
      "final_reserve"
    )
  )

  testthat::expect_equal(out$prefunding_ratio, 0.6)
  testthat::expect_equal(
    out$final_reserve,
    out$increased_basic_reserve
  )
})


testthat::test_that("vectorized AG 38 calculation returns a data frame", {
  out <- ag38_reserve_ul(
    basic_reserve = c(10000, 12000),
    deficiency_reserve = 0,
    excess_payment = c(60000, 80000),
    nsp_required = 100000,
    valuation_nsp = c(150000, 160000),
    surrender_charge = 5000
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(nrow(out), 2)

  testthat::expect_named(
    out,
    c(
      "prefunding_ratio",
      "net_amount_additional",
      "reduced_deficiency_reserve",
      "step8_reserve",
      "increased_basic_reserve",
      "final_reserve"
    )
  )
})


testthat::test_that("universal life functions reject invalid inputs", {
  testthat::expect_error(
    coi_ul_typeB(
      B = 100000,
      qx = 1.2,
      iq = 0.03
    ),
    "values in \\[0, 1\\]"
  )

  testthat::expect_error(
    coi_ul_typeB(
      B = 100000,
      qx = 0.01,
      iq = -1
    ),
    "greater than -1"
  )

  testthat::expect_error(
    AV_path_ul_typeB(
      G = 5000,
      r = 1.2,
      e = 20,
      qx = 0.001,
      ic = 0.03,
      B = 100000
    ),
    "values in \\[0, 1\\]"
  )

  testthat::expect_error(
    iP_eiul(c(0, 100)),
    "must be positive"
  )

  testthat::expect_error(
    i_credit_eiul(
      i_raw = 0.05,
      floor = 0.10,
      cap = 0.05
    ),
    "floor cannot exceed cap"
  )

  testthat::expect_error(
    pxtau_ul(
      qd = 0.70,
      qw = 0.40,
      year_end_withdrawal = FALSE
    ),
    "must not exceed 1"
  )

  testthat::expect_error(
    rt_ul(
      AV = 100,
      GMF = 0
    ),
    "positive values"
  )

  testthat::expect_error(
    ag38_prefunding_ratio(
      excess_payment = 100,
      nsp_required = 0
    ),
    "positive values"
  )
})


testthat::test_that("universal life functions reject incompatible lengths", {
  testthat::expect_error(
    coi_ul_typeB(
      B = c(100000, 110000),
      qx = c(0.001, 0.002, 0.003),
      iq = 0.03
    ),
    "length 1 or the common length"
  )

  testthat::expect_error(
    AV_path_ul_typeA(
      G = c(5000, 5000),
      r = c(0.10, 0.10, 0.10),
      e = 20,
      qx = c(0.001, 0.002, 0.003),
      ic = 0.03,
      B = 100000
    ),
    "length 1 or the common length"
  )
})


