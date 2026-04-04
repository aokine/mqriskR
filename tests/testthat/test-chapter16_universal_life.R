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
