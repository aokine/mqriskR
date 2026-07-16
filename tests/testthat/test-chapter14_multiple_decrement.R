test_that("Axj_md reproduces Example 14.1 approximately", {
  q1 <- c(.02, .02, .02, .02, .02)
  q2 <- c(.03, .04, .05, .06, .00)
  q3 <- c(.00, .00, .00, .00, .98)
  qtau <- q1 + q2 + q3

  ptau <- numeric(length(qtau))
  ptau[1] <- 1
  for (k in 2:length(qtau)) {
    ptau[k] <- prod(1 - qtau[1:(k - 1)])
  }

  val <- Axj_md(qj = q1, ptau = ptau, i = 0.06, benefit = 1000)
  expect_equal(val, 75.33, tolerance = 1e-2)
})

test_that("Abarxj_md works for constant forces", {
  t <- seq(0, 200, by = 0.01)
  ptau <- exp(-0.012 * t)
  mu_ac <- rep(0.002, length(t))

  val <- Abarxj_md(
    t = t,
    ptau = ptau,
    muj = mu_ac,
    delta = 0.05,
    benefit = 2000
  )

  exact <- 2000 * 0.002 / (0.05 + 0.012)
  expect_equal(val, exact, tolerance = 1e-3)
})

test_that("AS_path reproduces Example 14.4 approximately", {
  r <- c(.05, .05, .05, .05, .05)
  e <- c(30, 30, 30, 30, 30)
  b1 <- c(1000, 1000, 1000, 1000, 1000)
  b2 <- c(50, 100, 300, 600, 0)
  b3 <- c(0, 0, 0, 0, 1000)
  q1 <- c(.02, .03, .04, .05, .06)
  q2 <- c(.30, .20, .20, .10, .00)
  p_tau <- c(.68, .77, .76, .85, .94)

  out <- AS_path(
    AS0 = 50,
    G = 200,
    r = r,
    e = e,
    b1 = b1,
    b2 = b2,
    q1 = q1,
    q2 = q2,
    p_tau = p_tau,
    i = 0.06,
    b3 = b3
  )

  expect_equal(out$AS[2], 275.90, tolerance = 1e-2)
  expect_equal(out$AS[3], 535.10, tolerance = 1e-2)
  expect_equal(out$AS[4], 837.90, tolerance = 1e-2)
  expect_equal(out$AS[5], 1115.03, tolerance = 1e-2)
  expect_equal(out$AS[6], 373.90, tolerance = 1e-2)
})

test_that("tp00_tp01_euler reproduces Example 14.10 approximately", {
  mu01 <- function(t) 0.10 * t + 0.20
  mu02 <- function(t) 0.20
  mu10 <- function(t) 0.50
  mu12 <- function(t) 0.125 * t + 0.20

  out <- tp00_tp01_euler(
    h = 0.10, n = 2.0,
    mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
  )

  i1 <- which(abs(out$t - 1.0) < 1e-12)
  i2 <- which(abs(out$t - 2.0) < 1e-12)

  expect_equal(round(out$tp00[i1], 3), 0.667)
  expect_equal(round(out$tp01[i1], 3), 0.144)
  expect_equal(round(out$tp00[i2], 3), 0.448)
  expect_equal(round(out$tp01[i2], 3), 0.186)
})

test_that("Pbar_trapz_ms reproduces Example 14.18 approximately", {
  mu01 <- function(t) 0.10 * t + 0.20
  mu02 <- function(t) 0.20
  mu10 <- function(t) 0.50
  mu12 <- function(t) 0.125 * t + 0.20

  ex1410 <- tp00_tp01_euler(
    h = 0.10, n = 2.0,
    mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
  )

  val <- Pbar_trapz_ms(
    t = ex1410$t,
    tp00 = ex1410$tp00,
    tp01 = ex1410$tp01,
    delta = 0.04,
    mu02 = mu02,
    mu12 = mu12,
    B02 = 1000,
    B12 = 1000,
    R = 1000
  )

  expect_equal(val, 446.95, tolerance = 1e-2)
})

test_that("thiele_dVdt_01 reproduces first derivatives in Example 14.12", {
  mu01 <- function(t) 0.10 * t + 0.20
  mu02 <- function(t) 0.20
  mu10 <- function(t) 0.50
  mu12 <- function(t) 0.125 * t + 0.20

  out <- thiele_dVdt_01(
    t = 2.00,
    V0 = 0,
    V1 = 0,
    delta = 0.04,
    Pbar = 446.95,
    B = 1000,
    R = 1000,
    mu01 = mu01,
    mu02 = mu02,
    mu10 = mu10,
    mu12 = mu12
  )

  expect_equal(unname(out["dV0"]), 246.95, tolerance = 1e-8)
  expect_equal(unname(out["dV1"]), -1450.00, tolerance = 1e-8)
})

test_that("thiele_path_01 reproduces Example 14.12 approximately", {
  mu01 <- function(t) 0.10 * t + 0.20
  mu02 <- function(t) 0.20
  mu10 <- function(t) 0.50
  mu12 <- function(t) 0.125 * t + 0.20

  out <- thiele_path_01(
    h = 0.10,
    n = 2.0,
    delta = 0.04,
    Pbar = 446.95,
    B = 1000,
    R = 1000,
    mu01 = mu01,
    mu02 = mu02,
    mu10 = mu10,
    mu12 = mu12
  )

  i190 <- which(abs(out$t - 1.90) < 1e-12)
  i100 <- which(abs(out$t - 1.00) < 1e-12)
  i000 <- which(abs(out$t - 0.00) < 1e-12)

  expect_equal(out$tV0[i190], -24.70, tolerance = 1e-2)
  expect_equal(out$tV1[i190], 145.00, tolerance = 1e-2)
  expect_equal(out$tV0[i100], -50.75, tolerance = 1e-2)
  expect_equal(out$tV1[i100], 917.99, tolerance = 1e-2)
  expect_equal(out$tV0[i000], -5.39, tolerance = 1e-2)
  expect_equal(out$tV1[i000], 1283.82, tolerance = 1e-2)
})

test_that("markov_nstep_prob agrees with direct matrix power", {
  P <- matrix(
    c(0.94, 0.03, 0.02, 0.01,
      0.50, 0.30, 0.18, 0.02,
      0.00, 0.00, 0.93, 0.07,
      0.00, 0.00, 0.00, 1.00),
    nrow = 4, byrow = TRUE
  )

  P3 <- P %*% P %*% P

  expect_equal(markov_nstep_prob(P, 3, 1, 1), P3[1, 1], tolerance = 1e-12)
  expect_equal(markov_nstep_prob(P, 3, 1, 3), P3[1, 3], tolerance = 1e-12)
})

test_that("gain_loss_md matches direct formula", {
  val <- gain_loss_md(
    Vt = 115.00,
    G = 16,
    r = 0,
    e = 3,
    i = 0.06,
    b1 = 1000,
    b2 = 110,
    s1 = 0,
    s2 = 0,
    q1 = 0.01,
    q2 = 0.10,
    Vt1 = 128.83
  )

  expected <- (115 + 16 - 3) * 1.06 - (1000 * 0.01 + 110 * 0.10 + (1 - 0.01 - 0.10) * 128.83)
  expect_equal(val, expected, tolerance = 1e-12)
})

test_that("gain_loss_md year-end Cause 2 version works", {
  val <- gain_loss_md(
    Vt = 115.00,
    G = 16,
    r = 0,
    e = 3,
    i = 0.06,
    b1 = 1000,
    b2 = 110,
    s1 = 0,
    s2 = 0,
    q1 = 0.01,
    q2 = 0.10,
    Vt1 = 128.83,
    year_end_cause2 = TRUE,
    q1prime = 0.01,
    q2prime = 0.10
  )

  expected <- (115 + 16 - 3) * 1.06 -
    (1000 * 0.01 + 110 * (1 - 0.01) * 0.10 + (1 - 0.01) * (1 - 0.10) * 128.83)

  expect_equal(val, expected, tolerance = 1e-12)
})




test_that("AS_path_md matches AS_path for two decrements", {

  r <- c(.05, .05, .05)
  e <- c(30, 30, 30)
  b1 <- c(1000, 1000, 1000)
  b2 <- c(50, 100, 300)
  q1 <- c(.02, .03, .04)
  q2 <- c(.30, .20, .20)
  p_tau <- c(.68, .77, .76)

  b_mat <- cbind(b1, b2)
  q_mat <- cbind(q1, q2)

  out_old <- AS_path(
    AS0 = 50,
    G = 200,
    r = r,
    e = e,
    b1 = b1,
    b2 = b2,
    q1 = q1,
    q2 = q2,
    p_tau = p_tau,
    i = 0.06
  )

  out_new <- AS_path_md(
    AS0 = 50,
    G = 200,
    r = r,
    e = e,
    b_mat = b_mat,
    q_mat = q_mat,
    p_tau = p_tau,
    i = 0.06
  )

  expect_equal(out_new$AS, out_old$AS, tolerance = 1e-10)
})


test_that("AS_path_md works with three decrement causes", {

  r <- c(.05, .05)
  e <- c(30, 30)

  b_mat <- matrix(
    c(1000, 50, 200,
      1000, 100, 300),
    nrow = 2, byrow = TRUE
  )

  q_mat <- matrix(
    c(.02, .03, .01,
      .03, .02, .02),
    nrow = 2, byrow = TRUE
  )

  p_tau <- 1 - rowSums(q_mat)

  out <- AS_path_md(
    AS0 = 50,
    G = 200,
    r = r,
    e = e,
    b_mat = b_mat,
    q_mat = q_mat,
    p_tau = p_tau,
    i = 0.06
  )

  # basic sanity checks
  expect_equal(nrow(out), 3)
  expect_true(all(is.finite(out$AS)))
})




testthat::test_that("Axj_md vectorizes over interest and benefit", {
  qj <- c(0.02, 0.03, 0.04)
  ptau <- c(1, 0.95, 0.88)

  out <- Axj_md(
    qj = qj,
    ptau = ptau,
    i = c(0.04, 0.06),
    benefit = 1000
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))
  testthat::expect_gt(out[1], out[2])
})


testthat::test_that("Axj_md rejects invalid inputs", {
  testthat::expect_error(
    Axj_md(
      qj = c(0.02, 1.2),
      ptau = c(1, 0.9),
      i = 0.05
    ),
    "\\[0, 1\\]"
  )

  testthat::expect_error(
    Axj_md(
      qj = c(0.02, 0.03),
      ptau = 1,
      i = 0.05
    ),
    "same length"
  )

  testthat::expect_error(
    Axj_md(
      qj = c(0.02, 0.03),
      ptau = c(1, 0.9),
      i = -1
    ),
    "greater than -1"
  )
})


testthat::test_that("Abarxj_md agrees with a constant-force closed form", {
  t <- seq(0, 20, by = 0.001)
  mu_total <- 0.03
  mu_j <- 0.01
  delta <- 0.04

  ptau <- exp(-mu_total * t)
  muj <- rep(mu_j, length(t))

  expected <- mu_j *
    (1 - exp(-(delta + mu_total) * 20)) /
    (delta + mu_total)

  testthat::expect_equal(
    Abarxj_md(
      t = t,
      ptau = ptau,
      muj = muj,
      delta = delta
    ),
    expected,
    tolerance = 1e-6
  )
})


testthat::test_that("Abarxj_md requires a strictly increasing grid", {
  testthat::expect_error(
    Abarxj_md(
      t = c(0, 2, 1),
      ptau = c(1, 0.9, 0.8),
      muj = c(0.01, 0.01, 0.01),
      delta = 0.04
    ),
    "strictly increasing"
  )
})


testthat::test_that("AS_path matches the general two-cause asset-share path", {
  b1 <- c(1000, 1000, 1000)
  b2 <- c(100, 105, 110)
  q1 <- c(0.01, 0.012, 0.014)
  q2 <- c(0.08, 0.07, 0.06)
  p_tau <- 1 - q1 - q2

  two_cause <- AS_path(
    AS0 = 0,
    G = 150,
    r = 0.05,
    e = 2,
    b1 = b1,
    b2 = b2,
    q1 = q1,
    q2 = q2,
    p_tau = p_tau,
    i = 0.04
  )

  general <- AS_path_md(
    AS0 = 0,
    G = 150,
    r = 0.05,
    e = 2,
    b_mat = cbind(b1, b2),
    q_mat = cbind(q1, q2),
    p_tau = p_tau,
    i = 0.04
  )

  testthat::expect_equal(two_cause, general)
})


testthat::test_that("AS_path_md supports yearly assumptions", {
  b_mat <- cbind(
    death = rep(1000, 3),
    withdrawal = c(100, 105, 110)
  )

  q_mat <- cbind(
    death = c(0.01, 0.012, 0.014),
    withdrawal = c(0.08, 0.07, 0.06)
  )

  out <- AS_path_md(
    AS0 = 0,
    G = c(150, 155, 160),
    r = c(0.05, 0.04, 0.03),
    e = 2,
    b_mat = b_mat,
    q_mat = q_mat,
    p_tau = 1 - rowSums(q_mat),
    i = c(0.03, 0.04, 0.05),
    b_surv = 0
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(out$k, 0:3)
  testthat::expect_true(all(is.finite(out$AS)))
})


testthat::test_that("AS_path_md rejects zero in-force probabilities", {
  testthat::expect_error(
    AS_path_md(
      AS0 = 0,
      G = 100,
      r = 0,
      e = 0,
      b_mat = matrix(1000, nrow = 1),
      q_mat = matrix(1, nrow = 1),
      p_tau = 0,
      i = 0.05
    ),
    "p_tau must be positive"
  )
})


testthat::test_that("Euler probability grid includes the terminal time", {
  mu01 <- function(t) 0.01
  mu02 <- function(t) 0.02
  mu10 <- function(t) 0.03
  mu12 <- function(t) 0.04

  out <- tp00_tp01_euler(
    h = 0.3,
    n = 1,
    mu01 = mu01,
    mu02 = mu02,
    mu10 = mu10,
    mu12 = mu12
  )

  testthat::expect_equal(tail(out$t, 1), 1)
  testthat::expect_equal(out$t, c(0, 0.3, 0.6, 0.9, 1))
})


testthat::test_that("Euler probability function handles zero term", {
  zero_rate <- function(t) 0

  out <- tp00_tp01_euler(
    h = 0.1,
    n = 0,
    mu01 = zero_rate,
    mu02 = zero_rate,
    mu10 = zero_rate,
    mu12 = zero_rate
  )

  testthat::expect_equal(out$t, 0)
  testthat::expect_equal(out$tp00, 1)
  testthat::expect_equal(out$tp01, 0)
  testthat::expect_equal(out$tp02, 0)
})


testthat::test_that("Pbar_trapz_ms rejects a zero premium denominator", {
  zero_rate <- function(t) 0

  testthat::expect_error(
    Pbar_trapz_ms(
      t = c(0, 1),
      tp00 = c(0, 0),
      tp01 = c(1, 1),
      delta = 0.04,
      mu02 = zero_rate,
      mu12 = zero_rate
    ),
    "denominator must be positive"
  )
})


testthat::test_that("Thiele derivative vectorizes", {
  mu01 <- function(t) rep(0.01, length(t))
  mu02 <- function(t) rep(0.02, length(t))
  mu10 <- function(t) rep(0.03, length(t))
  mu12 <- function(t) rep(0.04, length(t))

  out <- thiele_dVdt_01(
    t = c(0, 1),
    V0 = c(10, 20),
    V1 = c(30, 40),
    delta = 0.05,
    Pbar = 5,
    B = 1000,
    R = 50,
    mu01 = mu01,
    mu02 = mu02,
    mu10 = mu10,
    mu12 = mu12
  )

  testthat::expect_true(is.matrix(out))
  testthat::expect_equal(dim(out), c(2, 2))
  testthat::expect_equal(colnames(out), c("dV0", "dV1"))
})


testthat::test_that("backward Thiele path includes zero and terminal time", {
  mu01 <- function(t) 0.01
  mu02 <- function(t) 0.02
  mu10 <- function(t) 0.03
  mu12 <- function(t) 0.04

  out <- thiele_path_01(
    h = 0.3,
    n = 1,
    delta = 0.04,
    Pbar = 10,
    B = 1000,
    R = 100,
    mu01 = mu01,
    mu02 = mu02,
    mu10 = mu10,
    mu12 = mu12
  )

  testthat::expect_equal(out$t, c(0, 0.3, 0.6, 0.9, 1))
  testthat::expect_equal(tail(out$tV0, 1), 0)
  testthat::expect_equal(tail(out$tV1, 1), 0)
})


testthat::test_that("markov_nstep_prob handles zero steps", {
  P <- matrix(
    c(
      0.9, 0.1,
      0.0, 1.0
    ),
    nrow = 2,
    byrow = TRUE
  )

  testthat::expect_equal(
    markov_nstep_prob(P, n = 0, i = 1, j = 1),
    1
  )

  testthat::expect_equal(
    markov_nstep_prob(P, n = 0, i = 1, j = 2),
    0
  )
})


testthat::test_that("markov_nstep_prob validates transition matrices", {
  testthat::expect_error(
    markov_nstep_prob(
      matrix(c(0.8, 0.1, 0, 1), nrow = 2, byrow = TRUE),
      n = 1,
      i = 1,
      j = 1
    ),
    "row of P must sum to 1"
  )
})


testthat::test_that("gain_loss_md vectorizes", {
  out <- gain_loss_md(
    Vt = c(100, 110),
    G = 20,
    r = c(0.02, 0.03),
    e = 1,
    i = c(0.04, 0.05),
    b1 = 1000,
    b2 = 100,
    q1 = c(0.01, 0.012),
    q2 = c(0.08, 0.07),
    Vt1 = c(115, 125)
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("gain_loss_md supports the ordered decrement case", {
  out <- gain_loss_md(
    Vt = 100,
    G = 20,
    r = 0.02,
    e = 1,
    i = 0.04,
    b1 = 1000,
    b2 = 100,
    q1 = 0.01,
    q2 = 0.08,
    Vt1 = 115,
    year_end_cause2 = TRUE,
    q1prime = 0.01,
    q2prime = 0.08
  )

  expected <- (100 + 20 * 0.98 - 1) * 1.04 -
    (
      1000 * 0.01 +
        100 * (1 - 0.01) * 0.08 +
        (1 - 0.01) * (1 - 0.08) * 115
    )

  testthat::expect_equal(out, expected)
})
