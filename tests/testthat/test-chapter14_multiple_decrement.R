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
