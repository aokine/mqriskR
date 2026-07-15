test_that("qxtau and pxtau work", {
  expect_equal(qxtau(c(0.011, 0.100)), 0.111)
  expect_equal(pxtau(c(0.011, 0.100)), 0.889)
})

test_that("dxj and dxtau work", {
  expect_equal(dxj(1000, c(0.011, 0.100)), c(11, 100))
  expect_equal(dxtau(1000, c(0.011, 0.100)), 111)
})

test_that("md_table reproduces Example 13.1", {
  x <- 45:50
  qmat <- cbind(
    q1 = c(.011, .012, .013, .014, .015, .016),
    q2 = c(.100, .100, .100, .100, .100, .100)
  )

  tbl <- md_table(x = x, qxj = qmat, radix = 1000)

  expect_equal(round(tbl$qtau[1], 3), 0.111)
  expect_equal(round(tbl$ptau[1], 3), 0.889)
  expect_equal(round(tbl$ltau[2], 2), 889.00)
  expect_equal(round(tbl$d1[1], 2), 11.00)
  expect_equal(round(tbl$d2[1], 2), 100.00)
})

test_that("n-year probabilities from md_table match Example 13.2", {
  x <- 45:50
  qmat <- cbind(
    q1 = c(.011, .012, .013, .014, .015, .016),
    q2 = c(.100, .100, .100, .100, .100, .100)
  )

  tbl <- md_table(x = x, qxj = qmat, radix = 1000)

  expect_equal(round(npxtau_md(tbl, x = 46, n = 3), 4), 0.6979)
  expect_equal(nqxj_md(tbl, x = 47, n = 2, j = 2), (78.94 + 70.02) / 789.43, tolerance = 1e-4)
  expect_equal(round(nqxj_md(tbl, x = 46, n = 2, j = 1), 5), 0.02354)
})

test_that("constant-force formulas match Example 13.4", {
  mu <- c(0.10, 0.20)
  t <- 5

  expect_equal(round(tpx_tau_cf(mu, t), 6), round(exp(-0.30 * t), 6))
  expect_equal(round(tqxprimej_cf(0.10, t), 6), round(1 - exp(-0.10 * t), 6))
  expect_equal(round(tqxj_cf(mu, 1, t), 6), round((1/3) * (1 - exp(-0.30 * t)), 6))
  expect_equal(round(tqxj_cf(mu, 1, 100), 6), round(1/3, 6))
})

test_that("MUDD matches Example 13.5", {
  out <- qxprime_mudd(c(.20, .10))
  expect_equal(round(out[2], 5), 0.11210)
})

test_that("SUDD matches Example 13.6", {
  out <- qxprime_sudd(q1 = .20, q2 = .10)
  expect_equal(unname(out["q2prime"]), 0.11184, tolerance = 1e-4)
})

test_that("constant-force and MUDD inverse from dependent q agree", {
  qdep <- c(.20, .10)
  out_mudd <- qxprime_mudd(qdep)
  expect_equal(length(out_mudd), 2)

  qprime <- c(.20, .10)
  qdep_cf <- qx_dep_cf(qprime)
  expect_true(all(qdep_cf >= 0))
  expect_true(sum(qdep_cf) <= 1)
})


testthat::test_that("multiple-decrement table validates age sequences", {
  qmat <- cbind(
    q1 = c(0.01, 0.02, 0.03),
    q2 = c(0.05, 0.05, 0.05)
  )

  testthat::expect_error(
    md_table(c(40, 42, 43), qmat),
    "consecutive"
  )

  testthat::expect_error(
    md_table(c(40, 40, 41), qmat),
    "duplicate"
  )

  testthat::expect_error(
    md_table(c(40, 41.5, 42), qmat),
    "integer-like"
  )
})


testthat::test_that("multiple-decrement table retains expected structure", {
  tbl <- md_table(
    x = 40:42,
    qxj = cbind(
      withdrawal = c(0.01, 0.02, 0.03),
      retirement = c(0.10, 0.10, 0.10)
    ),
    radix = 1000
  )

  testthat::expect_s3_class(tbl, "md_table")
  testthat::expect_s3_class(tbl, "data.frame")

  testthat::expect_true(all(
    c(
      "x", "q1", "q2", "qtau", "ptau",
      "ltau", "d1", "d2", "dtau"
    ) %in% names(tbl)
  ))

  testthat::expect_equal(
    attr(tbl, "cause_names"),
    c("withdrawal", "retirement")
  )
})


testthat::test_that("table probabilities satisfy total-decrement identities", {
  tbl <- md_table(
    x = 40:45,
    qxj = cbind(
      q1 = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06),
      q2 = rep(0.10, 6)
    ),
    radix = 1000
  )

  total <- nqxtau_md(tbl, x = 40, n = 4)
  cause_1 <- nqxj_md(tbl, x = 40, n = 4, j = 1)
  cause_2 <- nqxj_md(tbl, x = 40, n = 4, j = 2)

  testthat::expect_equal(
    total,
    cause_1 + cause_2,
    tolerance = 1e-12
  )

  testthat::expect_equal(
    total,
    1 - npxtau_md(tbl, x = 40, n = 4),
    tolerance = 1e-12
  )
})


testthat::test_that("table functions handle zero terms", {
  tbl <- md_table(
    x = 40:42,
    qxj = cbind(
      q1 = c(0.01, 0.02, 0.03),
      q2 = c(0.10, 0.10, 0.10)
    )
  )

  testthat::expect_equal(
    npxtau_md(tbl, x = 42, n = 0),
    1
  )

  testthat::expect_equal(
    nqxtau_md(tbl, x = 42, n = 0),
    0
  )

  testthat::expect_equal(
    nqxj_md(tbl, x = 42, n = 0, j = 1),
    0
  )
})


testthat::test_that("table functions vectorize over x, n, and cause", {
  tbl <- md_table(
    x = 40:45,
    qxj = cbind(
      q1 = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06),
      q2 = rep(0.10, 6)
    )
  )

  survival <- npxtau_md(
    tbl,
    x = c(40, 42),
    n = c(2, 3)
  )

  causes <- nqxj_md(
    tbl,
    x = c(40, 42),
    n = c(2, 3),
    j = c(1, 2)
  )

  testthat::expect_length(survival, 2)
  testthat::expect_length(causes, 2)
  testthat::expect_true(all(survival >= 0 & survival <= 1))
  testthat::expect_true(all(causes >= 0 & causes <= 1))
})


testthat::test_that("table functions reject non-integer terms", {
  tbl <- md_table(
    x = 40:42,
    qxj = cbind(
      q1 = c(0.01, 0.02, 0.03),
      q2 = c(0.10, 0.10, 0.10)
    )
  )

  testthat::expect_error(
    npxtau_md(tbl, x = 40, n = 1.5),
    "integer-like"
  )

  testthat::expect_error(
    nqxj_md(tbl, x = 40, n = 1.5, j = 1),
    "integer-like"
  )
})


testthat::test_that("constant-force functions satisfy probability identities", {
  mu <- c(0.10, 0.20)
  t <- c(1, 5, 10)

  total_failure <- 1 - tpx_tau_cf(mu, t)

  cause_sum <- tqxj_cf(mu, j = 1, t = t) +
    tqxj_cf(mu, j = 2, t = t)

  testthat::expect_equal(
    cause_sum,
    total_failure,
    tolerance = 1e-12
  )
})


testthat::test_that("constant-force single-decrement function vectorizes", {
  out <- tpxprimej_cf(
    mu = c(0.10, 0.20),
    t = 5
  )

  testthat::expect_length(out, 2)

  testthat::expect_equal(
    out,
    exp(-c(0.10, 0.20) * 5)
  )
})


testthat::test_that("MUDD handles the zero-decrement boundary", {
  testthat::expect_equal(
    qxprime_mudd(c(0, 0)),
    c(0, 0)
  )

  testthat::expect_equal(
    tqxprime_mudd(c(0, 0), t = 0.5),
    c(0, 0)
  )
})


testthat::test_that("constant-force conversion rejects probability one", {
  testthat::expect_error(
    qx_dep_cf(c(1, 0.10)),
    "less than 1"
  )
})


testthat::test_that("SUDD forward and inverse conversions agree", {
  independent <- c(q1prime = 0.20, q2prime = 0.10)

  dependent <- qx_dep_sudd(
    q1prime = independent["q1prime"],
    q2prime = independent["q2prime"]
  )

  recovered <- qxprime_sudd(
    q1 = dependent["q1"],
    q2 = dependent["q2"]
  )

  testthat::expect_equal(
    unname(recovered),
    unname(independent),
    tolerance = 1e-10
  )
})


testthat::test_that("SUDD functions vectorize pairwise", {
  dependent <- qx_dep_sudd(
    q1prime = c(0.20, 0.30),
    q2prime = c(0.10, 0.15)
  )

  testthat::expect_true(is.matrix(dependent))
  testthat::expect_equal(dim(dependent), c(2, 2))
  testthat::expect_equal(
    colnames(dependent),
    c("q1", "q2")
  )

  recovered <- qxprime_sudd(
    q1 = dependent[, "q1"],
    q2 = dependent[, "q2"]
  )

  testthat::expect_true(is.matrix(recovered))

  testthat::expect_equal(
    recovered[, "q1prime"],
    c(0.20, 0.30),
    tolerance = 1e-10
  )

  testthat::expect_equal(
    recovered[, "q2prime"],
    c(0.10, 0.15),
    tolerance = 1e-10
  )
})
