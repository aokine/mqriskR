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
