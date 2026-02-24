testthat::test_that("Uniform (de Moivre) model basics", {
  omega <- 100

  # endpoints / basic values
  testthat::expect_equal(S0(0, "uniform", omega = omega), 1)
  testthat::expect_equal(F0(0, "uniform", omega = omega), 0)
  testthat::expect_equal(S0(omega, "uniform", omega = omega), 0)
  testthat::expect_equal(F0(omega, "uniform", omega = omega), 1)

  # formula S0(t) = (omega - t)/omega on [0, omega]
  testthat::expect_equal(S0(25, "uniform", omega = omega), (omega - 25) / omega)

  # density constant on (0, omega), 0 outside
  testthat::expect_equal(f0(1, "uniform", omega = omega), 1 / omega)
  testthat::expect_equal(f0(0, "uniform", omega = omega), 0)
  testthat::expect_equal(f0(omega, "uniform", omega = omega), 0)
  testthat::expect_equal(f0(omega + 1, "uniform", omega = omega), 0)

  # hazard is 1/(omega - t) for t < omega; Inf at/after omega
  testthat::expect_equal(hazard0(10, "uniform", omega = omega), 1 / (omega - 10))
  testthat::expect_true(is.infinite(hazard0(omega, "uniform", omega = omega)))

  # cumulative hazard Lambda(t) = log(omega/(omega - t)) for t < omega
  testthat::expect_equal(
    cumhaz0(10, "uniform", omega = omega),
    log(omega / (omega - 10)),
    tolerance = 1e-12
  )
  testthat::expect_true(is.infinite(cumhaz0(omega, "uniform", omega = omega)))
})

testthat::test_that("Exponential model basics", {
  lam <- 0.02
  t <- 10

  testthat::expect_equal(S0(t, "exponential", lambda = lam), exp(-lam * t), tolerance = 1e-12)
  testthat::expect_equal(F0(t, "exponential", lambda = lam), 1 - exp(-lam * t), tolerance = 1e-12)
  testthat::expect_equal(f0(t, "exponential", lambda = lam), lam * exp(-lam * t), tolerance = 1e-12)

  # constant hazard and cumhaz
  testthat::expect_equal(hazard0(t, "exponential", lambda = lam), lam, tolerance = 1e-12)
  testthat::expect_equal(cumhaz0(t, "exponential", lambda = lam), lam * t, tolerance = 1e-12)

  # identity: f0 = hazard0 * S0
  testthat::expect_equal(
    f0(t, "exponential", lambda = lam),
    hazard0(t, "exponential", lambda = lam) * S0(t, "exponential", lambda = lam),
    tolerance = 1e-12
  )
})

testthat::test_that("Gompertz model: hazard and cumhaz consistency", {
  B <- 0.0005
  c <- 1.08
  t <- 12

  # hazard is B*c^t
  testthat::expect_equal(hazard0(t, "gompertz", B = B, c = c), B * c^t, tolerance = 1e-12)

  # cumhaz is (B/log c)(c^t - 1)
  testthat::expect_equal(
    cumhaz0(t, "gompertz", B = B, c = c),
    (B / log(c)) * (c^t - 1),
    tolerance = 1e-12
  )

  # S0 should match exp(-cumhaz) for Gompertz
  testthat::expect_equal(
    S0(t, "gompertz", B = B, c = c),
    exp(-cumhaz0(t, "gompertz", B = B, c = c)),
    tolerance = 1e-12
  )
})

testthat::test_that("Makeham model: hazard and cumhaz consistency", {
  A <- 0.0002
  B <- 0.0004
  c <- 1.07
  t <- 15

  testthat::expect_equal(hazard0(t, "makeham", A = A, B = B, c = c), A + B * c^t, tolerance = 1e-12)
  testthat::expect_equal(
    cumhaz0(t, "makeham", A = A, B = B, c = c),
    A * t + (B / log(c)) * (c^t - 1),
    tolerance = 1e-12
  )

  testthat::expect_equal(
    S0(t, "makeham", A = A, B = B, c = c),
    exp(-cumhaz0(t, "makeham", A = A, B = B, c = c)),
    tolerance = 1e-12
  )
})

testthat::test_that("Weibull model: hazard/cumhaz consistency", {
  shape <- 2.3
  scale <- 50
  t <- 10

  testthat::expect_equal(
    cumhaz0(t, "weibull", shape = shape, scale = scale),
    (t / scale)^shape,
    tolerance = 1e-12
  )

  testthat::expect_equal(
    hazard0(t, "weibull", shape = shape, scale = scale),
    (shape / scale) * (t / scale)^(shape - 1),
    tolerance = 1e-12
  )

  testthat::expect_equal(
    S0(t, "weibull", shape = shape, scale = scale),
    exp(-cumhaz0(t, "weibull", shape = shape, scale = scale)),
    tolerance = 1e-12
  )
})

testthat::test_that("Conditional relationships: tpx, tqx, fx", {
  omega <- 120
  x <- 30
  t <- 10

  # tpx = S0(x+t)/S0(x)
  testthat::expect_equal(
    tpx(t, x, "uniform", omega = omega),
    S0(x + t, "uniform", omega = omega) / S0(x, "uniform", omega = omega),
    tolerance = 1e-12
  )

  # tqx = 1 - tpx
  testthat::expect_equal(
    tqx(t, x, "uniform", omega = omega),
    1 - tpx(t, x, "uniform", omega = omega),
    tolerance = 1e-12
  )

  # fx(t) = f0(x+t)/S0(x)
  testthat::expect_equal(
    fx(t, x, "uniform", omega = omega),
    f0(x + t, "uniform", omega = omega) / S0(x, "uniform", omega = omega),
    tolerance = 1e-12
  )
})

testthat::test_that("Vectorization over t (and x where relevant)", {
  omega <- 100
  t <- c(0, 10, 25)
  x <- c(20, 30, 40)

  testthat::expect_length(S0(t, "uniform", omega = omega), length(t))
  testthat::expect_length(F0(t, "uniform", omega = omega), length(t))
  testthat::expect_length(f0(t, "uniform", omega = omega), length(t))
  testthat::expect_length(hazard0(t, "uniform", omega = omega), length(t))
  testthat::expect_length(cumhaz0(t, "uniform", omega = omega), length(t))

  testthat::expect_length(tpx(t, x, "uniform", omega = omega), length(t))
  testthat::expect_length(tqx(t, x, "uniform", omega = omega), length(t))
  testthat::expect_length(fx(t, x, "uniform", omega = omega), length(t))
})

testthat::test_that("Monotonicity / sanity checks", {
  lam <- 0.03
  t <- c(0, 5, 10, 20)

  # S0 decreasing; F0 increasing
  S <- S0(t, "exponential", lambda = lam)
  F <- F0(t, "exponential", lambda = lam)
  testthat::expect_true(all(diff(S) <= 0))
  testthat::expect_true(all(diff(F) >= 0))

  # hazard nonnegative
  testthat::expect_true(all(hazard0(t, "exponential", lambda = lam) >= 0))
})

testthat::test_that("ex_complete closed forms: uniform and exponential", {
  # uniform: e_x^o = (omega - x)/2 for x < omega; 0 for x >= omega
  omega <- 100
  testthat::expect_equal(ex_complete(20, "uniform", omega = omega), (omega - 20) / 2, tolerance = 1e-10)
  testthat::expect_equal(ex_complete(100, "uniform", omega = omega), 0)

  # exponential: e_x^o = 1/lambda (memoryless)
  lam <- 0.04
  testthat::expect_equal(ex_complete(0, "exponential", lambda = lam), 1 / lam, tolerance = 1e-10)
  testthat::expect_equal(ex_complete(50, "exponential", lambda = lam), 1 / lam, tolerance = 1e-10)
})

testthat::test_that("ex_curtate basic properties (uniform and exponential)", {
  # uniform: curtate should be finite and <= (omega-x) roughly
  omega <- 120
  x <- 30
  ec <- ex_curtate(x, "uniform", omega = omega, k_max = 5000)
  testthat::expect_true(is.finite(ec))
  testthat::expect_true(ec >= 0)
  testthat::expect_true(ec <= (omega - x) + 1) # loose sanity bound

  # exponential: curtate expectation should be close to sum_{k>=1} exp(-lambda k)
  lam <- 0.05
  approx_geom <- sum(exp(-lam * (1:2000)))  # high truncation
  ec2 <- ex_curtate(10, "exponential", lambda = lam, k_max = 5000)
  testthat::expect_equal(ec2, approx_geom, tolerance = 1e-6)
})

testthat::test_that("Input validation errors", {
  # negative t or x
  testthat::expect_error(S0(-1, "exponential", lambda = 0.1))
  testthat::expect_error(tpx(1, -2, "exponential", lambda = 0.1))
  testthat::expect_error(tpx(-1, 2, "exponential", lambda = 0.1))

  # unknown model
  testthat::expect_error(S0(1, "notamodel", lambda = 0.1))

  # missing parameters
  testthat::expect_error(S0(1, "uniform"))
  testthat::expect_error(S0(1, "exponential"))
  testthat::expect_error(S0(1, "gompertz", B = 0.1))        # missing c
  testthat::expect_error(S0(1, "makeham", A = 0.1, B = 0.1))# missing c

  # invalid parameters
  testthat::expect_error(S0(1, "uniform", omega = -1))
  testthat::expect_error(S0(1, "exponential", lambda = 0))
  testthat::expect_error(S0(1, "gompertz", B = 0, c = 1.1))
  testthat::expect_error(S0(1, "gompertz", B = 0.1, c = 1))
  testthat::expect_error(S0(1, "weibull", shape = 0, scale = 10))
  testthat::expect_error(S0(1, "weibull", shape = 2, scale = 0))
})
