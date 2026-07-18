testthat::test_that("salary scale follows compound growth", {
  out <- salary_scale(
    k = 30:34,
    g = 0.04,
    base_age = 30
  )

  testthat::expect_equal(
    out,
    1.04^(0:4),
    tolerance = 1e-12
  )
})


testthat::test_that("salary scale allows valid negative growth", {
  out <- salary_scale(
    k = 30:32,
    g = -0.02,
    base_age = 30
  )

  testthat::expect_equal(
    out,
    0.98^(0:2),
    tolerance = 1e-12
  )

  testthat::expect_error(
    salary_scale(30:32, g = -1),
    "greater than -1"
  )
})


testthat::test_that("defined contribution accumulation uses beginning-year contributions", {
  out <- AVz_dc(
    x = 30,
    z = 33,
    Sx = 100,
    c = 0.10,
    i = 0.05,
    g = 0
  )

  expected <- 10 * (
    1.05^3 +
      1.05^2 +
      1.05
  )

  testthat::expect_equal(
    out,
    expected,
    tolerance = 1e-12
  )
})


testthat::test_that("salary scale and growth inputs give equivalent DC values", {
  scale <- salary_scale(
    k = 30:34,
    g = 0.04,
    base_age = 30
  )

  from_growth <- AVz_dc(
    x = 30,
    z = 35,
    Sx = 50000,
    c = 0.10,
    i = 0.05,
    g = 0.04
  )

  from_scale <- AVz_dc(
    x = 30,
    z = 35,
    Sx = 50000,
    c = 0.10,
    i = 0.05,
    s = scale
  )

  testthat::expect_equal(
    from_growth,
    from_scale,
    tolerance = 1e-10
  )
})


testthat::test_that("defined contribution replacement ratio is internally consistent", {
  AV <- AVz_dc(
    x = 30,
    z = 35,
    Sx = 50000,
    c = 0.10,
    i = 0.05,
    g = 0.04
  )

  income <- Income_dc(AV, adue_z = 12)
  final_salary <- 50000 * 1.04^4

  expected <- income / final_salary

  testthat::expect_equal(
    replacement_ratio_dc(
      x = 30,
      z = 35,
      Sx = 50000,
      c = 0.10,
      i = 0.05,
      adue_z = 12,
      g = 0.04
    ),
    expected,
    tolerance = 1e-12
  )
})


testthat::test_that("target contribution rate reproduces the target ratio", {
  contribution <- contribution_rate_target(
    x = 30,
    z = 65,
    Sx = 60000,
    RR_target = 0.50,
    i = 0.06,
    adue_z = 11,
    g = 0.04
  )

  achieved <- replacement_ratio_dc(
    x = 30,
    z = 65,
    Sx = 60000,
    c = contribution,
    i = 0.06,
    adue_z = 11,
    g = 0.04
  )

  testthat::expect_equal(
    achieved,
    0.50,
    tolerance = 1e-10
  )
})


testthat::test_that("income and DB replacement ratios vectorize", {
  testthat::expect_equal(
    Income_dc(
      AVz = c(120000, 240000),
      adue_z = 12
    ),
    c(10000, 20000)
  )

  testthat::expect_equal(
    replacement_ratio_db(
      benefit = c(30000, 40000),
      salary = 100000
    ),
    c(0.30, 0.40)
  )
})


testthat::test_that("final average salary benefit uses the final salary years", {
  out <- PAB_fas(
    x = 30,
    z = 35,
    CASx = 100000,
    p = 2,
    fas_years = 3,
    g = 0.04
  )

  salaries <- 100000 * 1.04^(0:4)
  expected_fas <- mean(tail(salaries, 3))
  expected <- 0.02 * 5 * expected_fas

  testthat::expect_equal(
    out,
    expected,
    tolerance = 1e-10
  )
})


testthat::test_that("final average salary benefit includes past service", {
  without_past <- PAB_fas(
    x = 40,
    z = 45,
    CASx = 100000,
    p = 2,
    fas_years = 3,
    past_service = 0,
    g = 0.03
  )

  with_past <- PAB_fas(
    x = 40,
    z = 45,
    CASx = 100000,
    p = 2,
    fas_years = 3,
    past_service = 10,
    g = 0.03
  )

  testthat::expect_true(with_past > without_past)
})


testthat::test_that("career average projected benefit uses salary total", {
  salaries <- 100000 * 1.04^(0:4)

  out <- PAB_cae(
    x = 30,
    z = 35,
    CASx = 100000,
    p = 1,
    past_salary_total = 200000,
    g = 0.04
  )

  expected <- 0.01 * (
    200000 + sum(salaries)
  )

  testthat::expect_equal(
    out,
    expected,
    tolerance = 1e-10
  )
})


testthat::test_that("accrued benefit formulas use actual salary history", {
  salary_history <- c(
    100000,
    104000,
    108160
  )

  testthat::expect_equal(
    AB_cae(
      salary_history = salary_history,
      p = 1
    ),
    0.01 * sum(salary_history)
  )

  testthat::expect_equal(
    AB_fas(
      salary_history = salary_history,
      p = 2,
      fas_years = 2
    ),
    0.02 * 3 * mean(tail(salary_history, 2))
  )
})


testthat::test_that("normal retirement APV vectorizes consistently", {
  out <- APV_NR_db(
    PABz = c(10000, 20000),
    v_to_ret = 0.50,
    p_surv = c(0.80, 0.90),
    adue_ret = 12
  )

  expected <- c(10000, 20000) *
    0.50 *
    c(0.80, 0.90) *
    12

  testthat::expect_equal(
    out,
    expected
  )
})


testthat::test_that("TUC formulas agree with the normal retirement APV", {
  normal_cost <- NC_TUC_db(
    accrual_benefit = 1500,
    v_to_ret = 0.50,
    p_surv = 0.90,
    adue_ret = 12
  )

  liability <- AAL_TUC_db(
    accrued_benefit = 12000,
    v_to_ret = 0.50,
    p_surv = 0.90,
    adue_ret = 12
  )

  testthat::expect_equal(
    normal_cost,
    APV_NR_db(
      PABz = 1500,
      v_to_ret = 0.50,
      p_surv = 0.90,
      adue_ret = 12
    )
  )

  testthat::expect_equal(
    liability,
    APV_NR_db(
      PABz = 12000,
      v_to_ret = 0.50,
      p_surv = 0.90,
      adue_ret = 12
    )
  )
})


testthat::test_that("PUC cost and liability allocate projected benefit by service", {
  normal_cost <- NC_PUC_db(
    projected_benefit = 30000,
    total_service = 30,
    v_to_ret = 0.50,
    p_surv = 0.90,
    adue_ret = 12
  )

  liability <- AAL_PUC_db(
    projected_benefit = 30000,
    past_service = 10,
    total_service = 30,
    v_to_ret = 0.50,
    p_surv = 0.90,
    adue_ret = 12
  )

  testthat::expect_equal(
    normal_cost,
    APV_NR_db(
      PABz = 1000,
      v_to_ret = 0.50,
      p_surv = 0.90,
      adue_ret = 12
    )
  )

  testthat::expect_equal(
    liability,
    APV_NR_db(
      PABz = 10000,
      v_to_ret = 0.50,
      p_surv = 0.90,
      adue_ret = 12
    )
  )
})


testthat::test_that("PUC functions vectorize over members", {
  out <- AAL_PUC_db(
    projected_benefit = c(30000, 40000),
    past_service = c(10, 15),
    total_service = c(30, 40),
    v_to_ret = 0.50,
    p_surv = 0.90,
    adue_ret = 12
  )

  testthat::expect_length(out, 2)
  testthat::expect_true(all(is.finite(out)))
})


testthat::test_that("salary projections require integer projection years", {
  testthat::expect_error(
    AVz_dc(
      x = 30,
      z = 65.5,
      Sx = 50000,
      c = 0.10,
      i = 0.05,
      g = 0.04
    ),
    "integer number of years"
  )
})


testthat::test_that("salary projection requires exactly one salary assumption", {
  testthat::expect_error(
    AVz_dc(
      x = 30,
      z = 35,
      Sx = 50000,
      c = 0.10,
      i = 0.05
    ),
    "Supply one of g or s"
  )

  testthat::expect_error(
    AVz_dc(
      x = 30,
      z = 35,
      Sx = 50000,
      c = 0.10,
      i = 0.05,
      g = 0.04,
      s = rep(1, 5)
    ),
    "Supply only one"
  )
})


testthat::test_that("salary scale vector must match projection period", {
  testthat::expect_error(
    AVz_dc(
      x = 30,
      z = 35,
      Sx = 50000,
      c = 0.10,
      i = 0.05,
      s = rep(1, 4)
    ),
    "length z - x"
  )
})


testthat::test_that("pension functions reject invalid rates and percentages", {
  testthat::expect_error(
    AVz_dc(
      x = 30,
      z = 35,
      Sx = 50000,
      c = 1.20,
      i = 0.05,
      g = 0.04
    ),
    "between 0 and 1"
  )

  testthat::expect_error(
    AVz_dc(
      x = 30,
      z = 35,
      Sx = 50000,
      c = 0.10,
      i = -1,
      g = 0.04
    ),
    "greater than -1"
  )

  testthat::expect_error(
    PAB_fas(
      x = 30,
      z = 35,
      CASx = 50000,
      p = -1,
      g = 0.04
    ),
    "nonnegative"
  )
})


testthat::test_that("final average period must be a positive integer", {
  testthat::expect_error(
    PAB_fas(
      x = 30,
      z = 35,
      CASx = 50000,
      p = 2,
      fas_years = 2.5,
      g = 0.04
    ),
    "must be an integer"
  )

  testthat::expect_error(
    PAB_fas(
      x = 30,
      z = 35,
      CASx = 50000,
      p = 2,
      fas_years = 6,
      g = 0.04
    ),
    "cannot exceed"
  )
})


testthat::test_that("past service cannot exceed total service", {
  testthat::expect_error(
    AAL_PUC_db(
      projected_benefit = 30000,
      past_service = 31,
      total_service = 30,
      v_to_ret = 0.50,
      p_surv = 0.90,
      adue_ret = 12
    ),
    "cannot exceed"
  )
})


testthat::test_that("elementwise pension functions reject incompatible lengths", {
  testthat::expect_error(
    Income_dc(
      AVz = c(100000, 200000),
      adue_z = c(10, 11, 12)
    ),
    "length 1 or the common length"
  )

  testthat::expect_error(
    APV_NR_db(
      PABz = c(10000, 20000),
      v_to_ret = c(0.5, 0.6, 0.7),
      p_surv = 0.9,
      adue_ret = 12
    ),
    "length 1 or the common length"
  )
})
