# =========================================================
# Pension mathematics functions
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.pension_check_numeric <- function(x, name, allow_empty = FALSE) {
  if (!is.numeric(x)) {
    stop(name, " must be numeric.", call. = FALSE)
  }

  if (!allow_empty && length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop(name, " must contain finite values.", call. = FALSE)
  }

  as.numeric(x)
}

#' @noRd
.pension_check_nonnegative <- function(x, name) {
  x <- .pension_check_numeric(x, name)

  if (any(x < 0)) {
    stop(name, " must contain nonnegative values.", call. = FALSE)
  }

  x
}

#' @noRd
.pension_check_positive <- function(x, name) {
  x <- .pension_check_numeric(x, name)

  if (any(x <= 0)) {
    stop(name, " must contain positive values.", call. = FALSE)
  }

  x
}

#' @noRd
.pension_check_probability <- function(x, name) {
  x <- .pension_check_numeric(x, name)

  if (any(x < 0 | x > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }

  x
}

#' @noRd
.pension_check_rate <- function(x, name) {
  x <- .pension_check_numeric(x, name)

  if (any(x <= -1)) {
    stop(name, " must contain values greater than -1.", call. = FALSE)
  }

  x
}

#' @noRd
.pension_check_scalar <- function(
    x,
    name,
    nonnegative = FALSE,
    positive = FALSE,
    integer = FALSE) {
  if (positive) {
    x <- .pension_check_positive(x, name)
  } else if (nonnegative) {
    x <- .pension_check_nonnegative(x, name)
  } else {
    x <- .pension_check_numeric(x, name)
  }

  if (length(x) != 1L) {
    stop(name, " must be a numeric scalar.", call. = FALSE)
  }

  if (integer && abs(x - round(x)) > 1e-10) {
    stop(name, " must be an integer.", call. = FALSE)
  }

  if (integer) {
    x <- as.integer(round(x))
  }

  x
}

#' @noRd
.pension_recycle_common <- function(..., .names = NULL) {
  values <- list(...)
  lengths <- vapply(values, length, integer(1))

  if (any(lengths == 0L)) {
    stop("Inputs must have positive length.", call. = FALSE)
  }

  common_length <- max(lengths)

  if (is.null(.names)) {
    .names <- paste0("Argument ", seq_along(values))
  }

  if (length(.names) != length(values)) {
    stop(
      ".names must have the same length as the supplied arguments.",
      call. = FALSE
    )
  }

  for (j in seq_along(values)) {
    if (!lengths[j] %in% c(1L, common_length)) {
      stop(
        .names[j],
        " must have length 1 or the common length ",
        common_length,
        ".",
        call. = FALSE
      )
    }

    values[[j]] <- rep(values[[j]], length.out = common_length)
  }

  values
}

#' @noRd
.pension_salary_path <- function(x, z, Sx, g = NULL, s = NULL) {
  x <- .pension_check_scalar(x, "x")
  z <- .pension_check_scalar(z, "z")
  Sx <- .pension_check_scalar(Sx, "Sx", positive = TRUE)

  if (z <= x) {
    stop("z must be greater than x.", call. = FALSE)
  }

  n_raw <- z - x

  if (abs(n_raw - round(n_raw)) > 1e-10) {
    stop("z - x must be an integer number of years.", call. = FALSE)
  }

  n <- as.integer(round(n_raw))

  if (!is.null(g) && !is.null(s)) {
    stop("Supply only one of g or s.", call. = FALSE)
  }

  if (is.null(g) && is.null(s)) {
    stop("Supply one of g or s.", call. = FALSE)
  }

  if (!is.null(g)) {
    g <- .pension_check_scalar(g, "g")

    if (g <= -1) {
      stop("g must be greater than -1.", call. = FALSE)
    }

    return(Sx * (1 + g)^(0:(n - 1L)))
  }

  s <- .pension_check_positive(s, "s")

  if (length(s) != n) {
    stop("s must have length z - x.", call. = FALSE)
  }

  Sx * s / s[1L]
}

# -------------------------------------------------------------------------
# Salary scales and defined contribution plans
# -------------------------------------------------------------------------

#' Salary scale under constant annual growth
#'
#' Constructs salary-scale values under a constant annual growth rate.
#'
#' @param k Numeric vector of ages or durations.
#' @param g Annual salary growth rate greater than \code{-1}.
#' @param base_age Scalar age or duration at which the scale is normalized.
#' @param s_base Positive scalar salary-scale value at \code{base_age}.
#'
#' @return A numeric vector with the same length as \code{k}.
#'
#' @examples
#' salary_scale(k = 30:34, g = 0.04, base_age = 30)
#'
#' @export
salary_scale <- function(k, g, base_age = min(k), s_base = 1) {
  k <- .pension_check_numeric(k, "k")
  g <- .pension_check_scalar(g, "g")
  base_age <- .pension_check_scalar(base_age, "base_age")
  s_base <- .pension_check_scalar(
    s_base,
    "s_base",
    positive = TRUE
  )

  if (g <= -1) {
    stop("g must be greater than -1.", call. = FALSE)
  }

  s_base * (1 + g)^(k - base_age)
}

#' Accumulated value of defined contribution plan contributions
#'
#' Computes the accumulated value at retirement from contributions paid at
#' the beginning of each year and accumulated to retirement.
#'
#' @param x Scalar entry age.
#' @param z Scalar retirement age.
#' @param Sx Positive scalar salary at age \code{x}.
#' @param c Scalar contribution rate in \code{[0, 1]}.
#' @param i Scalar annual effective investment return greater than \code{-1}.
#' @param g Optional scalar annual salary growth rate greater than \code{-1}.
#' @param s Optional positive salary-scale vector of length \code{z - x}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' AVz_dc(
#'   x = 30,
#'   z = 65,
#'   Sx = 50000,
#'   c = 0.10,
#'   i = 0.05,
#'   g = 0.04
#' )
#'
#' @export
AVz_dc <- function(x, z, Sx, c, i, g = NULL, s = NULL) {
  c <- .pension_check_scalar(c, "c")

  if (c < 0 || c > 1) {
    stop("c must be between 0 and 1.", call. = FALSE)
  }

  i <- .pension_check_scalar(i, "i")

  if (i <= -1) {
    stop("i must be greater than -1.", call. = FALSE)
  }

  salary <- .pension_salary_path(
    x = x,
    z = z,
    Sx = Sx,
    g = g,
    s = s
  )

  n <- length(salary)
  times <- 0:(n - 1L)

  sum(c * salary * (1 + i)^(n - times))
}

#' Retirement income from a defined contribution accumulation
#'
#' Converts an accumulated account value into annual annuity-due income.
#' Scalar arguments are recycled to a common length.
#'
#' @param AVz Nonnegative accumulated value at retirement.
#' @param adue_z Positive whole-life annuity-due factor at retirement.
#'
#' @return A numeric vector.
#'
#' @examples
#' Income_dc(AVz = 824211.35, adue_z = 12)
#'
#' @export
Income_dc <- function(AVz, adue_z) {
  AVz <- .pension_check_nonnegative(AVz, "AVz")
  adue_z <- .pension_check_positive(adue_z, "adue_z")

  values <- .pension_recycle_common(
    AVz,
    adue_z,
    .names = c("AVz", "adue_z")
  )

  values[[1L]] / values[[2L]]
}

#' Replacement ratio for a defined contribution plan
#'
#' Computes annual retirement income divided by salary in the final
#' pre-retirement year.
#'
#' @inheritParams AVz_dc
#' @param adue_z Positive scalar whole-life annuity-due factor at retirement.
#'
#' @return A numeric scalar.
#'
#' @examples
#' replacement_ratio_dc(
#'   x = 30,
#'   z = 65,
#'   Sx = 50000,
#'   c = 0.10,
#'   i = 0.05,
#'   adue_z = 12,
#'   g = 0.04
#' )
#'
#' @export
replacement_ratio_dc <- function(
    x,
    z,
    Sx,
    c,
    i,
    adue_z,
    g = NULL,
    s = NULL) {
  adue_z <- .pension_check_scalar(
    adue_z,
    "adue_z",
    positive = TRUE
  )

  accumulated_value <- AVz_dc(
    x = x,
    z = z,
    Sx = Sx,
    c = c,
    i = i,
    g = g,
    s = s
  )

  salary <- .pension_salary_path(
    x = x,
    z = z,
    Sx = Sx,
    g = g,
    s = s
  )

  Income_dc(accumulated_value, adue_z) /
    tail(salary, 1L)
}

#' Target contribution rate for a defined contribution plan
#'
#' Calculates the contribution rate required to achieve a target replacement
#' ratio.
#'
#' @param x Scalar entry age.
#' @param z Scalar retirement age.
#' @param Sx Positive scalar salary at age \code{x}.
#' @param RR_target Scalar target replacement ratio in \code{[0, 1]}.
#' @param i Scalar annual effective investment return greater than \code{-1}.
#' @param adue_z Positive scalar whole-life annuity-due factor at retirement.
#' @param g Optional scalar annual salary growth rate greater than \code{-1}.
#' @param s Optional positive salary-scale vector of length \code{z - x}.
#'
#' @return A numeric scalar. The result may exceed one when the target cannot
#'   be achieved with a contribution rate no greater than 100 percent.
#'
#' @examples
#' contribution_rate_target(
#'   x = 30,
#'   z = 65,
#'   Sx = 60000,
#'   RR_target = 0.50,
#'   i = 0.06,
#'   adue_z = 11,
#'   g = 0.04
#' )
#'
#' @export
contribution_rate_target <- function(
    x,
    z,
    Sx,
    RR_target,
    i,
    adue_z,
    g = NULL,
    s = NULL) {
  RR_target <- .pension_check_scalar(
    RR_target,
    "RR_target"
  )

  if (RR_target < 0 || RR_target > 1) {
    stop(
      "RR_target must be between 0 and 1.",
      call. = FALSE
    )
  }

  i <- .pension_check_scalar(i, "i")

  if (i <= -1) {
    stop("i must be greater than -1.", call. = FALSE)
  }

  adue_z <- .pension_check_scalar(
    adue_z,
    "adue_z",
    positive = TRUE
  )

  salary <- .pension_salary_path(
    x = x,
    z = z,
    Sx = Sx,
    g = g,
    s = s
  )

  n <- length(salary)
  times <- 0:(n - 1L)

  numerator <-
    RR_target *
    tail(salary, 1L) *
    adue_z

  denominator <-
    sum(salary * (1 + i)^(n - times))

  numerator / denominator
}

# -------------------------------------------------------------------------
# Defined benefit formulas
# -------------------------------------------------------------------------

#' Projected annual benefit under a final-average-salary plan
#'
#' Projects the annual benefit using a final-average-salary formula.
#'
#' @param x Scalar current or entry age.
#' @param z Scalar retirement age.
#' @param CASx Positive scalar current annual salary.
#' @param p Nonnegative scalar accrual percentage, such as \code{2} for
#'   2 percent.
#' @param fas_years Positive integer number of years in the final salary
#'   average.
#' @param past_service Nonnegative scalar years of service already completed.
#' @param g Optional scalar annual salary growth rate greater than \code{-1}.
#' @param s Optional positive salary-scale vector of length \code{z - x}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' PAB_fas(
#'   x = 35,
#'   z = 65,
#'   CASx = 60000,
#'   p = 2,
#'   fas_years = 3,
#'   g = 0.04
#' )
#'
#' @export
PAB_fas <- function(
    x,
    z,
    CASx,
    p,
    fas_years = 3,
    past_service = 0,
    g = NULL,
    s = NULL) {
  p <- .pension_check_scalar(
    p,
    "p",
    nonnegative = TRUE
  )

  fas_years <- .pension_check_scalar(
    fas_years,
    "fas_years",
    positive = TRUE,
    integer = TRUE
  )

  past_service <- .pension_check_scalar(
    past_service,
    "past_service",
    nonnegative = TRUE
  )

  salary <- .pension_salary_path(
    x = x,
    z = z,
    Sx = CASx,
    g = g,
    s = s
  )

  n <- length(salary)

  if (fas_years > n) {
    stop(
      "fas_years cannot exceed z - x.",
      call. = FALSE
    )
  }

  final_average_salary <-
    mean(tail(salary, fas_years))

  total_service <- past_service + n

  0.01 * p * total_service * final_average_salary
}

#' Projected annual benefit under a career-average-earnings plan
#'
#' Projects the annual benefit using a career-average-earnings formula.
#'
#' @param x Scalar current or entry age.
#' @param z Scalar retirement age.
#' @param CASx Positive scalar current annual salary.
#' @param p Nonnegative scalar accrual percentage.
#' @param past_salary_total Nonnegative total of actual prior salaries.
#' @param g Optional scalar annual salary growth rate greater than \code{-1}.
#' @param s Optional positive salary-scale vector of length \code{z - x}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' PAB_cae(
#'   x = 30,
#'   z = 65,
#'   CASx = 100000,
#'   p = 1,
#'   g = 0.04
#' )
#'
#' @export
PAB_cae <- function(
    x,
    z,
    CASx,
    p,
    past_salary_total = 0,
    g = NULL,
    s = NULL) {
  p <- .pension_check_scalar(
    p,
    "p",
    nonnegative = TRUE
  )

  past_salary_total <- .pension_check_scalar(
    past_salary_total,
    "past_salary_total",
    nonnegative = TRUE
  )

  salary <- .pension_salary_path(
    x = x,
    z = z,
    Sx = CASx,
    g = g,
    s = s
  )

  projected_salary_total <-
    past_salary_total + sum(salary)

  0.01 * p * projected_salary_total
}

#' Replacement ratio for a defined benefit plan
#'
#' Computes annual benefit divided by a selected salary measure. Scalar
#' arguments are recycled to a common length.
#'
#' @param benefit Nonnegative annual retirement benefit.
#' @param salary Positive salary measure used in the denominator.
#'
#' @return A numeric vector.
#'
#' @examples
#' replacement_ratio_db(
#'   benefit = 108008.66,
#'   salary = 187119.09
#' )
#'
#' @export
replacement_ratio_db <- function(benefit, salary) {
  benefit <- .pension_check_nonnegative(
    benefit,
    "benefit"
  )

  salary <- .pension_check_positive(
    salary,
    "salary"
  )

  values <- .pension_recycle_common(
    benefit,
    salary,
    .names = c("benefit", "salary")
  )

  values[[1L]] / values[[2L]]
}

#' Accrued benefit under a final-average-salary plan
#'
#' Computes the accrued benefit using salary and service history through the
#' valuation date.
#'
#' @param salary_history Positive numeric vector of annual salaries.
#' @param p Nonnegative scalar accrual percentage.
#' @param fas_years Positive integer number of years in the salary average.
#'
#' @return A numeric scalar.
#'
#' @examples
#' AB_fas(
#'   salary_history = c(150000, 156000),
#'   p = 1,
#'   fas_years = 2
#' )
#'
#' @export
AB_fas <- function(salary_history, p, fas_years = 3) {
  salary_history <- .pension_check_positive(
    salary_history,
    "salary_history"
  )

  p <- .pension_check_scalar(
    p,
    "p",
    nonnegative = TRUE
  )

  fas_years <- .pension_check_scalar(
    fas_years,
    "fas_years",
    positive = TRUE,
    integer = TRUE
  )

  years <- length(salary_history)

  final_average_salary <- mean(
    tail(
      salary_history,
      min(fas_years, years)
    )
  )

  0.01 * p * years * final_average_salary
}

#' Accrued benefit under a career-average-earnings plan
#'
#' Computes the accrued benefit using salary history through the valuation
#' date.
#'
#' @param salary_history Positive numeric vector of annual salaries.
#' @param p Nonnegative scalar accrual percentage.
#'
#' @return A numeric scalar.
#'
#' @examples
#' AB_cae(
#'   salary_history = c(100000, 104000, 108160),
#'   p = 1
#' )
#'
#' @export
AB_cae <- function(salary_history, p) {
  salary_history <- .pension_check_positive(
    salary_history,
    "salary_history"
  )

  p <- .pension_check_scalar(
    p,
    "p",
    nonnegative = TRUE
  )

  0.01 * p * sum(salary_history)
}

# -------------------------------------------------------------------------
# Defined benefit valuation and cost methods
# -------------------------------------------------------------------------

#' Actuarial present value of a normal retirement benefit
#'
#' Computes the actuarial present value of a projected annual retirement
#' benefit. Scalar arguments are recycled to a common length.
#'
#' @param PABz Nonnegative projected annual benefit at retirement.
#' @param v_to_ret Nonnegative discount factor from the valuation date to
#'   retirement.
#' @param p_surv Survival or active-service probability to retirement in
#'   \code{[0, 1]}.
#' @param adue_ret Positive retirement annuity-due factor.
#'
#' @return A numeric vector.
#'
#' @examples
#' APV_NR_db(
#'   PABz = 108008.66,
#'   v_to_ret = 1 / 1.06^30,
#'   p_surv = 0.8,
#'   adue_ret = 12
#' )
#'
#' @export
APV_NR_db <- function(PABz, v_to_ret, p_surv, adue_ret) {
  PABz <- .pension_check_nonnegative(PABz, "PABz")
  v_to_ret <- .pension_check_nonnegative(
    v_to_ret,
    "v_to_ret"
  )
  p_surv <- .pension_check_probability(
    p_surv,
    "p_surv"
  )
  adue_ret <- .pension_check_positive(
    adue_ret,
    "adue_ret"
  )

  values <- .pension_recycle_common(
    PABz,
    v_to_ret,
    p_surv,
    adue_ret,
    .names = c(
      "PABz",
      "v_to_ret",
      "p_surv",
      "adue_ret"
    )
  )

  values[[1L]] *
    values[[2L]] *
    values[[3L]] *
    values[[4L]]
}

#' Entry Age Normal normal cost
#'
#' Computes Entry Age Normal normal cost as total benefit APV divided by an
#' active-service annuity-due factor.
#'
#' @param APV_total Nonnegative total actuarial present value of benefits.
#' @param adue_active Positive active-service annuity-due factor.
#'
#' @return A numeric vector.
#'
#' @examples
#' NC_EAN_db(APV_total = 25000, adue_active = 15)
#'
#' @export
NC_EAN_db <- function(APV_total, adue_active) {
  APV_total <- .pension_check_nonnegative(
    APV_total,
    "APV_total"
  )

  adue_active <- .pension_check_positive(
    adue_active,
    "adue_active"
  )

  values <- .pension_recycle_common(
    APV_total,
    adue_active,
    .names = c("APV_total", "adue_active")
  )

  values[[1L]] / values[[2L]]
}

#' Traditional Unit Credit normal cost
#'
#' Computes the actuarial present value of the benefit accrued during the
#' current year.
#'
#' @param accrual_benefit Nonnegative benefit accrued during the current year.
#' @inheritParams APV_NR_db
#'
#' @return A numeric vector.
#'
#' @examples
#' NC_TUC_db(
#'   accrual_benefit = 1560,
#'   v_to_ret = 0.5,
#'   p_surv = 0.9,
#'   adue_ret = 12
#' )
#'
#' @export
NC_TUC_db <- function(
    accrual_benefit,
    v_to_ret,
    p_surv,
    adue_ret) {
  APV_NR_db(
    PABz = accrual_benefit,
    v_to_ret = v_to_ret,
    p_surv = p_surv,
    adue_ret = adue_ret
  )
}

#' Traditional Unit Credit accrued liability
#'
#' Computes the actuarial present value of the benefit accrued through the
#' valuation date.
#'
#' @param accrued_benefit Nonnegative accrued benefit.
#' @inheritParams APV_NR_db
#'
#' @return A numeric vector.
#'
#' @examples
#' AAL_TUC_db(
#'   accrued_benefit = 12000,
#'   v_to_ret = 0.5,
#'   p_surv = 0.9,
#'   adue_ret = 12
#' )
#'
#' @export
AAL_TUC_db <- function(
    accrued_benefit,
    v_to_ret,
    p_surv,
    adue_ret) {
  APV_NR_db(
    PABz = accrued_benefit,
    v_to_ret = v_to_ret,
    p_surv = p_surv,
    adue_ret = adue_ret
  )
}

#' Projected Unit Credit normal cost
#'
#' Computes the actuarial present value of the portion of the projected
#' benefit attributed to the current year of service.
#'
#' @param projected_benefit Nonnegative projected annual benefit at retirement.
#' @param total_service Positive total service at retirement.
#' @inheritParams APV_NR_db
#'
#' @return A numeric vector.
#'
#' @examples
#' NC_PUC_db(
#'   projected_benefit = 30000,
#'   total_service = 30,
#'   v_to_ret = 0.5,
#'   p_surv = 0.9,
#'   adue_ret = 12
#' )
#'
#' @export
NC_PUC_db <- function(
    projected_benefit,
    total_service,
    v_to_ret,
    p_surv,
    adue_ret) {
  projected_benefit <- .pension_check_nonnegative(
    projected_benefit,
    "projected_benefit"
  )

  total_service <- .pension_check_positive(
    total_service,
    "total_service"
  )

  values <- .pension_recycle_common(
    projected_benefit,
    total_service,
    v_to_ret,
    p_surv,
    adue_ret,
    .names = c(
      "projected_benefit",
      "total_service",
      "v_to_ret",
      "p_surv",
      "adue_ret"
    )
  )

  APV_NR_db(
    PABz = values[[1L]] / values[[2L]],
    v_to_ret = values[[3L]],
    p_surv = values[[4L]],
    adue_ret = values[[5L]]
  )
}

#' Projected Unit Credit accrued liability
#'
#' Computes the actuarial present value of the portion of the projected
#' benefit attributed to past service.
#'
#' @param projected_benefit Nonnegative projected annual benefit at retirement.
#' @param past_service Nonnegative service completed through the valuation date.
#' @param total_service Positive total service at retirement.
#' @inheritParams APV_NR_db
#'
#' @return A numeric vector.
#'
#' @examples
#' AAL_PUC_db(
#'   projected_benefit = 30000,
#'   past_service = 10,
#'   total_service = 30,
#'   v_to_ret = 0.5,
#'   p_surv = 0.9,
#'   adue_ret = 12
#' )
#'
#' @export
AAL_PUC_db <- function(
    projected_benefit,
    past_service,
    total_service,
    v_to_ret,
    p_surv,
    adue_ret) {
  projected_benefit <- .pension_check_nonnegative(
    projected_benefit,
    "projected_benefit"
  )

  past_service <- .pension_check_nonnegative(
    past_service,
    "past_service"
  )

  total_service <- .pension_check_positive(
    total_service,
    "total_service"
  )

  values <- .pension_recycle_common(
    projected_benefit,
    past_service,
    total_service,
    v_to_ret,
    p_surv,
    adue_ret,
    .names = c(
      "projected_benefit",
      "past_service",
      "total_service",
      "v_to_ret",
      "p_surv",
      "adue_ret"
    )
  )

  projected_benefit <- values[[1L]]
  past_service <- values[[2L]]
  total_service <- values[[3L]]
  v_to_ret <- values[[4L]]
  p_surv <- values[[5L]]
  adue_ret <- values[[6L]]

  if (any(past_service > total_service)) {
    stop(
      "past_service cannot exceed total_service.",
      call. = FALSE
    )
  }

  APV_NR_db(
    PABz = projected_benefit *
      past_service /
      total_service,
    v_to_ret = v_to_ret,
    p_surv = p_surv,
    adue_ret = adue_ret
  )
}
