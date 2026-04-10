# =========================================================
# Chapter 18 helpers for mqriskR
# Pensions
# =========================================================
# Internal helpers ---------------------------------------------------------

.validate_scalar_numeric_ch18 <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be a numeric scalar.", call. = FALSE)
  }
}

.validate_nonneg_scalar_ch18 <- function(x, name) {
  .validate_scalar_numeric_ch18(x, name)
  if (x < 0) stop(name, " must be nonnegative.", call. = FALSE)
}

.validate_positive_scalar_ch18 <- function(x, name) {
  .validate_scalar_numeric_ch18(x, name)
  if (x <= 0) stop(name, " must be positive.", call. = FALSE)
}

.validate_prob_scalar_ch18 <- function(x, name) {
  .validate_scalar_numeric_ch18(x, name)
  if (x < 0 || x > 1) stop(name, " must be between 0 and 1.", call. = FALSE)
}

.salary_path_from_inputs_ch18 <- function(x, z, Sx, g = NULL, s = NULL) {
  .validate_scalar_numeric_ch18(x, "x")
  .validate_scalar_numeric_ch18(z, "z")
  .validate_positive_scalar_ch18(Sx, "Sx")

  if (z <= x) stop("z must be greater than x.", call. = FALSE)

  n <- z - x

  if (!is.null(g) && !is.null(s)) {
    stop("Supply only one of g or s.", call. = FALSE)
  }

  if (is.null(g) && is.null(s)) {
    stop("Supply one of g or s.", call. = FALSE)
  }

  if (!is.null(g)) {
    .validate_scalar_numeric_ch18(g, "g")
    return(Sx * (1 + g)^(0:(n - 1L)))
  }

  if (!is.numeric(s) || anyNA(s)) {
    stop("s must be a numeric vector.", call. = FALSE)
  }

  if (length(s) != n) {
    stop("s must have length z - x.", call. = FALSE)
  }

  if (any(s <= 0)) {
    stop("All elements of s must be positive.", call. = FALSE)
  }

  Sx * s / s[1L]
}

#' Salary scale under constant annual growth
#'
#' Constructs salary scale factors \eqn{s_k} under a constant annual growth rate.
#'
#' @param k Numeric vector of ages.
#' @param g Annual salary growth rate.
#' @param base_age Age at which the scale is normalized.
#' @param s_base Salary scale value at \code{base_age}.
#'
#' @return Numeric vector of salary scale factors.
#' @export
#'
#' @examples
#' salary_scale(k = 30:34, g = 0.04, base_age = 30)
salary_scale <- function(k, g, base_age = min(k), s_base = 1) {
  if (!is.numeric(k) || anyNA(k)) {
    stop("k must be a numeric vector.", call. = FALSE)
  }
  .validate_scalar_numeric_ch18(g, "g")
  .validate_scalar_numeric_ch18(base_age, "base_age")
  .validate_positive_scalar_ch18(s_base, "s_base")

  s_base * (1 + g)^(k - base_age)
}

#' Accumulated value of defined contribution plan contributions
#'
#' Computes the accumulated value at retirement age \eqn{z} for the defined
#' contribution model in Equation (18.1).
#'
#' @param x Entry age.
#' @param z Retirement age.
#' @param Sx Salary at age \code{x}.
#' @param c Contribution rate as a proportion of salary.
#' @param i Annual effective interest rate.
#' @param g Optional constant annual salary growth rate.
#' @param s Optional salary scale vector of length \code{z - x}.
#'
#' @return The accumulated value of contributions at age \code{z}.
#' @export
#'
#' @examples
#' AVz_dc(x = 30, z = 65, Sx = 50000, c = 0.10, i = 0.05, g = 0.04)
AVz_dc <- function(x, z, Sx, c, i, g = NULL, s = NULL) {
  .validate_prob_scalar_ch18(c, "c")
  .validate_scalar_numeric_ch18(i, "i")

  salary <- .salary_path_from_inputs_ch18(x = x, z = z, Sx = Sx, g = g, s = s)
  n <- z - x
  t <- 0:(n - 1L)

  sum(c * salary * (1 + i)^(n - t))
}

#' Retirement income from a defined contribution accumulation
#'
#' Converts a defined contribution accumulation to annual annuity-due income.
#'
#' @param AVz Accumulated value at retirement.
#' @param adue_z Whole life annuity-due factor at retirement age.
#'
#' @return Annual retirement income.
#' @export
#'
#' @examples
#' Income_dc(AVz = 824211.35, adue_z = 12)
Income_dc <- function(AVz, adue_z) {
  .validate_nonneg_scalar_ch18(AVz, "AVz")
  .validate_positive_scalar_ch18(adue_z, "adue_z")
  AVz / adue_z
}

#' Replacement ratio for a defined contribution plan
#'
#' Computes the replacement ratio defined in Equation (18.4).
#'
#' @param x Entry age.
#' @param z Retirement age.
#' @param Sx Salary at age \code{x}.
#' @param c Contribution rate.
#' @param i Annual effective interest rate.
#' @param adue_z Whole life annuity-due factor at age \code{z}.
#' @param g Optional constant annual salary growth rate.
#' @param s Optional salary scale vector of length \code{z - x}.
#'
#' @return Replacement ratio.
#' @export
#'
#' @examples
#' replacement_ratio_dc(
#'   x = 30, z = 65, Sx = 50000, c = 0.10, i = 0.05,
#'   adue_z = 12, g = 0.04
#' )
replacement_ratio_dc <- function(x, z, Sx, c, i, adue_z, g = NULL, s = NULL) {
  AVz <- AVz_dc(x = x, z = z, Sx = Sx, c = c, i = i, g = g, s = s)
  income <- Income_dc(AVz = AVz, adue_z = adue_z)
  final_salary <- .salary_path_from_inputs_ch18(x = x, z = z, Sx = Sx, g = g, s = s)[z - x]
  income / final_salary
}

#' Target contribution rate for a defined contribution plan
#'
#' Solves Equation (18.5) for the contribution rate required to achieve a
#' target replacement ratio.
#'
#' @param x Entry age.
#' @param z Retirement age.
#' @param Sx Salary at age \code{x}.
#' @param RR_target Target replacement ratio.
#' @param i Annual effective interest rate.
#' @param adue_z Whole life annuity-due factor at age \code{z}.
#' @param g Optional constant annual salary growth rate.
#' @param s Optional salary scale vector of length \code{z - x}.
#'
#' @return Required contribution rate.
#' @export
#'
#' @examples
#' contribution_rate_target(
#'   x = 30, z = 65, Sx = 60000, RR_target = 0.50,
#'   i = 0.06, adue_z = 11, g = 0.04
#' )
contribution_rate_target <- function(x, z, Sx, RR_target, i, adue_z, g = NULL, s = NULL) {
  .validate_prob_scalar_ch18(RR_target, "RR_target")
  .validate_scalar_numeric_ch18(i, "i")
  .validate_positive_scalar_ch18(adue_z, "adue_z")

  salary <- .salary_path_from_inputs_ch18(x = x, z = z, Sx = Sx, g = g, s = s)
  n <- z - x
  t <- 0:(n - 1L)

  numer <- RR_target * salary[n] * adue_z
  denom <- sum(salary * (1 + i)^(n - t))

  numer / denom
}

#' Projected annual benefit under a final average salary DB plan
#'
#' Computes the projected annual benefit for a final average salary plan
#' using Equation (18.6).
#'
#' @param x Current or entry age.
#' @param z Retirement age.
#' @param CASx Current annual salary at age \code{x}.
#' @param p Accrual percentage, e.g. \code{2} for 2 percent.
#' @param fas_years Number of years in the final average salary period.
#' @param past_service Past years of service already completed at age \code{x}.
#' @param g Optional constant annual salary growth rate.
#' @param s Optional salary scale vector of length \code{z - x}.
#'
#' @return Projected annual benefit.
#' @export
#'
#' @examples
#' PAB_fas(x = 35, z = 65, CASx = 60000, p = 2, fas_years = 3, g = 0.04)
PAB_fas <- function(x, z, CASx, p, fas_years = 3, past_service = 0, g = NULL, s = NULL) {
  .validate_scalar_numeric_ch18(p, "p")
  .validate_positive_scalar_ch18(fas_years, "fas_years")
  .validate_nonneg_scalar_ch18(past_service, "past_service")

  salary <- .salary_path_from_inputs_ch18(x = x, z = z, Sx = CASx, g = g, s = s)
  n <- z - x

  if (fas_years > n) {
    stop("fas_years cannot exceed z - x.", call. = FALSE)
  }

  FASz <- mean(tail(salary, fas_years))
  YOSz <- past_service + n

  0.01 * p * YOSz * FASz
}

#' Projected annual benefit under a career average earnings DB plan
#'
#' Computes the projected annual benefit for a CAE plan using Equation (18.9).
#'
#' @param x Current or entry age.
#' @param z Retirement age.
#' @param CASx Current annual salary at age \code{x}.
#' @param p Accrual percentage, e.g. \code{1} for 1 percent.
#' @param past_salary_total Optional total of actual past salaries.
#' @param g Optional constant annual salary growth rate.
#' @param s Optional salary scale vector of length \code{z - x}.
#'
#' @return Projected annual benefit.
#' @export
#'
#' @examples
#' PAB_cae(x = 30, z = 65, CASx = 100000, p = 1, g = 0.04)
PAB_cae <- function(x, z, CASx, p, past_salary_total = 0, g = NULL, s = NULL) {
  .validate_scalar_numeric_ch18(p, "p")
  .validate_nonneg_scalar_ch18(past_salary_total, "past_salary_total")

  salary <- .salary_path_from_inputs_ch18(x = x, z = z, Sx = CASx, g = g, s = s)
  PASz <- past_salary_total + sum(salary)

  0.01 * p * PASz
}

#' Replacement ratio for a defined benefit plan
#'
#' Computes a DB replacement ratio as benefit divided by a chosen salary measure.
#'
#' @param benefit Annual retirement benefit.
#' @param salary Salary measure used in the denominator.
#'
#' @return Replacement ratio.
#' @export
#'
#' @examples
#' replacement_ratio_db(benefit = 108008.66, salary = 187119.09)
replacement_ratio_db <- function(benefit, salary) {
  .validate_nonneg_scalar_ch18(benefit, "benefit")
  .validate_positive_scalar_ch18(salary, "salary")
  benefit / salary
}

#' Accrued benefit for a final average salary plan
#'
#' Computes the accrued benefit at the current date using service and salary
#' history only.
#'
#' @param salary_history Numeric vector of annual salaries to date.
#' @param p Accrual percentage, e.g. \code{2} for 2 percent.
#' @param fas_years Number of years in the final average salary average.
#'
#' @return Accrued benefit.
#' @export
#'
#' @examples
#' AB_fas(salary_history = c(150000, 156000), p = 1, fas_years = 2)
AB_fas <- function(salary_history, p, fas_years = 3) {
  if (!is.numeric(salary_history) || anyNA(salary_history) || any(salary_history <= 0)) {
    stop("salary_history must be a positive numeric vector.", call. = FALSE)
  }
  .validate_scalar_numeric_ch18(p, "p")
  .validate_positive_scalar_ch18(fas_years, "fas_years")

  yrs <- length(salary_history)
  FAS <- mean(tail(salary_history, min(fas_years, yrs)))

  0.01 * p * yrs * FAS
}

#' Accrued benefit for a career average earnings plan
#'
#' Computes the accrued benefit using actual salary history only.
#'
#' @param salary_history Numeric vector of annual salaries to date.
#' @param p Accrual percentage.
#'
#' @return Accrued benefit.
#' @export
#'
#' @examples
#' AB_cae(salary_history = c(100000, 104000, 108160), p = 1)
AB_cae <- function(salary_history, p) {
  if (!is.numeric(salary_history) || anyNA(salary_history) || any(salary_history <= 0)) {
    stop("salary_history must be a positive numeric vector.", call. = FALSE)
  }
  .validate_scalar_numeric_ch18(p, "p")

  0.01 * p * sum(salary_history)
}

#' APV of normal retirement benefit for a DB plan
#'
#' Computes Equation (18.10).
#'
#' @param PABz Projected annual benefit at retirement.
#' @param v_to_ret Discount factor from current age to retirement.
#' @param p_surv Active-service survival probability to retirement.
#' @param adue_ret Retirement annuity factor.
#'
#' @return Actuarial present value of the normal retirement benefit.
#' @export
#'
#' @examples
#' APV_NR_db(PABz = 108008.66, v_to_ret = 1 / 1.06^30, p_surv = 0.8, adue_ret = 12)
APV_NR_db <- function(PABz, v_to_ret, p_surv, adue_ret) {
  .validate_nonneg_scalar_ch18(PABz, "PABz")
  .validate_nonneg_scalar_ch18(v_to_ret, "v_to_ret")
  .validate_prob_scalar_ch18(p_surv, "p_surv")
  .validate_positive_scalar_ch18(adue_ret, "adue_ret")

  PABz * v_to_ret * p_surv * adue_ret
}

#' Entry Age Normal normal cost for a DB plan
#'
#' Computes Equation (18.13).
#'
#' @param APV_total Total actuarial present value of benefits.
#' @param adue_active Active-service annuity-due factor.
#'
#' @return Entry Age Normal normal cost.
#' @export
#'
#' @examples
#' NC_EAN_db(APV_total = 25000, adue_active = 15)
NC_EAN_db <- function(APV_total, adue_active) {
  .validate_nonneg_scalar_ch18(APV_total, "APV_total")
  .validate_positive_scalar_ch18(adue_active, "adue_active")

  APV_total / adue_active
}

#' Traditional Unit Credit normal cost for a DB plan
#'
#' Computes the TUC normal cost as the APV of the current year's accrual.
#'
#' @param accrual_benefit Benefit accrued in the current year.
#' @param v_to_ret Discount factor to retirement.
#' @param p_surv Active-service survival probability to retirement.
#' @param adue_ret Retirement annuity factor.
#'
#' @return TUC normal cost.
#' @export
#'
#' @examples
#' NC_TUC_db(accrual_benefit = 1560, v_to_ret = 0.5, p_surv = 0.9, adue_ret = 12)
NC_TUC_db <- function(accrual_benefit, v_to_ret, p_surv, adue_ret) {
  APV_NR_db(
    PABz = accrual_benefit,
    v_to_ret = v_to_ret,
    p_surv = p_surv,
    adue_ret = adue_ret
  )
}

#' Traditional Unit Credit accrued liability for a DB plan
#'
#' Computes the TUC accrued liability as the APV of the accrued benefit.
#'
#' @param accrued_benefit Accrued benefit at the valuation date.
#' @param v_to_ret Discount factor to retirement.
#' @param p_surv Active-service survival probability to retirement.
#' @param adue_ret Retirement annuity factor.
#'
#' @return TUC accrued liability.
#' @export
#'
#' @examples
#' AAL_TUC_db(accrued_benefit = 12000, v_to_ret = 0.5, p_surv = 0.9, adue_ret = 12)
AAL_TUC_db <- function(accrued_benefit, v_to_ret, p_surv, adue_ret) {
  APV_NR_db(
    PABz = accrued_benefit,
    v_to_ret = v_to_ret,
    p_surv = p_surv,
    adue_ret = adue_ret
  )
}

#' Projected Unit Credit normal cost for a DB plan
#'
#' Computes the PUC normal cost as the APV of the portion of the projected
#' benefit attributed to the current year of service.
#'
#' @param projected_benefit Projected benefit at retirement.
#' @param total_service Total service at retirement.
#' @param v_to_ret Discount factor to retirement.
#' @param p_surv Active-service survival probability to retirement.
#' @param adue_ret Retirement annuity factor.
#'
#' @return PUC normal cost.
#' @export
#'
#' @examples
#' NC_PUC_db(
#'   projected_benefit = 30000,
#'   total_service = 30,
#'   v_to_ret = 0.5,
#'   p_surv = 0.9,
#'   adue_ret = 12
#' )
NC_PUC_db <- function(projected_benefit, total_service, v_to_ret, p_surv, adue_ret) {
  .validate_nonneg_scalar_ch18(projected_benefit, "projected_benefit")
  .validate_positive_scalar_ch18(total_service, "total_service")

  APV_NR_db(
    PABz = projected_benefit / total_service,
    v_to_ret = v_to_ret,
    p_surv = p_surv,
    adue_ret = adue_ret
  )
}

#' Projected Unit Credit accrued liability for a DB plan
#'
#' Computes the PUC accrued liability as the APV of the portion of the
#' projected benefit attributed to past service.
#'
#' @param projected_benefit Projected benefit at retirement.
#' @param past_service Past service completed.
#' @param total_service Total service at retirement.
#' @param v_to_ret Discount factor to retirement.
#' @param p_surv Active-service survival probability to retirement.
#' @param adue_ret Retirement annuity factor.
#'
#' @return PUC accrued liability.
#' @export
#'
#' @examples
#' AAL_PUC_db(projected_benefit = 30000, past_service = 10, total_service = 30,
#' v_to_ret = 0.5, p_surv = 0.9, adue_ret = 12)
AAL_PUC_db <- function(projected_benefit, past_service, total_service, v_to_ret, p_surv, adue_ret) {
  .validate_nonneg_scalar_ch18(projected_benefit, "projected_benefit")
  .validate_nonneg_scalar_ch18(past_service, "past_service")
  .validate_positive_scalar_ch18(total_service, "total_service")

  if (past_service > total_service) {
    stop("past_service cannot exceed total_service.", call. = FALSE)
  }

  APV_NR_db(
    PABz = projected_benefit * (past_service / total_service),
    v_to_ret = v_to_ret,
    p_surv = p_surv,
    adue_ret = adue_ret
  )
}
