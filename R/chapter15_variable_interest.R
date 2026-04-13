# =========================================================
# Chapter 15 helpers for mqriskR
# Models with Variable Interest Rates
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.validate_prob_vector_ch15 <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite numeric vector.", name), call. = FALSE)
  }
  if (any(x < 0 | x > 1)) {
    stop(sprintf("%s must contain values in [0, 1].", name), call. = FALSE)
  }
  invisible(x)
}

.validate_nonneg_numeric_ch15 <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite numeric vector.", name), call. = FALSE)
  }
  if (any(x < 0)) {
    stop(sprintf("%s must be nonnegative.", name), call. = FALSE)
  }
  invisible(x)
}

.surv_start_from_qx_ch15 <- function(qx) {
  n <- length(qx)
  c(1, cumprod(1 - qx))[1:n]
}

.surv_end_from_qx_ch15 <- function(qx) {
  cumprod(1 - qx)
}

# -------------------------------------------------------------------------
# Variable-interest discount factors
# -------------------------------------------------------------------------

#' Discount factors under a variable annual interest scenario
#'
#' Computes the sequence of discount factors
#' \deqn{v_1,\; v_2,\; \dots,\; v_n}
#' where
#' \deqn{v_t = \prod_{k=1}^{t}(1+i_k)^{-1}.}
#'
#' This corresponds to the Chapter 15 notation \eqn{{}_j v^t} for a fixed
#' scenario \eqn{j}.
#'
#' @param i Numeric vector of annual effective interest rates
#'   \eqn{i_1, i_2, \dots, i_n}.
#'
#' @return Numeric vector of discount factors of the same length as \code{i}.
#'
#' @examples
#' vt_var(c(0.06, 0.07, 0.08))
#'
#' @export
vt_var <- function(i) {
  .validate_nonneg_numeric_ch15(i, "i")
  cumprod((1 + i)^(-1))
}

# -------------------------------------------------------------------------
# APVs with variable annual interest scenarios
# -------------------------------------------------------------------------

#' Variable-interest actuarial present value functions
#'
#' Chapter 15 functions for actuarial present values under variable annual
#' effective interest rates interpreted as a yearly scenario
#' \eqn{i_1, i_2, \dots, i_n}.
#'
#' @name chapter15_variable_interest_apv
NULL

#' Pure endowment APV under variable annual interest rates
#'
#' Computes the APV of an \eqn{n}-year pure endowment under a variable annual
#' interest scenario:
#' \deqn{
#' {}_nE_x = v_n \cdot {}_np_x.
#' }
#'
#' If a benefit amount is supplied, the function returns that benefit times the
#' APV factor.
#'
#' @param qx Numeric vector of one-year mortality rates
#'   \eqn{q_x, q_{x+1}, \dots, q_{x+n-1}}.
#' @param i Numeric vector of annual effective interest rates
#'   \eqn{i_1, i_2, \dots, i_n}.
#' @param benefit Benefit amount payable at the end of year \eqn{n}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- c(.03, .04, .05, .06, .07)
#' nEx_var(qx = qx, i = c(.06, .07, .08, .09, .10), benefit = 1000)
#'
#' @rdname chapter15_variable_interest_apv
#' @aliases nEx_var
#' @export
nEx_var <- function(qx, i, benefit = 1) {
  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(i, "i")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(i)) {
    stop("qx and i must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  benefit * prod(1 - qx) * tail(vt_var(i), 1L)
}

#' Term insurance APV under variable annual interest rates
#'
#' Computes the APV of an \eqn{n}-year term insurance with benefit paid at the
#' end of the year of death under a variable annual interest scenario:
#' \deqn{
#' A_{x:\overline{n}|}^1 = \sum_{t=1}^{n} v_t \cdot {}_{t-1}p_x \cdot q_{x+t-1}.
#' }
#'
#' If a benefit amount is supplied, the function returns that benefit times the
#' APV factor.
#'
#' @param qx Numeric vector of one-year mortality rates
#'   \eqn{q_x, q_{x+1}, \dots, q_{x+n-1}}.
#' @param i Numeric vector of annual effective interest rates
#'   \eqn{i_1, i_2, \dots, i_n}.
#' @param benefit Benefit amount payable at the end of the year of death.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- c(.03, .04, .05, .06, .07)
#' Axn1_var(qx = qx, i = c(.06, .07, .08, .09, .10))
#'
#' @rdname chapter15_variable_interest_apv
#' @aliases Axn1_var
#' @export
Axn1_var <- function(qx, i, benefit = 1) {
  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(i, "i")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(i)) {
    stop("qx and i must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  vt <- vt_var(i)
  surv_start <- .surv_start_from_qx_ch15(qx)

  benefit * sum(vt * surv_start * qx)
}

#' Endowment insurance APV under variable annual interest rates
#'
#' Computes the APV of an \eqn{n}-year endowment insurance under a variable
#' annual interest scenario:
#' \deqn{
#' A_{x:\overline{n}|} = A_{x:\overline{n}|}^1 + {}_nE_x.
#' }
#'
#' If a benefit amount is supplied, the function returns that benefit times the
#' APV factor.
#'
#' @param qx Numeric vector of one-year mortality rates
#'   \eqn{q_x, q_{x+1}, \dots, q_{x+n-1}}.
#' @param i Numeric vector of annual effective interest rates
#'   \eqn{i_1, i_2, \dots, i_n}.
#' @param benefit Benefit amount payable either at the end of the year of death
#'   during the term or at the end of year \eqn{n} on survival.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- c(.03, .04, .05, .06, .07)
#' Axn_var(qx = qx, i = c(.06, .07, .08, .09, .10))
#'
#' @rdname chapter15_variable_interest_apv
#' @aliases Axn_var
#' @export
Axn_var <- function(qx, i, benefit = 1) {
  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(i, "i")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(i)) {
    stop("qx and i must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  Axn1_var(qx = qx, i = i, benefit = benefit) +
    nEx_var(qx = qx, i = i, benefit = benefit)
}

#' Temporary annuity APV under variable annual interest rates
#'
#' Computes the APV of an \eqn{n}-year temporary life annuity under a variable
#' annual interest scenario.
#'
#' For an immediate annuity,
#' \deqn{
#' a_{x:\overline{n}|} = \sum_{t=1}^{n} v_t \cdot {}_tp_x.
#' }
#'
#' For an annuity-due,
#' \deqn{
#' \ddot{a}_{x:\overline{n}|} = \sum_{t=0}^{n-1} v_t \cdot {}_tp_x,
#' }
#' with \eqn{v_0 = 1}.
#'
#' @param qx Numeric vector of one-year mortality rates
#'   \eqn{q_x, q_{x+1}, \dots, q_{x+n-1}}.
#' @param i Numeric vector of annual effective interest rates
#'   \eqn{i_1, i_2, \dots, i_n}.
#' @param type Either \code{"immediate"} or \code{"due"}.
#' @param benefit Amount of each annuity payment.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- rep(.02, 5)
#' axn_var(qx = qx, i = c(.06, .05, .04, .03, .03), type = "immediate")
#' axn_var(qx = qx, i = c(.03, .04, .05, .06, .07), type = "due")
#'
#' @rdname chapter15_variable_interest_apv
#' @aliases axn_var
#' @export
axn_var <- function(qx, i, type = c("immediate", "due"), benefit = 1) {
  type <- match.arg(type)

  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(i, "i")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(i)) {
    stop("qx and i must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  vt <- vt_var(i)

  if (type == "immediate") {
    surv_end <- .surv_end_from_qx_ch15(qx)
    return(benefit * sum(vt * surv_end))
  }

  surv_start <- .surv_start_from_qx_ch15(qx)
  vt_due <- c(1, vt[-length(vt)])
  benefit * sum(vt_due * surv_start)
}

# -------------------------------------------------------------------------
# APVs with spot rates
# -------------------------------------------------------------------------

#' Spot-rate actuarial present value functions
#'
#' Chapter 15 functions for actuarial present values when discounting uses
#' spot rates by maturity. If \eqn{z_t} denotes the annual effective spot rate
#' for maturity \eqn{t}, then the discount factor is \eqn{(1+z_t)^{-t}}.
#'
#' @name chapter15_spot_interest_apv
NULL

#' Pure endowment APV using spot rates
#'
#' Computes
#' \deqn{
#' {}_nE_x = (1+z_n)^{-n}\cdot {}_np_x.
#' }
#'
#' @param qx Numeric vector of one-year mortality rates.
#' @param z Numeric vector of annual effective spot rates for maturities
#'   \eqn{1,\dots,n}.
#' @param benefit Benefit amount payable at the end of year \eqn{n}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- c(.02, .03, .04, .05, .06)
#' z  <- c(.03, .04, .05, .06, .07)
#' nEx_spot(qx, z, benefit = 1000000)
#'
#' @rdname chapter15_spot_interest_apv
#' @aliases nEx_spot
#' @export
nEx_spot <- function(qx, z, benefit = 1) {
  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(z, "z")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(z)) {
    stop("qx and z must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  n <- length(z)
  benefit * prod(1 - qx) / (1 + z[n])^n
}

#' Term insurance APV using spot rates
#'
#' Computes
#' \deqn{
#' A_{x:\overline{n}|}^1 = \sum_{t=1}^{n}(1+z_t)^{-t}\cdot {}_{t-1}p_x \cdot q_{x+t-1}.
#' }
#'
#' @param qx Numeric vector of one-year mortality rates.
#' @param z Numeric vector of annual effective spot rates for maturities
#'   \eqn{1,\dots,n}.
#' @param benefit Benefit amount payable at the end of the year of death.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- c(.02, .03, .04, .05, .06)
#' z  <- c(.03, .04, .05, .06, .07)
#' Axn1_spot(qx, z)
#'
#' @rdname chapter15_spot_interest_apv
#' @aliases Axn1_spot
#' @export
Axn1_spot <- function(qx, z, benefit = 1) {
  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(z, "z")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(z)) {
    stop("qx and z must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  n <- length(z)
  vt <- (1 + z)^(-seq_len(n))
  surv_start <- .surv_start_from_qx_ch15(qx)

  benefit * sum(vt * surv_start * qx)
}

#' Endowment insurance APV using spot rates
#'
#' Computes
#' \deqn{
#' A_{x:\overline{n}|} = A_{x:\overline{n}|}^1 + {}_nE_x
#' }
#' using spot-rate discount factors.
#'
#' @param qx Numeric vector of one-year mortality rates.
#' @param z Numeric vector of annual effective spot rates for maturities
#'   \eqn{1,\dots,n}.
#' @param benefit Benefit amount.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- c(.02, .03, .04, .05, .06)
#' z  <- c(.03, .04, .05, .06, .07)
#' Axn_spot(qx, z)
#'
#' @rdname chapter15_spot_interest_apv
#' @aliases Axn_spot
#' @export
Axn_spot <- function(qx, z, benefit = 1) {
  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(z, "z")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(z)) {
    stop("qx and z must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  Axn1_spot(qx = qx, z = z, benefit = benefit) +
    nEx_spot(qx = qx, z = z, benefit = benefit)
}

#' Temporary annuity using spot rates
#'
#' For an immediate annuity,
#' \deqn{
#' a_{x:\overline{n}|} = \sum_{t=1}^{n}(1+z_t)^{-t}\cdot {}_tp_x.
#' }
#'
#' For an annuity-due,
#' \deqn{
#' \ddot{a}_{x:\overline{n}|} = \sum_{t=0}^{n-1}(1+z_t)^{-t}\cdot {}_tp_x,
#' }
#' where the time-0 discount factor is 1.
#'
#' @param qx Numeric vector of one-year mortality rates.
#' @param z Numeric vector of annual effective spot rates for maturities
#'   \eqn{1,\dots,n}.
#' @param type Either \code{"immediate"} or \code{"due"}.
#' @param benefit Amount of each annuity payment.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qx <- c(.02, .03, .04, .05, .06)
#' z  <- c(.03, .04, .05, .06, .07)
#' axn_spot(qx, z, type = "due")
#'
#' @rdname chapter15_spot_interest_apv
#' @aliases axn_spot
#' @export
axn_spot <- function(qx, z, type = c("immediate", "due"), benefit = 1) {
  type <- match.arg(type)

  .validate_prob_vector_ch15(qx, "qx")
  .validate_nonneg_numeric_ch15(z, "z")
  .validate_nonneg_numeric_ch15(benefit, "benefit")

  if (length(qx) != length(z)) {
    stop("qx and z must have the same length.", call. = FALSE)
  }
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }

  n <- length(z)

  if (type == "immediate") {
    surv_end <- .surv_end_from_qx_ch15(qx)
    vt <- (1 + z)^(-seq_len(n))
    return(benefit * sum(vt * surv_end))
  }

  surv_start <- .surv_start_from_qx_ch15(qx)
  vt_due <- c(1, (1 + z[-n])^(-seq_len(n - 1L)))
  benefit * sum(vt_due * surv_start)
}

# -------------------------------------------------------------------------
# Spot-rate valuation and bootstrapping
# -------------------------------------------------------------------------

#' Present value of cash flows using spot rates
#'
#' Discounts deterministic cash flows using spot rates matched to their
#' maturities.
#'
#' For annual compounding, each positive-time cash flow at time \eqn{t} is
#' discounted by \eqn{(1+z_t)^{-t}}.
#'
#' For semiannual nominal compounding, each positive-time cash flow at time
#' \eqn{t} is discounted by \eqn{(1+z_t/2)^{-2t}}.
#'
#' Time-0 cash flows are left undiscounted.
#'
#' @param amounts Numeric vector of cash flow amounts.
#' @param times Numeric vector of payment times in years.
#' @param spot Numeric vector of spot rates matched elementwise to
#'   \code{times}. Use 0 for any time-0 entry.
#' @param compounding Either \code{"annual"} or \code{"semiannual"}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' pv_spot_cashflows(
#'   amounts = c(200000, 50000, 50000, 100000),
#'   times   = c(0, 0.5, 1, 2),
#'   spot    = c(0, 0.02440, 0.02601, 0.02936),
#'   compounding = "semiannual"
#' )
#'
#' @export
pv_spot_cashflows <- function(amounts, times, spot,
                              compounding = c("annual", "semiannual")) {
  compounding <- match.arg(compounding)

  if (!is.numeric(amounts) || !is.numeric(times) || !is.numeric(spot)) {
    stop("amounts, times, and spot must be numeric vectors.", call. = FALSE)
  }
  if (!(length(amounts) == length(times) && length(times) == length(spot))) {
    stop("amounts, times, and spot must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(amounts)) || any(!is.finite(times)) || any(!is.finite(spot))) {
    stop("amounts, times, and spot must be finite.", call. = FALSE)
  }
  if (any(times < 0) || any(spot < 0)) {
    stop("times and spot rates must be nonnegative.", call. = FALSE)
  }

  df <- rep(1, length(times))
  pos <- times > 0

  if (compounding == "annual") {
    df[pos] <- (1 + spot[pos])^(-times[pos])
  } else {
    df[pos] <- (1 + spot[pos] / 2)^(-2 * times[pos])
  }

  sum(amounts * df)
}

#' Bootstrap semiannual nominal spot rates from coupon-bond yields
#'
#' Bootstraps the semiannual nominal annual zero-coupon yields from par
#' coupon-bearing bond yields of the same maturities.
#'
#' Both coupon yields and spot yields are interpreted as nominal annual rates
#' convertible semiannually.
#'
#' @param maturity Numeric vector of maturities in years, typically
#'   \code{0.5, 1.0, 1.5, ...}, in increasing order.
#' @param coupon_yield Numeric vector of nominal annual coupon yields
#'   convertible semiannually.
#' @param par Par value of each bond.
#'
#' @return Numeric vector of semiannual nominal annual spot rates.
#'
#' @examples
#' maturity <- c(0.5, 1.0, 1.5, 2.0)
#' coupon_yield <- c(0.0244, 0.0260, 0.0276, 0.0293)
#' z_from_coupon_semi(maturity, coupon_yield)
#'
#' @export
z_from_coupon_semi <- function(maturity, coupon_yield, par = 1000) {
  if (!is.numeric(maturity) || !is.numeric(coupon_yield)) {
    stop("maturity and coupon_yield must be numeric vectors.", call. = FALSE)
  }
  if (length(maturity) != length(coupon_yield)) {
    stop("maturity and coupon_yield must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(maturity)) || any(!is.finite(coupon_yield))) {
    stop("maturity and coupon_yield must be finite.", call. = FALSE)
  }
  if (any(maturity <= 0) || is.unsorted(maturity, strictly = TRUE)) {
    stop("maturity must be positive and strictly increasing.", call. = FALSE)
  }
  if (any(coupon_yield < 0) || !is.numeric(par) || length(par) != 1L || par <= 0) {
    stop("coupon_yield must be nonnegative and par must be a positive scalar.", call. = FALSE)
  }

  periods <- round(2 * maturity)
  if (any(abs(2 * maturity - periods) > 1e-10)) {
    stop("For semiannual bootstrapping, maturities must be multiples of 0.5 years.", call. = FALSE)
  }

  z <- numeric(length(maturity))

  for (j in seq_along(maturity)) {
    cpn <- par * coupon_yield[j] / 2
    m <- periods[j]

    if (j == 1L) {
      z[j] <- coupon_yield[j]
    } else {
      pv_prior <- 0
      for (r in seq_len(m - 1L)) {
        idx <- which(periods == r)
        if (length(idx) != 1L) {
          stop("Missing earlier maturity needed for bootstrapping.", call. = FALSE)
        }
        pv_prior <- pv_prior + cpn / (1 + z[idx] / 2)^r
      }

      df_m <- (par + cpn) / (par - pv_prior)
      z[j] <- 2 * (df_m^(1 / m) - 1)
    }
  }

  z
}

#' Bootstrap annual spot rates from annual coupon-bond yields
#'
#' Bootstraps annual effective zero-coupon yields from par annual coupon-bearing
#' bond yields of the same maturities.
#'
#' @param maturity Integer vector of maturities in years, in increasing order.
#' @param coupon_yield Numeric vector of annual coupon yields.
#' @param par Par value of each bond.
#'
#' @return Numeric vector of annual effective spot rates.
#'
#' @examples
#' maturity <- 1:4
#' coupon_yield <- c(0.02, 0.04, 0.06, 0.08)
#' z_from_coupon_annual(maturity, coupon_yield)
#'
#' @export
z_from_coupon_annual <- function(maturity, coupon_yield, par = 1000) {
  if (!is.numeric(maturity) || !is.numeric(coupon_yield)) {
    stop("maturity and coupon_yield must be numeric vectors.", call. = FALSE)
  }
  if (length(maturity) != length(coupon_yield)) {
    stop("maturity and coupon_yield must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(maturity)) || any(!is.finite(coupon_yield))) {
    stop("maturity and coupon_yield must be finite.", call. = FALSE)
  }
  if (any(maturity <= 0) || is.unsorted(maturity, strictly = TRUE)) {
    stop("maturity must be positive and strictly increasing.", call. = FALSE)
  }
  if (any(abs(maturity - round(maturity)) > 1e-10)) {
    stop("For annual bootstrapping, maturities must be integers.", call. = FALSE)
  }
  if (any(coupon_yield < 0) || !is.numeric(par) || length(par) != 1L || par <= 0) {
    stop("coupon_yield must be nonnegative and par must be a positive scalar.", call. = FALSE)
  }

  maturity <- as.integer(round(maturity))
  z <- numeric(length(maturity))

  for (j in seq_along(maturity)) {
    cpn <- par * coupon_yield[j]
    m <- maturity[j]

    if (j == 1L) {
      z[j] <- coupon_yield[j]
    } else {
      pv_prior <- 0
      for (r in seq_len(m - 1L)) {
        idx <- which(maturity == r)
        if (length(idx) != 1L) {
          stop("Missing earlier maturity needed for bootstrapping.", call. = FALSE)
        }
        pv_prior <- pv_prior + cpn / (1 + z[idx])^r
      }

      df_m <- (par + cpn) / (par - pv_prior)
      z[j] <- df_m^(1 / m) - 1
    }
  }

  z
}

# -------------------------------------------------------------------------
# Forward and spot rate conversions
# -------------------------------------------------------------------------

#' Forward rate \eqn{f_{n,k}} from spot rates
#'
#' Computes the \eqn{n}-year forward \eqn{k}-year annual effective rate implied
#' by annual effective spot rates:
#' \deqn{
#' (1+z_{n+k})^{n+k} = (1+z_n)^n (1+f_{n,k})^k.
#' }
#'
#' The input vector \code{z} should contain annual effective spot rates for
#' maturities 1, 2, ..., length(z).
#'
#' @param z Numeric vector of annual effective spot rates.
#' @param n Forward start in years.
#' @param k Forward maturity in years.
#'
#' @return A numeric scalar.
#'
#' @examples
#' z <- c(0.03, 0.04, 0.05, 0.06, 0.07)
#' fnk_from_z(z, n = 1, k = 4)
#' fnk_from_z(z, n = 2, k = 2)
#'
#' @export
fnk_from_z <- function(z, n, k) {
  .validate_nonneg_numeric_ch15(z, "z")

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || abs(n - round(n)) > 1e-10) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0 || abs(k - round(k)) > 1e-10) {
    stop("k must be a positive integer.", call. = FALSE)
  }

  n <- as.integer(round(n))
  k <- as.integer(round(k))

  if (n + k > length(z)) {
    stop("Need spot rates through maturity n + k.", call. = FALSE)
  }

  zn <- if (n == 0L) 0 else z[n]
  znk <- z[n + k]

  ((1 + znk)^(n + k) / (1 + zn)^n)^(1 / k) - 1
}

#' Matrix of all determinable forward rates from spot rates
#'
#' Constructs an upper-left triangular matrix of annual effective forward rates
#' \eqn{f_{n,k}} implied by annual effective spot rates \eqn{z_1,\dots,z_m}.
#'
#' Rows correspond to \eqn{n = 1,\dots,m-1} and columns correspond to
#' \eqn{k = 1,\dots,m-1}. Entries that are not determinable are returned as
#' \code{NA}.
#'
#' @param z Numeric vector of annual effective spot rates.
#'
#' @return A numeric matrix.
#'
#' @examples
#' z <- c(0.03, 0.04, 0.05, 0.06, 0.07)
#' forward_matrix_from_z(z)
#'
#' @export
forward_matrix_from_z <- function(z) {
  .validate_nonneg_numeric_ch15(z, "z")

  m <- length(z)
  out <- matrix(NA_real_, nrow = m - 1L, ncol = m - 1L)

  for (n in 1:(m - 1L)) {
    for (k in 1:(m - n)) {
      out[n, k] <- fnk_from_z(z, n = n, k = k)
    }
  }

  rownames(out) <- paste0("n=", 1:(m - 1L))
  colnames(out) <- paste0("k=", 1:(m - 1L))
  out
}

#' Spot rates from forward one-year rates
#'
#' Converts annual effective forward one-year rates
#' \eqn{f_{0,1}, f_{1,1}, \dots, f_{n-1,1}}
#' into annual effective spot rates
#' \eqn{z_1, z_2, \dots, z_n}.
#'
#' Since
#' \deqn{
#' (1+z_n)^n = \prod_{j=0}^{n-1}(1+f_{j,1}),
#' }
#' the spot rates are recovered directly.
#'
#' @param fn1 Numeric vector of annual effective forward one-year rates.
#'
#' @return Numeric vector of annual effective spot rates.
#'
#' @examples
#' z_from_fn1(c(0.04, 0.05, 0.06, 0.07, 0.08))
#'
#' @export
z_from_fn1 <- function(fn1) {
  .validate_nonneg_numeric_ch15(fn1, "fn1")

  gross <- cumprod(1 + fn1)
  gross^(1 / seq_along(fn1)) - 1
}
