# life_table_fractional.R
#
# Fractional-age life table utilities for Chapter 6 of Models for Quantifying Risk.
# These functions work from a discrete life table and a within-year assumption:
#   - "udd"       : uniform distribution of deaths
#   - "cf"        : constant force
#   - "balducci"  : hyperbolic / Balducci
#
# Main quantities:
#   {}_t p_x, {}_t q_x, mu_{x+t}, f_x(t | T0 > x)

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.match_assumption <- function(assumption) {
  assumption <- tolower(assumption[1])
  ok <- assumption %in% c("udd", "cf", "balducci")
  if (!ok) {
    stop("assumption must be one of: 'udd', 'cf', 'balducci'.", call. = FALSE)
  }
  assumption
}

.recycle_xt <- function(x, t) {
  x <- as.numeric(x)
  t <- as.numeric(t)

  if (length(x) == 0L || length(t) == 0L) {
    stop("x and t must have positive length.", call. = FALSE)
  }

  if (length(x) == 1L && length(t) > 1L) x <- rep(x, length(t))
  if (length(t) == 1L && length(x) > 1L) t <- rep(t, length(x))

  if (length(x) != length(t)) {
    stop("x and t must have the same length, or one must have length 1.", call. = FALSE)
  }

  list(x = x, t = t)
}

.check_fractional_inputs <- function(tbl, x, t) {
  .validate_life_table(tbl)

  if (any(!is.finite(x)) || any(!is.finite(t))) {
    stop("x and t must be finite.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("x must be >= 0.", call. = FALSE)
  }
  if (any(t < 0) || any(t > 1)) {
    stop("t must satisfy 0 <= t <= 1.", call. = FALSE)
  }

  # Need q_x available at these integer ages
  qx_vals <- qx_tab(tbl, x)
  if (any(is.na(qx_vals))) {
    stop("Some x values are not available in the life table.", call. = FALSE)
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Fractional-age survival probabilities
# -------------------------------------------------------------------------

#' Fractional survival probability from a life table
#'
#' Computes \eqn{{}_t p_x} for \eqn{0 \le t \le 1} from a discrete life table
#' under one of the standard Chapter 6 assumptions:
#' UDD, constant force, or Balducci.
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of integer ages.
#' @param t Numeric vector of fractional durations with \eqn{0 \le t \le 1}.
#' @param assumption One of \code{"udd"}, \code{"cf"}, \code{"balducci"}.
#'
#' @return Numeric vector of \eqn{{}_t p_x} values.
#' @export
tpx_tab <- function(tbl, x, t, assumption = c("udd", "cf", "balducci")) {
  assumption <- .match_assumption(assumption)
  xt <- .recycle_xt(x, t)
  x <- xt$x
  t <- xt$t

  .check_fractional_inputs(tbl, x, t)

  qx <- qx_tab(tbl, x)
  px <- 1 - qx

  out <- switch(
    assumption,
    "udd" = 1 - t * qx,
    "cf" = px^t,
    "balducci" = (1 - qx) / (1 - (1 - t) * qx)
  )

  pmin(pmax(out, 0), 1)
}

#' Fractional failure probability from a life table
#'
#' Computes \eqn{{}_t q_x = 1 - {}_t p_x} for \eqn{0 \le t \le 1}.
#'
#' @inheritParams tpx_tab
#' @return Numeric vector of \eqn{{}_t q_x} values.
#' @export
tqx_tab <- function(tbl, x, t, assumption = c("udd", "cf", "balducci")) {
  1 - tpx_tab(tbl, x, t, assumption = assumption)
}

# -------------------------------------------------------------------------
# Force of mortality / hazard within the year
# -------------------------------------------------------------------------

#' Fractional force of mortality from a life table
#'
#' Computes \eqn{\mu_{x+t}} under UDD, constant force, or Balducci.
#'
#' @inheritParams tpx_tab
#' @return Numeric vector of \eqn{\mu_{x+t}} values.
#' @export
mux_tab <- function(tbl, x, t, assumption = c("udd", "cf", "balducci")) {
  assumption <- .match_assumption(assumption)
  xt <- .recycle_xt(x, t)
  x <- xt$x
  t <- xt$t

  .check_fractional_inputs(tbl, x, t)

  qx <- qx_tab(tbl, x)
  px <- 1 - qx

  out <- switch(
    assumption,
    "udd" = qx / (1 - t * qx),
    "cf" = -log(px),
    "balducci" = qx / (1 - (1 - t) * qx)
  )

  out
}

# -------------------------------------------------------------------------
# Fractional conditional density within the year
# -------------------------------------------------------------------------

#' Fractional conditional density from a life table
#'
#' Computes the conditional density
#' \eqn{f_x(t \mid T_0 > x) = {}_t p_x \mu_{x+t}}
#' for \eqn{0 < t < 1}.
#'
#' @inheritParams tpx_tab
#' @return Numeric vector of conditional density values.
#' @export
fx_tab <- function(tbl, x, t, assumption = c("udd", "cf", "balducci")) {
  assumption <- .match_assumption(assumption)
  xt <- .recycle_xt(x, t)
  x <- xt$x
  t <- xt$t

  .check_fractional_inputs(tbl, x, t)

  out <- tpx_tab(tbl, x, t, assumption = assumption) *
    mux_tab(tbl, x, t, assumption = assumption)

  pmax(out, 0)
}

# -------------------------------------------------------------------------
# Deferred probabilities: {}_{n|m} q_x and {}_{k|} q_x
# -------------------------------------------------------------------------

#' Deferred probability {}_{n|m} q_x from a life table
#'
#' Computes the probability that a life aged x survives n years and then
#' dies within the following m years:
#' \eqn{{}_{n|m} q_x = {}_n p_x \cdot {}_m q_{x+n}}
#'
#' This function is for integer n and m in the discrete tabular setting.
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of ages.
#' @param n Nonnegative integer deferred period.
#' @param m Nonnegative integer subsequent period.
#'
#' @return Numeric vector of \eqn{{}_{n|m} q_x} values.
#' @export
nmxq <- function(tbl, x, n, m) {
  .validate_life_table(tbl)

  len <- max(length(x), length(n), length(m))
  x <- rep_len(x, len)
  n <- rep_len(n, len)
  m <- rep_len(m, len)

  if (any(n < 0 | n != floor(n), na.rm = TRUE)) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }
  if (any(m < 0 | m != floor(m), na.rm = TRUE)) {
    stop("m must be a nonnegative integer.", call. = FALSE)
  }

  npx(tbl, x, n) * nqx(tbl, x + n, m)
}

#' Curtate death probability {}_{k|} q_x from a life table
#'
#' Computes
#' \eqn{{}_{k|} q_x = {}_k p_x - {}_{k+1} p_x}
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of ages.
#' @param k Nonnegative integer.
#'
#' @return Numeric vector of \eqn{{}_{k|} q_x} values.
#' @export
nkqx <- function(tbl, x, k) {
  .validate_life_table(tbl)

  len <- max(length(x), length(k))
  x <- rep_len(x, len)
  k <- rep_len(k, len)

  if (any(k < 0 | k != floor(k), na.rm = TRUE)) {
    stop("k must be a nonnegative integer.", call. = FALSE)
  }

  npx(tbl, x, k) - npx(tbl, x, k + 1)
}
