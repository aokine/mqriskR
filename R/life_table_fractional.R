# life_table_fractional.R
#
# Fractional-age life table utilities.
# These functions work from a discrete life table and a within-year assumption:
#   - "udd"       : uniform distribution of deaths
#   - "cf"        : constant force
#   - "balducci"  : hyperbolic / Balducci

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.match_assumption <- function(assumption) {
  assumption <- match.arg(assumption, choices = c("udd", "cf", "balducci"))
  assumption
}

.recycle_fractional_args <- function(...) {
  args <- lapply(list(...), as.numeric)
  lens <- vapply(args, length, integer(1))

  if (any(lens == 0L)) {
    stop("All arguments must have positive length.", call. = FALSE)
  }

  if (any(vapply(args, function(z) any(!is.finite(z)), logical(1)))) {
    stop("All arguments must contain finite numeric values.", call. = FALSE)
  }

  common_len <- max(lens)

  if (any(!(lens %in% c(1L, common_len)))) {
    stop("Arguments must have compatible lengths.", call. = FALSE)
  }

  lapply(args, function(z) {
    if (length(z) == 1L) rep(z, common_len) else z
  })
}

.check_fractional_inputs <- function(tbl, x, t) {
  .validate_life_table(tbl)

  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  if (any(t < 0) || any(t > 1)) {
    stop("t must satisfy 0 <= t <= 1.", call. = FALSE)
  }

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
#' under UDD, constant force, or Balducci.
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

  args <- .recycle_fractional_args(x, t)
  x <- args[[1L]]
  t <- args[[2L]]

  .check_fractional_inputs(tbl, x, t)

  qx <- qx_tab(tbl, x)
  px <- 1 - qx

  out <- switch(
    assumption,
    "udd" = 1 - t * qx,
    "cf" = px^t,
    "balducci" = {
      val <- px / (1 - (1 - t) * qx)
      val[t == 0] <- 1
      val
    }
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

  args <- .recycle_fractional_args(x, t)
  x <- args[[1L]]
  t <- args[[2L]]

  .check_fractional_inputs(tbl, x, t)

  qx <- qx_tab(tbl, x)
  px <- 1 - qx

  switch(
    assumption,
    "udd" = qx / (1 - t * qx),
    "cf" = -log(px),
    "balducci" = qx / (1 - (1 - t) * qx)
  )
}

# -------------------------------------------------------------------------
# Fractional conditional density within the year
# -------------------------------------------------------------------------

#' Fractional conditional density from a life table
#'
#' Computes the conditional density
#' \eqn{f_x(t \mid T_0 > x) = {}_t p_x \mu_{x+t}}.
#'
#' @inheritParams tpx_tab
#' @return Numeric vector of conditional density values.
#' @export
fx_tab <- function(tbl, x, t, assumption = c("udd", "cf", "balducci")) {
  assumption <- .match_assumption(assumption)

  args <- .recycle_fractional_args(x, t)
  x <- args[[1L]]
  t <- args[[2L]]

  .check_fractional_inputs(tbl, x, t)

  out <- tpx_tab(tbl, x, t, assumption = assumption) *
    mux_tab(tbl, x, t, assumption = assumption)

  pmax(out, 0)
}

# -------------------------------------------------------------------------
# Deferred probabilities
# -------------------------------------------------------------------------

#' Deferred death probability from a life table
#'
#' Computes the probability that a life aged x survives n years and then
#' dies within the following m years:
#' \eqn{{}_{n|m} q_x = {}_n p_x \cdot {}_m q_{x+n}}.
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

  args <- .recycle_fractional_args(x, n, m)
  x <- args[[1L]]
  n <- args[[2L]]
  m <- args[[3L]]

  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  if (any(n < 0 | n != floor(n))) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }

  if (any(m < 0 | m != floor(m))) {
    stop("m must be a nonnegative integer.", call. = FALSE)
  }

  npx(tbl, x, n) * nqx(tbl, x + n, m)
}

#' Curtate death probability from a life table
#'
#' Computes \eqn{{}_{k|} q_x = {}_k p_x - {}_{k+1} p_x}.
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of ages.
#' @param k Nonnegative integer.
#'
#' @return Numeric vector of \eqn{{}_{k|} q_x} values.
#' @export
nkqx <- function(tbl, x, k) {
  .validate_life_table(tbl)

  args <- .recycle_fractional_args(x, k)
  x <- args[[1L]]
  k <- args[[2L]]

  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  if (any(k < 0 | k != floor(k))) {
    stop("k must be a nonnegative integer.", call. = FALSE)
  }

  npx(tbl, x, k) - npx(tbl, x, k + 1)
}
