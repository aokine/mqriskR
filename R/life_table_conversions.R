# life_table_conversions.R
#
# Conversion utilities linking Chapter 5 survival-model notation
# and Chapter 6 life-table notation.
#
# Functions:
#   S0_to_lx()  : convert S0(x) values to l_x values using a radix
#   lx_to_S0()  : convert l_x values to S0(x) values
#   px_to_lx()  : construct l_x values from p_x values and a radix
#   qx_to_lx()  : construct l_x values from q_x values and a radix

# -------------------------------------------------------------------------
# Internal helper
# -------------------------------------------------------------------------

.check_radix <- function(radix) {
  if (!is.numeric(radix) || length(radix) != 1L || !is.finite(radix) || radix <= 0) {
    stop("radix must be a single positive finite number.", call. = FALSE)
  }
  radix
}

# -------------------------------------------------------------------------
# S0 -> lx
# -------------------------------------------------------------------------

#' Convert survival probabilities to life-table values
#'
#' Converts Chapter 5 survival function values \eqn{S_0(x)} into
#' Chapter 6 life-table values \eqn{l_x = l_0 S_0(x)} using a chosen radix.
#'
#' @param S0 Numeric vector of survival probabilities.
#' @param radix Positive radix \eqn{l_0}.
#'
#' @return Numeric vector of \eqn{l_x} values.
#' @export
S0_to_lx <- function(S0, radix = 100000) {
  radix <- .check_radix(radix)
  S0 <- as.numeric(S0)

  if (any(!is.finite(S0))) {
    stop("S0 must contain only finite values.", call. = FALSE)
  }
  if (any(S0 < 0 | S0 > 1)) {
    stop("S0 values must lie in [0, 1].", call. = FALSE)
  }
  if (is.unsorted(-S0)) {
    stop("S0 must be nonincreasing.", call. = FALSE)
  }

  radix * S0
}

# -------------------------------------------------------------------------
# lx -> S0
# -------------------------------------------------------------------------

#' Convert life-table values to survival probabilities
#'
#' Converts Chapter 6 life-table values \eqn{l_x} into Chapter 5
#' survival probabilities \eqn{S_0(x) = l_x / l_0}.
#'
#' @param lx Numeric vector of life-table survivor values.
#'
#' @return Numeric vector of \eqn{S_0(x)} values.
#' @export
lx_to_S0 <- function(lx) {
  lx <- as.numeric(lx)

  if (length(lx) == 0L) {
    stop("lx must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(lx))) {
    stop("lx must contain only finite values.", call. = FALSE)
  }
  if (any(lx < 0)) {
    stop("lx values must be nonnegative.", call. = FALSE)
  }
  if (is.unsorted(-lx)) {
    stop("lx must be nonincreasing.", call. = FALSE)
  }
  if (lx[1] <= 0) {
    stop("The radix l0 = lx[1] must be positive.", call. = FALSE)
  }

  lx / lx[1]
}

# -------------------------------------------------------------------------
# px -> lx
# -------------------------------------------------------------------------

#' Construct life-table values from p_x values
#'
#' Builds life-table survivor values recursively from
#' \eqn{l_{x+1} = l_x p_x}, starting from a chosen radix.
#'
#' @param px Numeric vector of one-year survival probabilities \eqn{p_x}.
#' @param radix Positive radix \eqn{l_0}.
#'
#' @return Numeric vector of \eqn{l_x} values of length \code{length(px)+1}.
#' @export
px_to_lx <- function(px, radix = 100000) {
  radix <- .check_radix(radix)
  px <- as.numeric(px)

  if (length(px) == 0L) {
    stop("px must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(px))) {
    stop("px must contain only finite values.", call. = FALSE)
  }
  if (any(px < 0 | px > 1)) {
    stop("px values must lie in [0, 1].", call. = FALSE)
  }

  lx <- numeric(length(px) + 1L)
  lx[1] <- radix

  for (i in seq_along(px)) {
    lx[i + 1L] <- lx[i] * px[i]
  }

  lx
}

# -------------------------------------------------------------------------
# qx -> lx
# -------------------------------------------------------------------------

#' Construct life-table values from q_x values
#'
#' Builds life-table survivor values recursively from
#' \eqn{l_{x+1} = l_x (1-q_x)}, starting from a chosen radix.
#'
#' @param qx Numeric vector of one-year death probabilities \eqn{q_x}.
#' @param radix Positive radix \eqn{l_0}.
#'
#' @return Numeric vector of \eqn{l_x} values of length \code{length(qx)+1}.
#' @export
qx_to_lx <- function(qx, radix = 100000) {
  radix <- .check_radix(radix)
  qx <- as.numeric(qx)

  if (length(qx) == 0L) {
    stop("qx must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(qx))) {
    stop("qx must contain only finite values.", call. = FALSE)
  }
  if (any(qx < 0 | qx > 1)) {
    stop("qx values must lie in [0, 1].", call. = FALSE)
  }

  lx <- numeric(length(qx) + 1L)
  lx[1] <- radix

  for (i in seq_along(qx)) {
    lx[i + 1L] <- lx[i] * (1 - qx[i])
  }

  lx
}
