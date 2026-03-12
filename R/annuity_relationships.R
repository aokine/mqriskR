#' Annuity-insurance relationships (Chapter 8)
#'
#' This file provides the core Chapter 8 identities linking annual and
#' continuous annuity functions to the corresponding insurance functions.
#'
#' Included identities:
#' \itemize{
#'   \item whole life immediate: \eqn{a_x = (v - A_x)/d}
#'   \item whole life due: \eqn{\ddot{a}_x = (1 - A_x)/d}
#'   \item whole life continuous: \eqn{\bar{a}_x = (1 - \bar{A}_x)/\delta}
#'   \item temporary immediate: \eqn{a_{x:\angl{n}} = (1 - A_{x:\angl{n}})/d - 1 + {}_nE_x}
#'   \item temporary due: \eqn{\ddot{a}_{x:\angl{n}} = (1 - A_{x:\angl{n}})/d}
#'   \item temporary continuous: \eqn{\bar{a}_{x:\angl{n}} = (1 - \bar{A}_{x:\angl{n}})/\delta}
#'   \item deferred immediate: \eqn{{}_{n\mid}a_x = {}_nE_x a_{x+n}}
#'   \item deferred due: \eqn{{}_{n\mid}\ddot{a}_x = {}_nE_x \ddot{a}_{x+n}}
#'   \item deferred continuous: \eqn{{}_{n\mid}\bar{a}_x = {}_nE_x \bar{a}_{x+n}}
#' }
#'
#' These are wrapper functions that evaluate the Chapter 8 relationships
#' using the Chapter 7 insurance functions already implemented in the package.
#'
#' @param x Age.
#' @param n Term in years.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters passed to the survival model.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return Numeric vector.
#' @name annuity_relationships
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_rel_i <- function(i) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite number greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.check_rel_n <- function(n) {
  n <- as.numeric(n)
  if (length(n) != 1L || !is.finite(n) || n < 0) {
    stop("n must be a single nonnegative number.", call. = FALSE)
  }
  n
}

#' @noRd
.recycle_rel_xn <- function(x, n) {
  x <- as.numeric(x)
  n <- as.numeric(n)

  if (length(x) == 0L || length(n) == 0L) {
    stop("x and n must have positive length.", call. = FALSE)
  }

  if (length(x) == 1L && length(n) > 1L) x <- rep(x, length(n))
  if (length(n) == 1L && length(x) > 1L) n <- rep(n, length(x))

  if (length(x) != length(n)) {
    stop("x and n must have the same length, or one must have length 1.", call. = FALSE)
  }

  list(x = x, n = n)
}

# -------------------------------------------------------------------------
# Whole life identities
# -------------------------------------------------------------------------

#' Whole life annuity-immediate from insurance identity
#'
#' Computes \eqn{a_x = (v - A_x)/d}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_ax <- function(x, i, model, ...) {
  .check_rel_i(i)
  d <- i / (1 + i)
  v <- 1 / (1 + i)

  (v - Ax(x, i, model = model, ...)) / d
}

#' Whole life annuity-due from insurance identity
#'
#' Computes \eqn{\ddot{a}_x = (1 - A_x)/d}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_adotx <- function(x, i, model, ...) {
  .check_rel_i(i)
  d <- i / (1 + i)

  (1 - Ax(x, i, model = model, ...)) / d
}

#' Continuous whole life annuity from insurance identity
#'
#' Computes \eqn{\bar{a}_x = (1 - \bar{A}_x)/\delta}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_abarx <- function(x, i, model, ...) {
  .check_rel_i(i)
  delta <- log(1 + i)

  (1 - Abarx(x, i, model = model, ...)) / delta
}

# -------------------------------------------------------------------------
# Temporary identities
# -------------------------------------------------------------------------

#' Temporary annuity-immediate from insurance identity
#'
#' Computes
#' \eqn{a_{x:\overline{n}|} = \ddot{a}_{x:\overline{n}|} - 1 + {}_nE_x}
#' together with
#' \eqn{\ddot{a}_{x:\overline{n}|} = (1 - A_{x:\overline{n}|})/d}.
#'
#' Hence
#' \eqn{a_{x:\overline{n}|} = (1 - A_{x:\overline{n}|})/d - 1 + {}_nE_x}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_axn <- function(x, n, i, model, ...) {
  .check_rel_i(i)
  d <- i / (1 + i)
  xn <- .recycle_rel_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity relationships, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    (1 - Axn(x[j], nn, i, model = model, ...)) / d - 1 +
      nEx(x[j], nn, i, model = model, ...)
  })
}

#' Temporary annuity-due from insurance identity
#'
#' Computes \eqn{\ddot{a}_{x:\overline{n}|} = (1 - A_{x:\overline{n}|})/d}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_adotxn <- function(x, n, i, model, ...) {
  .check_rel_i(i)
  d <- i / (1 + i)
  xn <- .recycle_rel_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity relationships, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    (1 - Axn(x[j], nn, i, model = model, ...)) / d
  })
}

#' Continuous temporary annuity from insurance identity
#'
#' Computes \eqn{\bar{a}_{x:\overline{n}|} = (1 - \bar{A}_{x:\overline{n}|})/\delta}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_abarxn <- function(x, n, i, model, ...) {
  .check_rel_i(i)
  delta <- log(1 + i)
  xn <- .recycle_rel_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])
    (1 - Abarxn(x[j], nn, i, model = model, ...)) / delta
  })
}

# -------------------------------------------------------------------------
# Deferred identities
# -------------------------------------------------------------------------

#' Deferred annuity-immediate from pure endowment identity
#'
#' Computes \eqn{{}_{n|}a_x = {}_nE_x\, a_{x+n}}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_nax <- function(x, n, i, model, ..., k_max = 5000, tol = 1e-12) {
  .check_rel_i(i)
  xn <- .recycle_rel_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity relationships, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    nEx(x[j], nn, i, model = model, ...) *
      ax(x[j] + nn, i, model = model, ..., k_max = k_max, tol = tol)
  })
}

#' Deferred annuity-due from pure endowment identity
#'
#' Computes \eqn{{}_{n|}\ddot{a}_x = {}_nE_x\, \ddot{a}_{x+n}}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_nadotx <- function(x, n, i, model, ..., k_max = 5000, tol = 1e-12) {
  .check_rel_i(i)
  xn <- .recycle_rel_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity relationships, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    nEx(x[j], nn, i, model = model, ...) *
      adotx(x[j] + nn, i, model = model, ..., k_max = k_max, tol = tol)
  })
}

#' Deferred continuous annuity from pure endowment identity
#'
#' Computes \eqn{{}_{n|}\bar{a}_x = {}_nE_x\, \bar{a}_{x+n}}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_nabarx <- function(x, n, i, model, ..., tol = 1e-10) {
  .check_rel_i(i)
  xn <- .recycle_rel_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])
    nEx(x[j], nn, i, model = model, ...) *
      abarx(x[j] + nn, i, model = model, ..., tol = tol)
  })
}
