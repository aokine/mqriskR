#' Discrete insurance models (Chapter 7)
#'
#' Discrete contingent payment / insurance functions matching Chapter 7 notation.
#'
#' These functions handle:
#' \itemize{
#'   \item whole life insurance: \eqn{A_x},
#'   \item term insurance: \eqn{A_{x:\angl{n}}^{1}},
#'   \item deferred insurance: \eqn{{}_{n\mid}A_x},
#'   \item pure endowment: \eqn{{}_nE_x},
#'   \item endowment insurance: \eqn{A_{x:\angl{n}}},
#'   \item second moments and variances.
#' }
#'
#' The functions may be evaluated from either:
#' \itemize{
#'   \item a life table object via \code{tbl = ...}, or
#'   \item a parametric survival model via \code{model = ...} and additional parameters.
#' }
#'
#' @name insurance_discrete
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

#' @noRd
.get_px1_discrete <- function(x, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(npx(tbl, x = x, n = 1))
  }
  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  tpx(1, x = x, model = model, ...)
}

#' @noRd
.get_qx1_discrete <- function(x, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(nqx(tbl, x = x, n = 1))
  }
  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  tqx(1, x = x, model = model, ...)
}

#' @noRd
.get_npx_discrete <- function(x, n, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(npx(tbl, x = x, n = n))
  }
  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  tpx(n, x = x, model = model, ...)
}

#' @noRd
.get_nqx_discrete <- function(x, n, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(nqx(tbl, x = x, n = n))
  }
  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  tqx(n, x = x, model = model, ...)
}

#' @noRd
.recycle_x_n <- function(x, n) {
  x <- as.numeric(x)
  n <- as.numeric(n)

  if (length(x) == 0L || length(n) == 0L) {
    stop("x and n must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(!is.finite(n))) stop("n must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be >= 0.", call. = FALSE)
  if (any(n < 0)) stop("n must be >= 0.", call. = FALSE)
  if (any(abs(n - round(n)) > 1e-10)) stop("n must contain integers only.", call. = FALSE)

  if (length(x) == 1L && length(n) > 1L) x <- rep(x, length(n))
  if (length(n) == 1L && length(x) > 1L) n <- rep(n, length(x))
  if (length(x) != length(n)) {
    stop("x and n must have compatible lengths.", call. = FALSE)
  }

  list(x = x, n = n)
}

#' @noRd
.max_k_discrete <- function(x, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  if (!is.null(tbl)) {
    k <- 0
    repeat {
      if (k >= k_max) return(k_max)
      surv <- .get_npx_discrete(x, k, tbl = tbl)
      if (surv <= tol) return(k)
      k <- k + 1
    }
  }

  for (k in 0:k_max) {
    surv <- .get_npx_discrete(x, k, model = model, ...)
    if (surv <= tol) return(k)
  }
  k_max
}

#' @noRd
.kq_vec <- function(x, k, tbl = NULL, model = NULL, ...) {
  # {}_{k|}q_x = {}_k p_x - {}_{k+1} p_x
  .get_npx_discrete(x, k, tbl = tbl, model = model, ...) -
    .get_npx_discrete(x, k + 1, tbl = tbl, model = model, ...)
}

# ------------------------------------------------------------
# First moments
# ------------------------------------------------------------

#' Whole life insurance APV
#'
#' Computes \eqn{A_x = \sum_{k=0}^\infty v^{k+1} {}_{k\mid}q_x}.
#'
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#' @param tol Numerical tolerance for truncating infinite sums.
#' @param k_max Maximum number of terms in the sum.
#'
#' @return Numeric vector of APVs.
#' @export
Ax <- function(x, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  x <- as.numeric(x)
  i <- as.numeric(i)

  if (any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be >= 0.", call. = FALSE)
  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  v <- 1 / (1 + i)

  sapply(x, function(xx) {
    kk_max <- .max_k_discrete(xx, tbl = tbl, model = model, ..., tol = tol, k_max = k_max)
    k <- 0:kk_max
    probs <- vapply(k, function(kk) .kq_vec(xx, kk, tbl = tbl, model = model, ...), numeric(1))
    sum(v^(k + 1) * probs)
  })
}

#' Term insurance APV
#'
#' Computes \eqn{A_{x:\angl{n}}^{1} = \sum_{k=0}^{n-1} v^{k+1} {}_{k\mid}q_x}.
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Axn1 <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  xn <- .recycle_x_n(x, n)
  x <- xn$x
  n <- xn$n

  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  v <- 1 / (1 + i)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]
    if (nn == 0) return(0)

    k <- 0:(nn - 1)
    probs <- vapply(k, function(kk) .kq_vec(xx, kk, tbl = tbl, model = model, ...), numeric(1))
    sum(v^(k + 1) * probs)
  })
}

#' Pure endowment APV
#'
#' Computes \eqn{{}_nE_x = v^n \, {}_n p_x}.
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
nEx <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  xn <- .recycle_x_n(x, n)
  x <- xn$x
  n <- xn$n

  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  discount(i, n) * .get_npx_discrete(x, n, tbl = tbl, model = model, ...)
}

#' Deferred insurance APV
#'
#' Computes \eqn{{}_{n\mid}A_x = {}_nE_x \cdot A_{x+n}}.
#'
#' @param x Age.
#' @param n Deferral period.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#' @param tol Numerical tolerance for truncating infinite sums.
#' @param k_max Maximum number of terms in the sum.
#'
#' @return Numeric vector of APVs.
#' @export
nAx <- function(x, n, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  xn <- .recycle_x_n(x, n)
  x <- xn$x
  n <- xn$n

  nEx(x, n, i, tbl = tbl, model = model, ...) *
    Ax(x + n, i, tbl = tbl, model = model, ..., tol = tol, k_max = k_max)
}

#' Endowment insurance APV
#'
#' Computes \eqn{A_{x:\angl{n}} = A_{x:\angl{n}}^{1} + {}_nE_x}.
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Axn <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  Axn1(x, n, i, tbl = tbl, model = model, ...) +
    nEx(x, n, i, tbl = tbl, model = model, ...)
}

# ------------------------------------------------------------
# Second moments
# ------------------------------------------------------------

#' Second moment of whole life insurance PV
#'
#' Computes \eqn{{}^{2}A_x} by evaluating \eqn{A_x} at doubled force.
#'
#' @inheritParams Ax
#' @return Numeric vector of second moments.
#' @export
A2x <- function(x, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  Ax(x, i = double_force_i(i), tbl = tbl, model = model, ..., tol = tol, k_max = k_max)
}

#' Second moment of term insurance PV
#'
#' Computes \eqn{{}^{2}A_{x:\angl{n}}^{1}} by evaluating
#' \eqn{A_{x:\angl{n}}^{1}} at doubled force.
#'
#' @inheritParams Axn1
#' @return Numeric vector of second moments.
#' @export
A2xn1 <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  Axn1(x, n, i = double_force_i(i), tbl = tbl, model = model, ...)
}

#' Second moment of pure endowment PV
#'
#' Computes \eqn{{}^{2}{}_nE_x = (v')^n {}_n p_x}.
#'
#' @inheritParams nEx
#' @return Numeric vector of second moments.
#' @export
A2nEx <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  nEx(x, n, i = double_force_i(i), tbl = tbl, model = model, ...)
}

#' Second moment of deferred insurance PV
#'
#' Computes \eqn{{}^{2}{}_{n\mid}A_x} by evaluating
#' \eqn{{}_{n\mid}A_x} at doubled force.
#'
#' @inheritParams nAx
#' @return Numeric vector of second moments.
#' @export
A2nAx <- function(x, n, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  nAx(x, n, i = double_force_i(i), tbl = tbl, model = model, ..., tol = tol, k_max = k_max)
}

#' Second moment of endowment insurance PV
#'
#' Computes \eqn{{}^{2}A_{x:\angl{n}}}.
#'
#' @inheritParams Axn
#' @return Numeric vector of second moments.
#' @export
A2xn <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  A2xn1(x, n, i, tbl = tbl, model = model, ...) +
    A2nEx(x, n, i, tbl = tbl, model = model, ...)
}

# ------------------------------------------------------------
# Variances
# ------------------------------------------------------------

#' Variance of whole life insurance PV
#'
#' Computes \eqn{\mathrm{Var}(Z_x) = {}^{2}A_x - A_x^2}.
#'
#' @inheritParams Ax
#' @return Numeric vector of variances.
#' @export
var_Ax <- function(x, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  A2x(x, i, tbl = tbl, model = model, ..., tol = tol, k_max = k_max) -
    Ax(x, i, tbl = tbl, model = model, ..., tol = tol, k_max = k_max)^2
}

#' Variance of term insurance PV
#'
#' Computes \eqn{\mathrm{Var}(Z_{x:\angl{n}}^{1}) = {}^{2}A_{x:\angl{n}}^{1} - (A_{x:\angl{n}}^{1})^2}.
#'
#' @inheritParams Axn1
#' @return Numeric vector of variances.
#' @export
var_Axn1 <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  A2xn1(x, n, i, tbl = tbl, model = model, ...) -
    Axn1(x, n, i, tbl = tbl, model = model, ...)^2
}

#' Variance of pure endowment PV
#'
#' @inheritParams nEx
#' @return Numeric vector of variances.
#' @export
var_nEx <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  A2nEx(x, n, i, tbl = tbl, model = model, ...) -
    nEx(x, n, i, tbl = tbl, model = model, ...)^2
}

#' Variance of deferred insurance PV
#'
#' @inheritParams nAx
#' @return Numeric vector of variances.
#' @export
var_nAx <- function(x, n, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  A2nAx(x, n, i, tbl = tbl, model = model, ..., tol = tol, k_max = k_max) -
    nAx(x, n, i, tbl = tbl, model = model, ..., tol = tol, k_max = k_max)^2
}

#' Variance of endowment insurance PV
#'
#' @inheritParams Axn
#' @return Numeric vector of variances.
#' @export
var_Axn <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  A2xn(x, n, i, tbl = tbl, model = model, ...) -
    Axn(x, n, i, tbl = tbl, model = model, ...)^2
}

# ------------------------------------------------------------
# Covariances from Chapter 7 decompositions
# ------------------------------------------------------------

#' Covariance of term and deferred insurance PVs
#'
#' Computes
#' \eqn{\mathrm{Cov}(Z_{x:\angl{n}}^{1}, {}_{n\mid}Z_x) = -A_{x:\angl{n}}^{1} \cdot {}_{n\mid}A_x}.
#'
#' @inheritParams Axn1
#' @return Numeric vector of covariances.
#' @export
cov_term_deferred <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  -Axn1(x, n, i, tbl = tbl, model = model, ...) *
    nAx(x, n, i, tbl = tbl, model = model, ...)
}

#' Covariance of term insurance and pure endowment PVs
#'
#' Computes
#' \eqn{\mathrm{Cov}(Z_{x:\angl{n}}^{1}, Z_{x:\angl{n}}^{1\text{(pure endow)}}) =
#' -A_{x:\angl{n}}^{1} \cdot {}_nE_x}.
#'
#' @inheritParams Axn1
#' @return Numeric vector of covariances.
#' @export
cov_term_endow <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  -Axn1(x, n, i, tbl = tbl, model = model, ...) *
    nEx(x, n, i, tbl = tbl, model = model, ...)
}
