#' Discrete insurance models
#'
#' Discrete contingent payment and insurance functions.
#'
#' @name insurance_discrete
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

.get_npx_discrete <- function(x, n, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) return(npx(tbl, x = x, n = n))
  if (is.null(model)) stop("Supply either tbl or model.", call. = FALSE)
  tpx(n, x = x, model = model, ...)
}

.get_nqx_discrete <- function(x, n, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) return(nqx(tbl, x = x, n = n))
  if (is.null(model)) stop("Supply either tbl or model.", call. = FALSE)
  tqx(n, x = x, model = model, ...)
}

.recycle_x_n <- function(x, n) {
  x <- as.numeric(x)
  n <- as.numeric(n)

  if (length(x) == 0L || length(n) == 0L) {
    stop("x and n must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(!is.finite(n))) stop("n must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be nonnegative.", call. = FALSE)
  if (any(n < 0)) stop("n must be nonnegative.", call. = FALSE)
  if (any(abs(n - round(n)) > 1e-10)) {
    stop("n must contain integers only.", call. = FALSE)
  }

  len <- max(length(x), length(n))

  if (!length(x) %in% c(1L, len) || !length(n) %in% c(1L, len)) {
    stop("x and n must have compatible lengths.", call. = FALSE)
  }

  list(
    x = if (length(x) == 1L) rep(x, len) else x,
    n = if (length(n) == 1L) rep(n, len) else n
  )
}

.check_i_discrete <- function(i) {
  i <- as.numeric(i)

  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value greater than -1.", call. = FALSE)
  }

  i
}

.max_k_discrete <- function(x, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  for (k in 0:k_max) {
    surv <- .get_npx_discrete(x, k, tbl = tbl, model = model, ...)

    if (is.na(surv)) return(max(k - 1L, 0L))
    if (surv <= tol) return(k)
  }

  k_max
}

.kq_vec <- function(x, k, tbl = NULL, model = NULL, ...) {
  p_k <- .get_npx_discrete(x, k, tbl = tbl, model = model, ...)
  p_k1 <- .get_npx_discrete(x, k + 1, tbl = tbl, model = model, ...)

  if (is.na(p_k)) p_k <- 0
  if (is.na(p_k1)) p_k1 <- 0

  p_k - p_k1
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
  i <- .check_i_discrete(i)

  if (length(x) == 0L || any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be nonnegative.", call. = FALSE)

  v <- 1 / (1 + i)

  sapply(x, function(xx) {
    kk_max <- .max_k_discrete(xx, tbl = tbl, model = model, ..., tol = tol, k_max = k_max)
    k <- 0:kk_max
    probs <- vapply(k, function(kk) .kq_vec(xx, kk, tbl = tbl, model = model, ...), numeric(1))
    sum(v^(k + 1) * probs, na.rm = TRUE)
  })
}

#' Term insurance APV
#'
#' Computes \eqn{A_{x:\overline{n}|}^{1} = \sum_{k=0}^{n-1} v^{k+1} {}_{k\mid}q_x}.
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
  i <- .check_i_discrete(i)

  v <- 1 / (1 + i)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]

    if (nn == 0) return(0)

    k <- 0:(nn - 1)
    probs <- vapply(k, function(kk) .kq_vec(xx, kk, tbl = tbl, model = model, ...), numeric(1))

    sum(v^(k + 1) * probs, na.rm = TRUE)
  })
}

#' Pure endowment APV
#'
#' Computes \eqn{{}_nE_x = v^n {}_n p_x}.
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
  i <- .check_i_discrete(i)

  p <- .get_npx_discrete(x, n, tbl = tbl, model = model, ...)
  p[is.na(p)] <- 0

  discount(i, n) * p
}

#' Deferred insurance APV
#'
#' Computes \eqn{{}_{n\mid}A_x = {}_nE_x A_{x+n}}.
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
#' Computes \eqn{A_{x:\overline{n}|} = A_{x:\overline{n}|}^{1} + {}_nE_x}.
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
#' Computes \eqn{{}^{2}A_{x:\overline{n}|}^{1}} by evaluating
#' \eqn{A_{x:\overline{n}|}^{1}} at doubled force.
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
#' Computes \eqn{{}^{2}A_{x:\overline{n}|}}.
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
#' Computes \eqn{\mathrm{Var}(Z_{x:\overline{n}|}^{1}) =
#' {}^{2}A_{x:\overline{n}|}^{1} - (A_{x:\overline{n}|}^{1})^2}.
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
# Covariances
# ------------------------------------------------------------

#' Covariance of term and deferred insurance PVs
#'
#' Computes
#' \eqn{\mathrm{Cov}(Z_{x:\overline{n}|}^{1}, {}_{n\mid}Z_x) =
#' -A_{x:\overline{n}|}^{1} \cdot {}_{n\mid}A_x}.
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
#' \eqn{\mathrm{Cov}(Z_{x:\overline{n}|}^{1},
#' Z_{x:\overline{n}|}^{1\text{(pure endow)}}) =
#' -A_{x:\overline{n}|}^{1} \cdot {}_nE_x}.
#'
#' @inheritParams Axn1
#' @return Numeric vector of covariances.
#' @export
cov_term_endow <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  -Axn1(x, n, i, tbl = tbl, model = model, ...) *
    nEx(x, n, i, tbl = tbl, model = model, ...)
}
