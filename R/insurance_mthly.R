#' m-thly insurance models (Chapter 7)
#'
#' Exact m-thly contingent payment / insurance functions matching Chapter 7 notation.
#'
#' These functions handle:
#' \itemize{
#'   \item m-thly whole life insurance: \eqn{A_x^{(m)}},
#'   \item m-thly term insurance: \eqn{A_{x:\angl{n}}^{1(m)}},
#'   \item m-thly deferred insurance: \eqn{{}_{n\mid}A_x^{(m)}},
#'   \item m-thly endowment insurance: \eqn{A_{x:\angl{n}}^{(m)}},
#'   \item second moments and variances.
#' }
#'
#' These functions are evaluated exactly from the parametric survival-model
#' framework. For UDD approximations from annual life tables, use the
#' corresponding *_udd helpers in insurance_utils.R.
#'
#' @name insurance_mthly
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

#' @noRd
.check_mthly_inputs <- function(x, i, m) {
  x <- as.numeric(x)
  i <- as.numeric(i)

  if (any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be >= 0.", call. = FALSE)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value > -1.", call. = FALSE)
  }
  if (length(m) != 1L || !is.finite(m) || m <= 0 || m != as.integer(m)) {
    stop("m must be a positive integer.", call. = FALSE)
  }

  list(x = x, i = i, m = m)
}

#' @noRd
.recycle_x_n_m <- function(x, n) {
  x <- as.numeric(x)
  n <- as.numeric(n)

  if (length(x) == 0L || length(n) == 0L) {
    stop("x and n must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(!is.finite(n))) stop("n must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be >= 0.", call. = FALSE)
  if (any(n < 0)) stop("n must be >= 0.", call. = FALSE)

  if (length(x) == 1L && length(n) > 1L) x <- rep(x, length(n))
  if (length(n) == 1L && length(x) > 1L) n <- rep(n, length(x))
  if (length(x) != length(n)) {
    stop("x and n must have compatible lengths.", call. = FALSE)
  }

  list(x = x, n = n)
}

#' @noRd
.get_upper_m <- function(x, n = NULL, model, ...) {
  if (tolower(model) == "uniform") {
    dots <- list(...)
    if (is.null(dots$omega)) {
      stop("For model = 'uniform', omega must be supplied.", call. = FALSE)
    }
    upper <- dots$omega - x
    if (!is.null(n)) upper <- pmin(upper, n)
    return(pmax(upper, 0))
  }

  if (is.null(n)) return(Inf)
  n
}

#' @noRd
.pure_endowment_m <- function(x, n, i, model, ...) {
  discount(i, n) * tpx(n, x = x, model = model, ...)
}

#' @noRd
.max_j_m <- function(x, m, model, ..., tol = 1e-12, j_max = 100000L) {
  upper <- .get_upper_m(x, model = model, ...)
  if (is.finite(upper)) {
    return(max(ceiling(upper * m) - 1L, 0L))
  }

  dt <- 1 / m
  for (j in 0:j_max) {
    s_j <- tpx(j * dt, x = x, model = model, ...)
    if (!is.finite(s_j) || s_j <= tol) return(j)
  }
  j_max
}

# ------------------------------------------------------------
# First moments
# ------------------------------------------------------------

#' m-thly whole life insurance APV
#'
#' Computes
#' \eqn{A_x^{(m)} = \sum_{j=0}^\infty v^{(j+1)/m}\Pr(j/m < T_x \le (j+1)/m)}.
#'
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#' @param tol Numerical tolerance for truncating the infinite sum.
#' @param j_max Maximum number of m-thly intervals in the sum.
#'
#' @return Numeric vector of APVs.
#' @export
Ax_m <- function(x, i, m, model, ..., tol = 1e-12, j_max = 100000L) {
  chk <- .check_mthly_inputs(x, i, m)
  x <- chk$x

  dt <- 1 / m

  sapply(x, function(xx) {
    jj_max <- .max_j_m(xx, m, model = model, ..., tol = tol, j_max = j_max)
    j <- 0:jj_max

    p_start <- tpx(j * dt, x = xx, model = model, ...)
    p_end   <- tpx((j + 1) * dt, x = xx, model = model, ...)
    probs   <- p_start - p_end

    sum(discount(i, (j + 1) * dt) * probs)
  })
}

#' m-thly term insurance APV
#'
#' Computes
#' \eqn{A_{x:\angl{n}}^{1(m)} = \sum_{j=0}^{mn-1} v^{(j+1)/m}\Pr(j/m < T_x \le (j+1)/m)}.
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Axn1_m <- function(x, n, i, m, model, ...) {
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  chk <- .check_mthly_inputs(x[1], i, m)
  m <- chk$m
  dt <- 1 / m

  sapply(seq_along(x), function(k) {
    xx <- x[k]
    nn <- n[k]

    upper <- .get_upper_m(xx, n = nn, model = model, ...)
    if (upper <= 0) return(0)

    jj_max <- max(ceiling(upper * m) - 1L, 0L)
    j <- 0:jj_max

    p_start <- tpx(j * dt, x = xx, model = model, ...)
    p_end   <- tpx((j + 1) * dt, x = xx, model = model, ...)
    probs   <- p_start - p_end

    sum(discount(i, (j + 1) * dt) * probs)
  })
}

#' m-thly deferred insurance APV
#'
#' Computes
#' \eqn{{}_{n\mid}A_x^{(m)} = v^n\,{}_np_x\,A_{x+n}^{(m)}}.
#'
#' @param x Age.
#' @param n Deferral period.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#' @param tol Numerical tolerance for truncating the infinite sum.
#' @param j_max Maximum number of m-thly intervals in the sum.
#'
#' @return Numeric vector of APVs.
#' @export
nAx_m <- function(x, n, i, m, model, ..., tol = 1e-12, j_max = 100000L) {
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  .pure_endowment_m(x, n, i, model = model, ...) *
    Ax_m(x + n, i, m, model = model, ..., tol = tol, j_max = j_max)
}

#' m-thly endowment insurance APV
#'
#' Computes
#' \eqn{A_{x:\angl{n}}^{(m)} = A_{x:\angl{n}}^{1(m)} + v^n\,{}_np_x}.
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Axn_m <- function(x, n, i, m, model, ...) {
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  Axn1_m(x, n, i, m, model = model, ...) +
    .pure_endowment_m(x, n, i, model = model, ...)
}

# ------------------------------------------------------------
# Second moments
# ------------------------------------------------------------

#' Second moment of m-thly whole life insurance PV
#'
#' Computes \eqn{{}^{2}A_x^{(m)}} by evaluating \eqn{A_x^{(m)}}
#' at doubled force.
#'
#' @inheritParams Ax_m
#' @return Numeric vector of second moments.
#' @export
A2x_m <- function(x, i, m, model, ..., tol = 1e-12, j_max = 100000L) {
  Ax_m(x, i = double_force_i(i), m = m, model = model, ..., tol = tol, j_max = j_max)
}

#' Second moment of m-thly term insurance PV
#'
#' @inheritParams Axn1_m
#' @return Numeric vector of second moments.
#' @export
A2xn1_m <- function(x, n, i, m, model, ...) {
  Axn1_m(x, n, i = double_force_i(i), m = m, model = model, ...)
}

#' Second moment of m-thly deferred insurance PV
#'
#' @inheritParams nAx_m
#' @return Numeric vector of second moments.
#' @export
A2nAx_m <- function(x, n, i, m, model, ..., tol = 1e-12, j_max = 100000L) {
  nAx_m(x, n, i = double_force_i(i), m = m, model = model, ..., tol = tol, j_max = j_max)
}

#' Second moment of m-thly endowment insurance PV
#'
#' @inheritParams Axn_m
#' @return Numeric vector of second moments.
#' @export
A2xn_m <- function(x, n, i, m, model, ...) {
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  A2xn1_m(x, n, i, m, model = model, ...) +
    .pure_endowment_m(x, n, double_force_i(i), model = model, ...)
}

# ------------------------------------------------------------
# Variances
# ------------------------------------------------------------

#' Variance of m-thly whole life insurance PV
#'
#' Computes \eqn{\mathrm{Var}(Z_x^{(m)}) = {}^{2}A_x^{(m)} - (A_x^{(m)})^2}.
#'
#' @inheritParams Ax_m
#' @return Numeric vector of variances.
#' @export
var_Ax_m <- function(x, i, m, model, ..., tol = 1e-12, j_max = 100000L) {
  A2x_m(x, i, m, model = model, ..., tol = tol, j_max = j_max) -
    Ax_m(x, i, m, model = model, ..., tol = tol, j_max = j_max)^2
}

#' Variance of m-thly term insurance PV
#'
#' @inheritParams Axn1_m
#' @return Numeric vector of variances.
#' @export
var_Axn1_m <- function(x, n, i, m, model, ...) {
  A2xn1_m(x, n, i, m, model = model, ...) -
    Axn1_m(x, n, i, m, model = model, ...)^2
}

#' Variance of m-thly deferred insurance PV
#'
#' @inheritParams nAx_m
#' @return Numeric vector of variances.
#' @export
var_nAx_m <- function(x, n, i, m, model, ..., tol = 1e-12, j_max = 100000L) {
  A2nAx_m(x, n, i, m, model = model, ..., tol = tol, j_max = j_max) -
    nAx_m(x, n, i, m, model = model, ..., tol = tol, j_max = j_max)^2
}

#' Variance of m-thly endowment insurance PV
#'
#' @inheritParams Axn_m
#' @return Numeric vector of variances.
#' @export
var_Axn_m <- function(x, n, i, m, model, ...) {
  A2xn_m(x, n, i, m, model = model, ...) -
    Axn_m(x, n, i, m, model = model, ...)^2
}
