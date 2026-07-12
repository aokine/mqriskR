#' m-thly insurance models
#'
#' Exact m-thly contingent payment and insurance functions.
#'
#' @name insurance_mthly
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

.check_m_mthly <- function(m) {
  m <- as.numeric(m)

  if (length(m) != 1L || !is.finite(m) || m <= 0 ||
      abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }

  as.integer(round(m))
}

.check_i_mthly <- function(i) {
  i <- as.numeric(i)

  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }

  i
}

.recycle_mthly_vectors <- function(...) {
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

.check_x_mthly <- function(x) {
  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  invisible(TRUE)
}

.check_n_mthly <- function(n) {
  if (any(n < 0)) {
    stop("n must be nonnegative.", call. = FALSE)
  }

  invisible(TRUE)
}

# Kept for backward compatibility with other package files.
.recycle_x_n_m <- function(x, n) {
  args <- .recycle_mthly_vectors(x, n)
  x <- args[[1L]]
  n <- args[[2L]]

  .check_x_mthly(x)
  .check_n_mthly(n)

  list(x = x, n = n)
}

.get_upper_m <- function(x, n = NULL, model, ...) {
  if (tolower(model) == "uniform") {
    dots <- list(...)

    if (is.null(dots$omega)) {
      stop("For model = 'uniform', omega must be supplied.", call. = FALSE)
    }

    upper <- dots$omega - x

    if (!is.null(n)) {
      upper <- pmin(upper, n)
    }

    return(pmax(upper, 0))
  }

  if (is.null(n)) {
    return(Inf)
  }

  n
}

.pure_endowment_m <- function(x, n, i, model, ...) {
  discount(i, n) * tpx(n, x = x, model = model, ...)
}

.max_j_m <- function(x, m, model, ..., tol = 1e-12, j_max = 100000L) {
  upper <- .get_upper_m(x, model = model, ...)

  if (is.finite(upper)) {
    return(max(ceiling(upper * m) - 1L, 0L))
  }

  dt <- 1 / m

  for (j in 0:j_max) {
    s_j <- tpx(j * dt, x = x, model = model, ...)
    if (!is.finite(s_j) || s_j <= tol) {
      return(j)
    }
  }

  j_max
}

# ------------------------------------------------------------
# First moments
# ------------------------------------------------------------

#' m-thly whole life insurance APV
#'
#' Computes
#' \eqn{A_x^{(m)} =
#' \sum_{j=0}^{\infty} v^{(j+1)/m}
#' \Pr(j/m < T_x \le (j+1)/m)}.
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
  m <- .check_m_mthly(m)

  args <- .recycle_mthly_vectors(x, .check_i_mthly(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_x_mthly(x)

  dt <- 1 / m

  sapply(seq_along(x), function(k) {
    xx <- x[k]
    ii <- i[k]

    jj_max <- .max_j_m(xx, m, model = model, ..., tol = tol, j_max = j_max)
    j <- 0:jj_max

    p_start <- tpx(j * dt, x = xx, model = model, ...)
    p_end <- tpx((j + 1) * dt, x = xx, model = model, ...)
    probs <- p_start - p_end

    sum(discount(ii, (j + 1) * dt) * probs, na.rm = TRUE)
  })
}

#' m-thly term insurance APV
#'
#' Computes
#' \eqn{A_{x:\overline{n}|}^{1(m)} =
#' \sum_{j=0}^{mn-1} v^{(j+1)/m}
#' \Pr(j/m < T_x \le (j+1)/m)}.
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
  m <- .check_m_mthly(m)

  args <- .recycle_mthly_vectors(x, n, .check_i_mthly(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_mthly(x)
  .check_n_mthly(n)

  dt <- 1 / m

  sapply(seq_along(x), function(k) {
    xx <- x[k]
    nn <- n[k]
    ii <- i[k]

    upper <- .get_upper_m(xx, n = nn, model = model, ...)

    if (upper <= 0) {
      return(0)
    }

    jj_max <- max(ceiling(upper * m) - 1L, 0L)
    j <- 0:jj_max

    p_start <- tpx(j * dt, x = xx, model = model, ...)
    p_end <- tpx((j + 1) * dt, x = xx, model = model, ...)
    probs <- p_start - p_end

    sum(discount(ii, (j + 1) * dt) * probs, na.rm = TRUE)
  })
}

#' m-thly deferred insurance APV
#'
#' Computes
#' \eqn{{}_{n\mid}A_x^{(m)} = v^n {}_np_x A_{x+n}^{(m)}}.
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
  m <- .check_m_mthly(m)

  args <- .recycle_mthly_vectors(x, n, .check_i_mthly(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_mthly(x)
  .check_n_mthly(n)

  .pure_endowment_m(x, n, i, model = model, ...) *
    Ax_m(x + n, i, m, model = model, ..., tol = tol, j_max = j_max)
}

#' m-thly endowment insurance APV
#'
#' Computes
#' \eqn{A_{x:\overline{n}|}^{(m)} =
#' A_{x:\overline{n}|}^{1(m)} + v^n {}_np_x}.
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
  m <- .check_m_mthly(m)

  args <- .recycle_mthly_vectors(x, n, .check_i_mthly(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_mthly(x)
  .check_n_mthly(n)

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
  Ax_m(x, i = double_force_i(i), m = m, model = model, ...,
       tol = tol, j_max = j_max)
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
  nAx_m(x, n, i = double_force_i(i), m = m, model = model, ...,
        tol = tol, j_max = j_max)
}

#' Second moment of m-thly endowment insurance PV
#'
#' @inheritParams Axn_m
#' @return Numeric vector of second moments.
#' @export
A2xn_m <- function(x, n, i, m, model, ...) {
  m <- .check_m_mthly(m)

  args <- .recycle_mthly_vectors(x, n, .check_i_mthly(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  A2xn1_m(x, n, i, m, model = model, ...) +
    .pure_endowment_m(x, n, double_force_i(i), model = model, ...)
}

# ------------------------------------------------------------
# Variances
# ------------------------------------------------------------

#' Variance of m-thly whole life insurance PV
#'
#' Computes \eqn{\mathrm{Var}(Z_x^{(m)}) =
#' {}^{2}A_x^{(m)} - (A_x^{(m)})^2}.
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
