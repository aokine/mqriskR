#' Insurance models with varying benefits (Chapter 7)
#'
#' Functions for discrete and continuous contingent payment models with
#' varying benefits, matching Chapter 7 notation.
#'
#' Included functions:
#' \itemize{
#'   \item \eqn{(IA)_x}
#'   \item \eqn{(IA)_{x:\angl{n}}^{1}}
#'   \item \eqn{(DA)_{x:\angl{n}}^{1}}
#'   \item \eqn{(\bar{I}\bar{A})_x}
#'   \item \eqn{(I\bar{A})_x}
#'   \item \eqn{(\bar{I}\bar{A})_{x:\angl{n}}^{1}}
#'   \item \eqn{(\bar{D}\bar{A})_{x:\angl{n}}^{1}}
#'   \item \eqn{(D\bar{A})_{x:\angl{n}}^{1}}
#' }
#'
#' @name insurance_varying_benefits
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

#' @noRd
.discrete_prob_k <- function(x, k, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(npx(tbl, x, k) - npx(tbl, x, k + 1))
  }
  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  tpx(k, x = x, model = model, ...) - tpx(k + 1, x = x, model = model, ...)
}

#' @noRd
.max_k_varying <- function(x, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  if (!is.null(tbl)) {
    k <- 0
    repeat {
      if (k >= k_max) return(k_max)
      surv <- npx(tbl, x, k)
      if (surv <= tol) return(k)
      k <- k + 1
    }
  }

  for (k in 0:k_max) {
    surv <- tpx(k, x = x, model = model, ...)
    if (!is.finite(surv) || surv <= tol) return(k)
  }
  k_max
}

#' @noRd
.get_upper_varying_cont <- function(x, n = NULL, model, ...) {
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
.integrate_varying_cont <- function(f, lower, upper) {
  stats::integrate(
    f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-9,
    subdivisions = 1000L
  )$value
}

# ------------------------------------------------------------
# Discrete varying-benefit models
# ------------------------------------------------------------

#' Increasing whole life insurance
#'
#' Computes
#' \deqn{(IA)_x = \sum_{k=0}^{\infty} (k+1) v^{k+1} \Pr(K_x = k).}
#'
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional model parameters.
#' @param tol Numerical tolerance for truncation.
#' @param k_max Maximum number of terms.
#'
#' @return Numeric vector.
#' @export
IAx <- function(x, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  x <- as.numeric(x)
  if (any(!is.finite(x) | x < 0)) stop("x must be non-negative and finite.", call. = FALSE)
  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  v <- 1 / (1 + i)

  sapply(x, function(xx) {
    kk_max <- .max_k_varying(xx, tbl = tbl, model = model, ..., tol = tol, k_max = k_max)
    k <- 0:kk_max
    probs <- vapply(k, function(kk) .discrete_prob_k(xx, kk, tbl = tbl, model = model, ...), numeric(1))
    sum((k + 1) * v^(k + 1) * probs)
  })
}

#' Increasing n-year term insurance
#'
#' Computes
#' \deqn{(IA)_{x:\angl{n}}^{1} = \sum_{k=0}^{n-1} (k+1) v^{k+1} \Pr(K_x = k).}
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @export
IAxn1 <- function(x, n, i, tbl = NULL, model = NULL, ...) {
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
    probs <- vapply(k, function(kk) .discrete_prob_k(xx, kk, tbl = tbl, model = model, ...), numeric(1))
    sum((k + 1) * v^(k + 1) * probs)
  })
}

#' Decreasing n-year term insurance
#'
#' Computes
#' \deqn{(DA)_{x:\angl{n}}^{1} = \sum_{k=0}^{n-1} (n-k) v^{k+1} \Pr(K_x = k).}
#'
#' @inheritParams IAxn1
#'
#' @return Numeric vector.
#' @export
DAxn1 <- function(x, n, i, tbl = NULL, model = NULL, ...) {
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
    probs <- vapply(k, function(kk) .discrete_prob_k(xx, kk, tbl = tbl, model = model, ...), numeric(1))
    sum((nn - k) * v^(k + 1) * probs)
  })
}

# ------------------------------------------------------------
# Continuous varying-benefit models
# ------------------------------------------------------------

#' Fully continuous increasing whole life insurance
#'
#' Computes
#' \deqn{(\bar{I}\bar{A})_x = \int_0^\infty t\, v^t\, {}_tp_x\, \mu_{x+t}\, dt.}
#'
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param model Parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @export
IbarAbarx <- function(x, i, model, ...) {
  x <- as.numeric(x)
  if (any(!is.finite(x) | x < 0)) stop("x must be non-negative and finite.", call. = FALSE)
  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  delta <- interest_convert(i = i)$delta

  sapply(x, function(xx) {
    upper <- .get_upper_varying_cont(xx, model = model, ...)
    if (upper <= 0) return(0)

    f <- function(t) {
      val <- t * exp(-delta * t) *
        tpx(t, x = xx, model = model, ...) *
        hazard0(xx + t, model = model, ...)
      val[!is.finite(val)] <- 0
      val
    }

    .integrate_varying_cont(f, 0, upper)
  })
}

#' Piecewise-continuous increasing whole life insurance
#'
#' Computes
#' \deqn{(I\bar{A})_x = \int_0^\infty \lfloor t+1 \rfloor\, v^t\, {}_tp_x\, \mu_{x+t}\, dt.}
#'
#' @inheritParams IbarAbarx
#'
#' @return Numeric vector.
#' @export
IAbarx <- function(x, i, model, ...) {
  x <- as.numeric(x)
  if (any(!is.finite(x) | x < 0)) stop("x must be non-negative and finite.", call. = FALSE)
  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  delta <- interest_convert(i = i)$delta

  sapply(x, function(xx) {
    upper <- .get_upper_varying_cont(xx, model = model, ...)
    if (upper <= 0) return(0)

    f <- function(t) {
      val <- floor(t + 1) * exp(-delta * t) *
        tpx(t, x = xx, model = model, ...) *
        hazard0(xx + t, model = model, ...)
      val[!is.finite(val)] <- 0
      val
    }

    .integrate_varying_cont(f, 0, upper)
  })
}

#' Fully continuous increasing n-year term insurance
#'
#' Computes
#' \deqn{(\bar{I}\bar{A})_{x:\angl{n}}^{1} = \int_0^n t\, v^t\, {}_tp_x\, \mu_{x+t}\, dt.}
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param model Parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @export
IbarAbarxn1 <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_cont(x, n)
  x <- xn$x
  n <- xn$n

  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  delta <- interest_convert(i = i)$delta

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]

    upper <- .get_upper_varying_cont(xx, n = nn, model = model, ...)
    if (upper <= 0) return(0)

    f <- function(t) {
      val <- t * exp(-delta * t) *
        tpx(t, x = xx, model = model, ...) *
        hazard0(xx + t, model = model, ...)
      val[!is.finite(val)] <- 0
      val
    }

    .integrate_varying_cont(f, 0, upper)
  })
}

#' Fully continuous decreasing n-year term insurance
#'
#' Computes
#' \deqn{(\bar{D}\bar{A})_{x:\angl{n}}^{1} = \int_0^n (n-t)\, v^t\, {}_tp_x\, \mu_{x+t}\, dt.}
#'
#' @inheritParams IbarAbarxn1
#'
#' @return Numeric vector.
#' @export
DbarAbarxn1 <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_cont(x, n)
  x <- xn$x
  n <- xn$n

  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  delta <- interest_convert(i = i)$delta

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]

    upper <- .get_upper_varying_cont(xx, n = nn, model = model, ...)
    if (upper <= 0) return(0)

    f <- function(t) {
      val <- (nn - t) * exp(-delta * t) *
        tpx(t, x = xx, model = model, ...) *
        hazard0(xx + t, model = model, ...)
      val[!is.finite(val)] <- 0
      val
    }

    .integrate_varying_cont(f, 0, upper)
  })
}

#' Piecewise-continuous decreasing n-year term insurance
#'
#' Computes
#' \deqn{(D\bar{A})_{x:\angl{n}}^{1} = \int_0^n \lfloor n+1-t \rfloor\, v^t\, {}_tp_x\, \mu_{x+t}\, dt.}
#'
#' @inheritParams IbarAbarxn1
#'
#' @return Numeric vector.
#' @export
DAbarxn1 <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_cont(x, n)
  x <- xn$x
  n <- xn$n

  if (length(i) != 1L || !is.finite(i) || i <= -1) stop("i must be a single value > -1.", call. = FALSE)

  delta <- interest_convert(i = i)$delta

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]

    upper <- .get_upper_varying_cont(xx, n = nn, model = model, ...)
    if (upper <= 0) return(0)

    f <- function(t) {
      val <- floor(nn + 1 - t) * exp(-delta * t) *
        tpx(t, x = xx, model = model, ...) *
        hazard0(xx + t, model = model, ...)
      val[!is.finite(val)] <- 0
      val
    }

    .integrate_varying_cont(f, 0, upper)
  })
}
