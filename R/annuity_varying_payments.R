#' Varying-payment annuity functions (Chapter 8)
#'
#' Chapter 8 non-level annuity functions for increasing and decreasing
#' life annuities.
#'
#' The functions implemented here match the notation in Section 8.6:
#' \itemize{
#'   \item \code{Iax()} = \eqn{(Ia)_x}
#'   \item \code{Iaxn()} = \eqn{(Ia)_{x:\angl{n}}}
#'   \item \code{Daxn()} = \eqn{(Da)_{x:\angl{n}}}
#'   \item \code{Iadotx()} = \eqn{(I\ddot{a})_x}
#'   \item \code{Iadotxn()} = \eqn{(I\ddot{a})_{x:\angl{n}}}
#'   \item \code{Dadotxn()} = \eqn{(D\ddot{a})_{x:\angl{n}}}
#'   \item \code{Iabarx()} = \eqn{(\bar{I}\bar{a})_x}
#'   \item \code{Iabarxn()} = \eqn{(\bar{I}\bar{a})_{x:\angl{n}}}
#'   \item \code{Dabarxn()} = \eqn{(\bar{D}\bar{a})_{x:\angl{n}}}
#' }
#'
#' @name annuity_varying_payments
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_annuity_vp_i <- function(i) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite number greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.check_annuity_vp_x <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stop("x must be a single nonnegative finite number.", call. = FALSE)
  }
  x
}

#' @noRd
.check_annuity_vp_n <- function(n) {
  n <- as.numeric(n)
  if (length(n) != 1L || !is.finite(n) || n < 0) {
    stop("n must be a single nonnegative finite number.", call. = FALSE)
  }
  n
}

#' @noRd
.recycle_x_n_vp <- function(x, n) {
  x <- as.numeric(x)
  n <- as.numeric(n)

  if (length(x) == 0L || length(n) == 0L) {
    stop("x and n must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("x must contain nonnegative finite values.", call. = FALSE)
  }
  if (any(!is.finite(n)) || any(n < 0)) {
    stop("n must contain nonnegative finite values.", call. = FALSE)
  }

  if (length(x) == 1L && length(n) > 1L) x <- rep(x, length(n))
  if (length(n) == 1L && length(x) > 1L) n <- rep(n, length(x))

  if (length(x) != length(n)) {
    stop("x and n must have the same length, or one must have length 1.", call. = FALSE)
  }

  list(x = x, n = n)
}

#' @noRd
.max_term_vp <- function(x, model, ...) {
  p <- list(...)
  model <- tolower(model)

  if (model == "uniform") {
    omega <- p$omega
    if (is.null(omega)) {
      stop("For model = 'uniform', omega must be supplied.", call. = FALSE)
    }
    return(max(0, floor(omega - x)))
  }

  Inf
}

#' @noRd
.integrate_vp <- function(f, lower, upper) {
  stats::integrate(
    f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-9,
    subdivisions = 1000L
  )$value
}

#' @noRd
.upper_cont_vp <- function(x, n = NULL, model, ..., tol = 1e-10) {
  max_term <- .max_term_vp(x, model, ...)
  if (is.finite(max_term)) {
    if (is.null(n)) return(max_term)
    return(min(max_term, n))
  }

  if (!is.null(n)) return(n)

  upper <- 10
  repeat {
    surv <- tpx(upper, x = x, model = model, ...)
    if (!is.finite(surv) || surv <= tol || upper >= 5000) break
    upper <- upper * 2
  }
  upper
}

# -------------------------------------------------------------------------
# Shared documentation block
# -------------------------------------------------------------------------

#' Increasing and decreasing life annuities
#'
#' Increasing and decreasing annuity functions from Chapter 8.
#'
#' @aliases Iax Iaxn Daxn
#' @aliases Iadotx Iadotxn Dadotxn
#' @aliases Iabarx Iabarxn Dabarxn
#' @rdname annuity_varying_payments
NULL

# -------------------------------------------------------------------------
# Annual immediate
# -------------------------------------------------------------------------

#' Increasing whole life annuity-immediate
#'
#' Computes
#' \deqn{(Ia)_x = \sum_{t=1}^{\infty} t \, v^t \, {}_tp_x.}
#'
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return Numeric vector.
#' @rdname annuity_varying_payments
#' @export
Iax <- function(x, i, model, ..., k_max = 5000, tol = 1e-12) {
  i <- .check_annuity_vp_i(i)
  x <- as.numeric(x)
  v <- 1 / (1 + i)

  sapply(x, function(xx) {
    xx <- .check_annuity_vp_x(xx)
    max_term <- .max_term_vp(xx, model, ...)

    times <- if (is.finite(max_term)) {
      if (max_term <= 0) integer(0) else 1:max_term
    } else {
      1:k_max
    }

    if (length(times) == 0L) return(0)

    terms <- times * (v^times) * tpx(times, x = xx, model = model, ...)

    if (is.finite(max_term)) return(sum(terms))

    acc <- 0
    small_run <- 0L
    for (j in seq_along(terms)) {
      acc <- acc + terms[j]
      if (terms[j] < tol) {
        small_run <- small_run + 1L
        if (small_run >= 10L) break
      } else {
        small_run <- 0L
      }
    }
    acc
  })
}

#' Increasing temporary annuity-immediate
#'
#' Computes
#' \deqn{(Ia)_{x:\angl{n}} = \sum_{t=1}^{n} t \, v^t \, {}_tp_x.}
#'
#' @param n Term in years.
#' @return Numeric vector.
#' @rdname annuity_varying_payments
#' @export
Iaxn <- function(x, n, i, model, ...) {
  i <- .check_annuity_vp_i(i)
  xn <- .recycle_x_n_vp(x, n)
  x <- xn$x
  n <- xn$n
  v <- 1 / (1 + i)

  sapply(seq_along(x), function(j) {
    xx <- .check_annuity_vp_x(x[j])
    nn <- .check_annuity_vp_n(n[j])

    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))
    if (nn == 0L) return(0)

    times <- 1:nn
    sum(times * (v^times) * tpx(times, x = xx, model = model, ...))
  })
}

#' Decreasing temporary annuity-immediate
#'
#' Computes
#' \deqn{(Da)_{x:\angl{n}} = \sum_{t=1}^{n} (n+1-t)\, v^t \, {}_tp_x.}
#'
#' @return Numeric vector.
#' @rdname annuity_varying_payments
#' @export
Daxn <- function(x, n, i, model, ...) {
  i <- .check_annuity_vp_i(i)
  xn <- .recycle_x_n_vp(x, n)
  x <- xn$x
  n <- xn$n
  v <- 1 / (1 + i)

  sapply(seq_along(x), function(j) {
    xx <- .check_annuity_vp_x(x[j])
    nn <- .check_annuity_vp_n(n[j])

    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))
    if (nn == 0L) return(0)

    times <- 1:nn
    sum((nn + 1 - times) * (v^times) * tpx(times, x = xx, model = model, ...))
  })
}

# -------------------------------------------------------------------------
# Annual due
# -------------------------------------------------------------------------

#' Increasing whole life annuity-due
#'
#' Computes
#' \deqn{(I\ddot{a})_x = \sum_{t=0}^{\infty} (t+1)\, v^t \, {}_tp_x.}
#'
#' @rdname annuity_varying_payments
#' @export
Iadotx <- function(x, i, model, ..., k_max = 5000, tol = 1e-12) {
  i <- .check_annuity_vp_i(i)
  x <- as.numeric(x)
  v <- 1 / (1 + i)

  sapply(x, function(xx) {
    xx <- .check_annuity_vp_x(xx)
    max_term <- .max_term_vp(xx, model, ...)

    times <- if (is.finite(max_term)) {
      0:max_term
    } else {
      0:k_max
    }

    terms <- (times + 1) * (v^times) * tpx(times, x = xx, model = model, ...)

    if (is.finite(max_term)) return(sum(terms))

    acc <- 0
    small_run <- 0L
    for (j in seq_along(terms)) {
      acc <- acc + terms[j]
      if (terms[j] < tol) {
        small_run <- small_run + 1L
        if (small_run >= 10L) break
      } else {
        small_run <- 0L
      }
    }
    acc
  })
}

#' Increasing temporary annuity-due
#'
#' Computes
#' \deqn{(I\ddot{a})_{x:\angl{n}} = \sum_{t=0}^{n-1} (t+1)\, v^t \, {}_tp_x.}
#'
#' @rdname annuity_varying_payments
#' @export
Iadotxn <- function(x, n, i, model, ...) {
  i <- .check_annuity_vp_i(i)
  xn <- .recycle_x_n_vp(x, n)
  x <- xn$x
  n <- xn$n
  v <- 1 / (1 + i)

  sapply(seq_along(x), function(j) {
    xx <- .check_annuity_vp_x(x[j])
    nn <- .check_annuity_vp_n(n[j])

    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))
    if (nn == 0L) return(0)

    times <- 0:(nn - 1L)
    sum((times + 1) * (v^times) * tpx(times, x = xx, model = model, ...))
  })
}

#' Decreasing temporary annuity-due
#'
#' Computes
#' \deqn{(D\ddot{a})_{x:\angl{n}} = \sum_{t=0}^{n-1} (n-t)\, v^t \, {}_tp_x.}
#'
#' @rdname annuity_varying_payments
#' @export
Dadotxn <- function(x, n, i, model, ...) {
  i <- .check_annuity_vp_i(i)
  xn <- .recycle_x_n_vp(x, n)
  x <- xn$x
  n <- xn$n
  v <- 1 / (1 + i)

  sapply(seq_along(x), function(j) {
    xx <- .check_annuity_vp_x(x[j])
    nn <- .check_annuity_vp_n(n[j])

    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))
    if (nn == 0L) return(0)

    times <- 0:(nn - 1L)
    sum((nn - times) * (v^times) * tpx(times, x = xx, model = model, ...))
  })
}

# -------------------------------------------------------------------------
# Continuous
# -------------------------------------------------------------------------

#' Increasing continuous whole life annuity
#'
#' Computes
#' \deqn{(\bar{I}\bar{a})_x = \int_0^\infty t\,v^t\,{}_tp_x\,dt.}
#'
#' @rdname annuity_varying_payments
#' @export
Iabarx <- function(x, i, model, ..., tol = 1e-10) {
  i <- .check_annuity_vp_i(i)
  x <- as.numeric(x)
  delta <- log(1 + i)

  sapply(x, function(xx) {
    xx <- .check_annuity_vp_x(xx)
    upper <- .upper_cont_vp(xx, model = model, ..., tol = tol)
    if (upper <= 0) return(0)

    f <- function(t) t * exp(-delta * t) * tpx(t, x = xx, model = model, ...)
    .integrate_vp(f, 0, upper)
  })
}

#' Increasing continuous temporary annuity
#'
#' Computes
#' \deqn{(\bar{I}\bar{a})_{x:\angl{n}} = \int_0^n t\,v^t\,{}_tp_x\,dt.}
#'
#' @rdname annuity_varying_payments
#' @export
Iabarxn <- function(x, n, i, model, ...) {
  i <- .check_annuity_vp_i(i)
  xn <- .recycle_x_n_vp(x, n)
  x <- xn$x
  n <- xn$n
  delta <- log(1 + i)

  sapply(seq_along(x), function(j) {
    xx <- .check_annuity_vp_x(x[j])
    nn <- .check_annuity_vp_n(n[j])
    if (nn == 0) return(0)

    f <- function(t) t * exp(-delta * t) * tpx(t, x = xx, model = model, ...)
    .integrate_vp(f, 0, nn)
  })
}

#' Decreasing continuous temporary annuity
#'
#' Computes
#' \deqn{(\bar{D}\bar{a})_{x:\angl{n}} = \int_0^n (n-t)\,v^t\,{}_tp_x\,dt.}
#'
#' @rdname annuity_varying_payments
#' @export
Dabarxn <- function(x, n, i, model, ...) {
  i <- .check_annuity_vp_i(i)
  xn <- .recycle_x_n_vp(x, n)
  x <- xn$x
  n <- xn$n
  delta <- log(1 + i)

  sapply(seq_along(x), function(j) {
    xx <- .check_annuity_vp_x(x[j])
    nn <- .check_annuity_vp_n(n[j])
    if (nn == 0) return(0)

    f <- function(t) (nn - t) * exp(-delta * t) * tpx(t, x = xx, model = model, ...)
    .integrate_vp(f, 0, nn)
  })
}
