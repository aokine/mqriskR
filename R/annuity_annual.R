#' Annual annuity functions (Chapter 8)
#'
#' Annual whole life, temporary, deferred, and actuarial accumulated value
#' annuity functions in immediate, due, and continuous forms.
#'
#' Naming convention follows Chapter 8 notation:
#' \itemize{
#'   \item \code{ax()} = \eqn{a_x}
#'   \item \code{adotx()} = \eqn{\ddot{a}_x}
#'   \item \code{abarx()} = \eqn{\bar{a}_x}
#'   \item \code{axn()} = \eqn{a_{x:\overline{n}|}}
#'   \item \code{adotxn()} = \eqn{\ddot{a}_{x:\overline{n}|}}
#'   \item \code{abarxn()} = \eqn{\bar{a}_{x:\overline{n}|}}
#'   \item \code{nax()} = \eqn{{}_{n|}a_x}
#'   \item \code{nadotx()} = \eqn{{}_{n|}\ddot{a}_x}
#'   \item \code{nabarx()} = \eqn{{}_{n|}\bar{a}_x}
#'   \item \code{sxn()} = \eqn{s_{x:\overline{n}|}}
#'   \item \code{sdotxn()} = \eqn{\ddot{s}_{x:\overline{n}|}}
#'   \item \code{sbarxn()} = \eqn{\bar{s}_{x:\overline{n}|}}
#' }
#'
#' These functions work directly from the Chapter 5 survival model functions
#' and the Chapter 7 pure endowment function \code{nEx()}.
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
#' @name annuity_annual
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_annuity_i <- function(i) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite number greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.check_n_scalar <- function(n) {
  n <- as.numeric(n)
  if (length(n) != 1L || !is.finite(n) || n < 0) {
    stop("n must be a single nonnegative number.", call. = FALSE)
  }
  n
}

#' @noRd
.recycle_xn <- function(x, n) {
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

#' @noRd
.max_term_from_model <- function(x, model, ...) {
  p <- list(...)
  model <- tolower(model)

  if (model == "uniform") {
    omega <- p$omega
    return(max(0, floor(omega - x)))
  }

  Inf
}

#' @noRd
.sum_annuity_terms <- function(x, i, model, ..., k_max = 5000, tol = 1e-12,
                               due = FALSE, temporary = NULL, start = 1) {
  # General annual annuity summation:
  # immediate whole life: due = FALSE, temporary = NULL, start = 1
  # due whole life      : due = TRUE,  temporary = NULL, start = 0
  # temporary immediate : due = FALSE, temporary = n,    start = 1
  # temporary due       : due = TRUE,  temporary = n,    start = 0

  .check_annuity_i(i)
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stop("x must be a single nonnegative finite number.", call. = FALSE)
  }

  if (!is.null(temporary)) {
    temporary <- .check_n_scalar(temporary)
    if (abs(temporary - round(temporary)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    temporary <- as.integer(round(temporary))
  }

  v <- 1 / (1 + i)

  if (due) {
    if (!is.null(temporary)) {
      if (temporary == 0L) return(0)
      times <- 0:(temporary - 1L)
    } else {
      max_term <- .max_term_from_model(x, model, ...)
      if (is.finite(max_term)) {
        times <- 0:max_term
      } else {
        times <- 0:k_max
      }
    }
  } else {
    if (!is.null(temporary)) {
      if (temporary == 0L) return(0)
      times <- 1:temporary
    } else {
      max_term <- .max_term_from_model(x, model, ...)
      if (is.finite(max_term)) {
        if (max_term <= 0) return(0)
        times <- 1:max_term
      } else {
        times <- 1:k_max
      }
    }
  }

  surv <- tpx(times, x = x, model = model, ...)
  terms <- v^times * surv

  if (!is.null(temporary) || is.finite(.max_term_from_model(x, model, ...))) {
    return(sum(terms))
  }

  acc <- 0
  small_run <- 0
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
}

#' @noRd
.integrate_annuity <- function(f, lower, upper) {
  stats::integrate(
    f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-9,
    subdivisions = 1000L
  )$value
}

# -------------------------------------------------------------------------
# Whole life annuities
# -------------------------------------------------------------------------

#' Whole life annuity-immediate
#'
#' Computes the annual whole life annuity-immediate
#' \eqn{a_x = \sum_{t=1}^{\infty} v^t \, {}_t p_x}.
#'
#' @rdname annuity_annual
#' @export
ax <- function(x, i, model, ..., k_max = 5000, tol = 1e-12) {
  .check_annuity_i(i)
  x <- as.numeric(x)
  sapply(x, function(xx) {
    .sum_annuity_terms(xx, i, model, ..., k_max = k_max, tol = tol, due = FALSE)
  })
}

#' Whole life annuity-due
#'
#' Computes the annual whole life annuity-due
#' \eqn{\ddot{a}_x = \sum_{t=0}^{\infty} v^t \, {}_t p_x = 1 + a_x}.
#'
#' @rdname annuity_annual
#' @export
adotx <- function(x, i, model, ..., k_max = 5000, tol = 1e-12) {
  .check_annuity_i(i)
  x <- as.numeric(x)
  sapply(x, function(xx) {
    .sum_annuity_terms(xx, i, model, ..., k_max = k_max, tol = tol, due = TRUE)
  })
}

#' Continuous whole life annuity
#'
#' Computes the continuous whole life annuity
#' \eqn{\bar{a}_x = \int_0^{\infty} v^t \, {}_t p_x \, dt}.
#'
#' @rdname annuity_annual
#' @export
abarx <- function(x, i, model, ..., tol = 1e-10) {
  .check_annuity_i(i)
  x <- as.numeric(x)
  delta <- log(1 + i)

  sapply(x, function(xx) {
    if (!is.finite(xx) || xx < 0) {
      stop("x must be nonnegative and finite.", call. = FALSE)
    }

    max_term <- .max_term_from_model(xx, model, ...)

    if (is.finite(max_term)) {
      if (max_term <= 0) return(0)
      f <- function(t) exp(-delta * t) * tpx(t, x = xx, model = model, ...)
      return(.integrate_annuity(f, 0, max_term))
    }

    if (tolower(model) == "exponential") {
      lambda <- list(...)$lambda
      return(1 / (delta + lambda))
    }

    upper <- 10
    while (upper < 5000 && tpx(upper, xx, model = model, ...) > tol) {
      upper <- upper * 2
    }

    f <- function(t) exp(-delta * t) * tpx(t, x = xx, model = model, ...)
    .integrate_annuity(f, 0, upper)
  })
}

# -------------------------------------------------------------------------
# Temporary annuities
# -------------------------------------------------------------------------

#' Temporary annuity-immediate
#'
#' Computes the annual temporary annuity-immediate
#' \eqn{a_{x:\overline{n}|} = \sum_{t=1}^{n} v^t \, {}_t p_x}.
#'
#' @rdname annuity_annual
#' @export
axn <- function(x, n, i, model, ...) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    .sum_annuity_terms(
      x[j], i, model, ...,
      due = FALSE,
      temporary = as.integer(round(nn))
    )
  })
}

#' Temporary annuity-due
#'
#' Computes the annual temporary annuity-due
#' \eqn{\ddot{a}_{x:\overline{n}|} = \sum_{t=0}^{n-1} v^t \, {}_t p_x}.
#'
#' @rdname annuity_annual
#' @export
adotxn <- function(x, n, i, model, ...) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    .sum_annuity_terms(
      x[j], i, model, ...,
      due = TRUE,
      temporary = as.integer(round(nn))
    )
  })
}

#' Continuous temporary annuity
#'
#' Computes the continuous temporary annuity
#' \eqn{\bar{a}_{x:\overline{n}|} = \int_0^n v^t \, {}_t p_x \, dt}.
#'
#' @rdname annuity_annual
#' @export
abarxn <- function(x, n, i, model, ...) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n
  delta <- log(1 + i)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- .check_n_scalar(n[j])
    f <- function(t) exp(-delta * t) * tpx(t, x = xx, model = model, ...)
    .integrate_annuity(f, 0, nn)
  })
}

# -------------------------------------------------------------------------
# Deferred annuities
# -------------------------------------------------------------------------

#' Deferred whole life annuity-immediate
#'
#' Computes the annual deferred whole life annuity-immediate
#' \eqn{{}_{n|}a_x = {}_nE_x \, a_{x+n}}.
#'
#' @rdname annuity_annual
#' @export
nax <- function(x, n, i, model, ..., k_max = 5000, tol = 1e-12) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    nEx(x[j], nn, i, model = model, ...) *
      ax(x[j] + nn, i, model = model, ..., k_max = k_max, tol = tol)
  })
}

#' Deferred whole life annuity-due
#'
#' Computes the annual deferred whole life annuity-due
#' \eqn{{}_{n|}\ddot{a}_x = {}_nE_x \, \ddot{a}_{x+n}}.
#'
#' @rdname annuity_annual
#' @export
nadotx <- function(x, n, i, model, ..., k_max = 5000, tol = 1e-12) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    nEx(x[j], nn, i, model = model, ...) *
      adotx(x[j] + nn, i, model = model, ..., k_max = k_max, tol = tol)
  })
}

#' Deferred continuous whole life annuity
#'
#' Computes the continuous deferred whole life annuity
#' \eqn{{}_{n|}\bar{a}_x = {}_nE_x \, \bar{a}_{x+n}}.
#'
#' @rdname annuity_annual
#' @export
nabarx <- function(x, n, i, model, ..., tol = 1e-10) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    nEx(x[j], nn, i, model = model, ...) *
      abarx(x[j] + nn, i, model = model, ..., tol = tol)
  })
}

# -------------------------------------------------------------------------
# Actuarial accumulated values
# -------------------------------------------------------------------------

#' Temporary annuity-immediate actuarial accumulated value
#'
#' Computes \eqn{s_{x:\overline{n}|} = a_{x:\overline{n}|} / {}_nE_x}.
#'
#' @rdname annuity_annual
#' @export
sxn <- function(x, n, i, model, ...) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    ex <- nEx(x[j], nn, i, model = model, ...)
    if (ex <= 0) return(Inf)
    axn(x[j], nn, i, model = model, ...) / ex
  })
}

#' Temporary annuity-due actuarial accumulated value
#'
#' Computes \eqn{\ddot{s}_{x:\overline{n}|} = \ddot{a}_{x:\overline{n}|} / {}_nE_x}.
#'
#' @rdname annuity_annual
#' @export
sdotxn <- function(x, n, i, model, ...) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    if (abs(nn - round(nn)) > 1e-10) {
      stop("For annual annuity functions, n must be an integer.", call. = FALSE)
    }
    nn <- as.integer(round(nn))

    ex <- nEx(x[j], nn, i, model = model, ...)
    if (ex <= 0) return(Inf)
    adotxn(x[j], nn, i, model = model, ...) / ex
  })
}

#' Continuous temporary annuity actuarial accumulated value
#'
#' Computes \eqn{\bar{s}_{x:\overline{n}|} = \bar{a}_{x:\overline{n}|} / {}_nE_x}.
#'
#' @rdname annuity_annual
#' @export
sbarxn <- function(x, n, i, model, ...) {
  .check_annuity_i(i)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    ex <- nEx(x[j], nn, i, model = model, ...)
    if (ex <= 0) return(Inf)
    abarxn(x[j], nn, i, model = model, ...) / ex
  })
}
