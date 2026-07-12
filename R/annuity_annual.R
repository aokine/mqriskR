#' Annual annuity functions
#'
#' Annual whole life, temporary, deferred, and actuarial accumulated value
#' annuity functions in immediate, due, and continuous forms.
#'
#' @param x Age.
#' @param n Term in years.
#' @param i Effective annual interest rate.
#' @param model Optional survival model name.
#' @param ... Additional model parameters.
#' @param tbl Optional life table object for annual discrete annuity functions.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return
#' Numeric vector containing the requested annuity present value or actuarial
#' accumulated value.
#'
#' @name annuity_annual
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_annuity_i <- function(i) {
  i <- as.numeric(i)

  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }

  i
}

.check_n_scalar <- function(n) {
  n <- as.numeric(n)

  if (length(n) != 1L || !is.finite(n) || n < 0) {
    stop("n must be a single nonnegative number.", call. = FALSE)
  }

  n
}

.check_integer_term <- function(n) {
  n <- .check_n_scalar(n)

  if (abs(n - round(n)) > 1e-10) {
    stop("For annual annuity functions, n must be an integer.", call. = FALSE)
  }

  as.integer(round(n))
}

.recycle_annuity_vectors <- function(...) {
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

.recycle_xn <- function(x, n) {
  args <- .recycle_annuity_vectors(x, n)
  list(x = args[[1L]], n = args[[2L]])
}

.check_age_annuity <- function(x) {
  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  invisible(TRUE)
}

.get_npx_annuity <- function(x, n, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(npx(tbl, x = x, n = n))
  }

  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }

  tpx(n, x = x, model = model, ...)
}

.max_term_from_source <- function(x, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(max(0, floor(max(tbl$x, na.rm = TRUE) - x)))
  }

  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }

  dots <- list(...)
  model <- tolower(model)

  if (model == "uniform") {
    return(max(0, floor(dots$omega - x)))
  }

  Inf
}

.sum_annuity_terms <- function(x, i, tbl = NULL, model = NULL, ...,
                               k_max = 5000, tol = 1e-12,
                               due = FALSE, temporary = NULL) {
  x <- as.numeric(x)
  i <- .check_annuity_i(i)

  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stop("x must be a single nonnegative finite number.", call. = FALSE)
  }

  if (!is.null(temporary)) {
    temporary <- .check_integer_term(temporary)
  }

  v <- 1 / (1 + i)

  if (due) {
    if (!is.null(temporary)) {
      if (temporary == 0L) return(0)
      times <- 0:(temporary - 1L)
    } else {
      max_term <- .max_term_from_source(x, tbl = tbl, model = model, ...)
      times <- if (is.finite(max_term)) 0:max_term else 0:k_max
    }
  } else {
    if (!is.null(temporary)) {
      if (temporary == 0L) return(0)
      times <- 1:temporary
    } else {
      max_term <- .max_term_from_source(x, tbl = tbl, model = model, ...)
      if (is.finite(max_term)) {
        if (max_term <= 0) return(0)
        times <- 1:max_term
      } else {
        times <- 1:k_max
      }
    }
  }

  surv <- .get_npx_annuity(x, times, tbl = tbl, model = model, ...)
  surv[is.na(surv)] <- 0

  terms <- v^times * surv

  if (!is.null(temporary) ||
      is.finite(.max_term_from_source(x, tbl = tbl, model = model, ...))) {
    return(sum(terms, na.rm = TRUE))
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
#' Computes \eqn{a_x = \sum_{t=1}^{\infty} v^t {}_t p_x}.
#'
#' @rdname annuity_annual
#' @export
ax <- function(x, i, model = NULL, ..., tbl = NULL,
               k_max = 5000, tol = 1e-12) {
  args <- .recycle_annuity_vectors(x, .check_annuity_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    .sum_annuity_terms(
      x[j],
      i[j],
      tbl = tbl,
      model = model,
      ...,
      k_max = k_max,
      tol = tol,
      due = FALSE
    )
  })
}

#' Whole life annuity-due
#'
#' Computes \eqn{\ddot{a}_x = \sum_{t=0}^{\infty} v^t {}_t p_x}.
#'
#' @rdname annuity_annual
#' @export
adotx <- function(x, i, model = NULL, ..., tbl = NULL,
                  k_max = 5000, tol = 1e-12) {
  args <- .recycle_annuity_vectors(x, .check_annuity_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    .sum_annuity_terms(
      x[j],
      i[j],
      tbl = tbl,
      model = model,
      ...,
      k_max = k_max,
      tol = tol,
      due = TRUE
    )
  })
}

#' Continuous whole life annuity
#'
#' Computes \eqn{\bar{a}_x = \int_0^{\infty} v^t {}_t p_x dt}.
#'
#' @rdname annuity_annual
#' @export
abarx <- function(x, i, model, ..., tol = 1e-10) {
  args <- .recycle_annuity_vectors(x, .check_annuity_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    delta <- log(1 + i[j])

    max_term <- .max_term_from_source(xx, model = model, ...)

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
#' Computes \eqn{a_{x:\overline{n}|} = \sum_{t=1}^{n} v^t {}_t p_x}.
#'
#' @rdname annuity_annual
#' @export
axn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_integer_term(n[j])

    .sum_annuity_terms(
      x[j],
      i[j],
      tbl = tbl,
      model = model,
      ...,
      due = FALSE,
      temporary = nn
    )
  })
}

#' Temporary annuity-due
#'
#' Computes \eqn{\ddot{a}_{x:\overline{n}|} =
#' \sum_{t=0}^{n-1} v^t {}_t p_x}.
#'
#' @rdname annuity_annual
#' @export
adotxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_integer_term(n[j])

    .sum_annuity_terms(
      x[j],
      i[j],
      tbl = tbl,
      model = model,
      ...,
      due = TRUE,
      temporary = nn
    )
  })
}

#' Continuous temporary annuity
#'
#' Computes \eqn{\bar{a}_{x:\overline{n}|} =
#' \int_0^n v^t {}_t p_x dt}.
#'
#' @rdname annuity_annual
#' @export
abarxn <- function(x, n, i, model, ...) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    delta <- log(1 + i[j])

    f <- function(t) {
      exp(-delta * t) * tpx(t, x = x[j], model = model, ...)
    }

    .integrate_annuity(f, 0, nn)
  })
}

# -------------------------------------------------------------------------
# Deferred annuities
# -------------------------------------------------------------------------

#' Deferred whole life annuity-immediate
#'
#' Computes \eqn{{}_{n|}a_x = {}_nE_x a_{x+n}}.
#'
#' @rdname annuity_annual
#' @export
nax <- function(x, n, i, model = NULL, ..., tbl = NULL,
                k_max = 5000, tol = 1e-12) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_integer_term(n[j])

    nEx(x[j], nn, i[j], tbl = tbl, model = model, ...) *
      ax(
        x[j] + nn,
        i[j],
        tbl = tbl,
        model = model,
        ...,
        k_max = k_max,
        tol = tol
      )
  })
}

#' Deferred whole life annuity-due
#'
#' Computes \eqn{{}_{n|}\ddot{a}_x = {}_nE_x \ddot{a}_{x+n}}.
#'
#' @rdname annuity_annual
#' @export
nadotx <- function(x, n, i, model = NULL, ..., tbl = NULL,
                   k_max = 5000, tol = 1e-12) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_integer_term(n[j])

    nEx(x[j], nn, i[j], tbl = tbl, model = model, ...) *
      adotx(
        x[j] + nn,
        i[j],
        tbl = tbl,
        model = model,
        ...,
        k_max = k_max,
        tol = tol
      )
  })
}

#' Deferred continuous whole life annuity
#'
#' Computes \eqn{{}_{n|}\bar{a}_x = {}_nE_x \bar{a}_{x+n}}.
#'
#' @rdname annuity_annual
#' @export
nabarx <- function(x, n, i, model, ..., tol = 1e-10) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])

    nEx(x[j], nn, i[j], model = model, ...) *
      abarx(x[j] + nn, i[j], model = model, ..., tol = tol)
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
sxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_integer_term(n[j])
    ex <- nEx(x[j], nn, i[j], tbl = tbl, model = model, ...)

    if (is.na(ex) || ex <= 0) return(Inf)

    axn(x[j], nn, i[j], tbl = tbl, model = model, ...) / ex
  })
}

#' Temporary annuity-due actuarial accumulated value
#'
#' Computes \eqn{\ddot{s}_{x:\overline{n}|} =
#' \ddot{a}_{x:\overline{n}|} / {}_nE_x}.
#'
#' @rdname annuity_annual
#' @export
sdotxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_integer_term(n[j])
    ex <- nEx(x[j], nn, i[j], tbl = tbl, model = model, ...)

    if (is.na(ex) || ex <= 0) return(Inf)

    adotxn(x[j], nn, i[j], tbl = tbl, model = model, ...) / ex
  })
}

#' Continuous temporary annuity actuarial accumulated value
#'
#' Computes \eqn{\bar{s}_{x:\overline{n}|} =
#' \bar{a}_{x:\overline{n}|} / {}_nE_x}.
#'
#' @rdname annuity_annual
#' @export
sbarxn <- function(x, n, i, model, ...) {
  args <- .recycle_annuity_vectors(x, n, .check_annuity_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_age_annuity(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_n_scalar(n[j])
    ex <- nEx(x[j], nn, i[j], model = model, ...)

    if (is.na(ex) || ex <= 0) return(Inf)

    abarxn(x[j], nn, i[j], model = model, ...) / ex
  })
}
