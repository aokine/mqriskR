#' Varying-payment annuity functions
#'
#' Increasing and decreasing life annuity functions.
#'
#' @param x Age.
#' @param n Term in years.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object for discrete functions.
#' @param model Optional survival model name.
#' @param ... Additional model parameters.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return
#' Numeric vector containing the requested increasing or decreasing annuity value.
#'
#' @name annuity_varying_payments
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_annuity_vp_i <- function(i) {
  i <- as.numeric(i)
  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }
  i
}

.recycle_vp_vectors <- function(...) {
  args <- lapply(list(...), as.numeric)
  lens <- vapply(args, length, integer(1))

  if (any(lens == 0L)) stop("All arguments must have positive length.", call. = FALSE)
  if (any(vapply(args, function(z) any(!is.finite(z)), logical(1)))) {
    stop("All arguments must contain finite numeric values.", call. = FALSE)
  }

  common_len <- max(lens)
  if (any(!(lens %in% c(1L, common_len)))) {
    stop("Arguments must have compatible lengths.", call. = FALSE)
  }

  lapply(args, function(z) if (length(z) == 1L) rep(z, common_len) else z)
}

.check_vp_x <- function(x) {
  if (any(x < 0)) stop("x must be nonnegative.", call. = FALSE)
  invisible(TRUE)
}

.check_vp_n <- function(n) {
  if (any(n < 0)) stop("n must be nonnegative.", call. = FALSE)
  invisible(TRUE)
}

.check_vp_integer_n <- function(n) {
  if (length(n) != 1L || !is.finite(n) || n < 0) {
    stop("n must be a single nonnegative finite number.", call. = FALSE)
  }
  if (abs(n - round(n)) > 1e-10) {
    stop("For annual annuity functions, n must be an integer.", call. = FALSE)
  }
  as.integer(round(n))
}

.get_vp_npx <- function(x, n, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(npx(tbl, x = x, n = n))
  }
  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  tpx(n, x = x, model = model, ...)
}

.max_term_vp <- function(x, tbl = NULL, model = NULL, ...) {
  if (!is.null(tbl)) {
    return(max(0, floor(max(tbl$x, na.rm = TRUE) - x)))
  }

  if (is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }

  p <- list(...)
  model <- tolower(model)

  if (model == "uniform") {
    if (is.null(p$omega)) {
      stop("For model = 'uniform', omega must be supplied.", call. = FALSE)
    }
    return(max(0, floor(p$omega - x)))
  }

  Inf
}

.integrate_vp <- function(f, lower, upper) {
  stats::integrate(
    f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-9,
    subdivisions = 1000L
  )$value
}

.upper_cont_vp <- function(x, n = NULL, model, ..., tol = 1e-10) {
  max_term <- .max_term_vp(x, model = model, ...)
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
# Annual immediate
# -------------------------------------------------------------------------

#' Increasing whole life annuity-immediate
#'
#' @rdname annuity_varying_payments
#' @export
Iax <- function(x, i, model = NULL, ..., tbl = NULL,
                k_max = 5000, tol = 1e-12) {
  args <- .recycle_vp_vectors(x, .check_annuity_vp_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_vp_x(x)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    v <- 1 / (1 + i[j])
    max_term <- .max_term_vp(xx, tbl = tbl, model = model, ...)

    times <- if (is.finite(max_term)) {
      if (max_term <= 0) integer(0) else 1:max_term
    } else {
      1:k_max
    }

    if (length(times) == 0L) return(0)

    surv <- .get_vp_npx(xx, times, tbl = tbl, model = model, ...)
    surv[is.na(surv)] <- 0

    terms <- times * v^times * surv

    if (is.finite(max_term)) return(sum(terms))

    acc <- 0
    small_run <- 0L
    for (jj in seq_along(terms)) {
      acc <- acc + terms[jj]
      if (terms[jj] < tol) {
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
#' @rdname annuity_varying_payments
#' @export
Iaxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_vp_vectors(x, n, .check_annuity_vp_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_vp_x(x)
  .check_vp_n(n)

  sapply(seq_along(x), function(j) {
    nn <- .check_vp_integer_n(n[j])
    if (nn == 0L) return(0)

    v <- 1 / (1 + i[j])
    times <- 1:nn

    surv <- .get_vp_npx(x[j], times, tbl = tbl, model = model, ...)
    surv[is.na(surv)] <- 0

    sum(times * v^times * surv)
  })
}

#' Decreasing temporary annuity-immediate
#'
#' @rdname annuity_varying_payments
#' @export
Daxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_vp_vectors(x, n, .check_annuity_vp_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_vp_x(x)
  .check_vp_n(n)

  sapply(seq_along(x), function(j) {
    nn <- .check_vp_integer_n(n[j])
    if (nn == 0L) return(0)

    v <- 1 / (1 + i[j])
    times <- 1:nn

    surv <- .get_vp_npx(x[j], times, tbl = tbl, model = model, ...)
    surv[is.na(surv)] <- 0

    sum((nn + 1 - times) * v^times * surv)
  })
}

# -------------------------------------------------------------------------
# Annual due
# -------------------------------------------------------------------------

#' Increasing whole life annuity-due
#'
#' @rdname annuity_varying_payments
#' @export
Iadotx <- function(x, i, model = NULL, ..., tbl = NULL,
                   k_max = 5000, tol = 1e-12) {
  args <- .recycle_vp_vectors(x, .check_annuity_vp_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_vp_x(x)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    v <- 1 / (1 + i[j])
    max_term <- .max_term_vp(xx, tbl = tbl, model = model, ...)

    times <- if (is.finite(max_term)) 0:max_term else 0:k_max

    surv <- .get_vp_npx(xx, times, tbl = tbl, model = model, ...)
    surv[is.na(surv)] <- 0

    terms <- (times + 1) * v^times * surv

    if (is.finite(max_term)) return(sum(terms))

    acc <- 0
    small_run <- 0L
    for (jj in seq_along(terms)) {
      acc <- acc + terms[jj]
      if (terms[jj] < tol) {
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
#' @rdname annuity_varying_payments
#' @export
Iadotxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_vp_vectors(x, n, .check_annuity_vp_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_vp_x(x)
  .check_vp_n(n)

  sapply(seq_along(x), function(j) {
    nn <- .check_vp_integer_n(n[j])
    if (nn == 0L) return(0)

    v <- 1 / (1 + i[j])
    times <- 0:(nn - 1L)

    surv <- .get_vp_npx(x[j], times, tbl = tbl, model = model, ...)
    surv[is.na(surv)] <- 0

    sum((times + 1) * v^times * surv)
  })
}

#' Decreasing temporary annuity-due
#'
#' @rdname annuity_varying_payments
#' @export
Dadotxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_vp_vectors(x, n, .check_annuity_vp_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_vp_x(x)
  .check_vp_n(n)

  sapply(seq_along(x), function(j) {
    nn <- .check_vp_integer_n(n[j])
    if (nn == 0L) return(0)

    v <- 1 / (1 + i[j])
    times <- 0:(nn - 1L)

    surv <- .get_vp_npx(x[j], times, tbl = tbl, model = model, ...)
    surv[is.na(surv)] <- 0

    sum((nn - times) * v^times * surv)
  })
}

# -------------------------------------------------------------------------
# Continuous
# -------------------------------------------------------------------------

#' Increasing continuous whole life annuity
#'
#' @rdname annuity_varying_payments
#' @export
Iabarx <- function(x, i, model, ..., tol = 1e-10) {
  args <- .recycle_vp_vectors(x, .check_annuity_vp_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_vp_x(x)

  sapply(seq_along(x), function(j) {
    delta <- log(1 + i[j])
    upper <- .upper_cont_vp(x[j], model = model, ..., tol = tol)

    if (upper <= 0) return(0)

    f <- function(t) {
      t * exp(-delta * t) * tpx(t, x = x[j], model = model, ...)
    }

    .integrate_vp(f, 0, upper)
  })
}

#' Increasing continuous temporary annuity
#'
#' @rdname annuity_varying_payments
#' @export
Iabarxn <- function(x, n, i, model, ...) {
  args <- .recycle_vp_vectors(x, n, .check_annuity_vp_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_vp_x(x)
  .check_vp_n(n)

  sapply(seq_along(x), function(j) {
    if (n[j] == 0) return(0)

    delta <- log(1 + i[j])

    f <- function(t) {
      t * exp(-delta * t) * tpx(t, x = x[j], model = model, ...)
    }

    .integrate_vp(f, 0, n[j])
  })
}

#' Decreasing continuous temporary annuity
#'
#' @rdname annuity_varying_payments
#' @export
Dabarxn <- function(x, n, i, model, ...) {
  args <- .recycle_vp_vectors(x, n, .check_annuity_vp_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_vp_x(x)
  .check_vp_n(n)

  sapply(seq_along(x), function(j) {
    if (n[j] == 0) return(0)

    delta <- log(1 + i[j])

    f <- function(t) {
      (n[j] - t) * exp(-delta * t) * tpx(t, x = x[j], model = model, ...)
    }

    .integrate_vp(f, 0, n[j])
  })
}
