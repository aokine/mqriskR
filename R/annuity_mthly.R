#' m-thly contingent annuity functions
#'
#' Life annuities payable m-thly.
#'
#' @param x Age.
#' @param n Term or deferral period in years.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param model Survival model name.
#' @param ... Additional parameters passed to the survival model.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return
#' Numeric vector containing the requested m-thly annuity value or actuarial
#' accumulated value.
#'
#' @name annuity_mthly
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_ann_m <- function(m) {
  m <- as.numeric(m)

  if (length(m) != 1L || !is.finite(m) || m <= 0 ||
      abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }

  as.integer(round(m))
}

.check_i_ann_m <- function(i) {
  i <- as.numeric(i)

  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }

  i
}

.recycle_ann_m_vectors <- function(...) {
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

  lapply(args, function(z) if (length(z) == 1L) rep(z, common_len) else z)
}

.check_x_ann_m <- function(x) {
  if (any(x < 0)) stop("x must be nonnegative.", call. = FALSE)
  invisible(TRUE)
}

.check_n_ann_m <- function(n) {
  if (any(n < 0)) stop("n must be nonnegative.", call. = FALSE)
  invisible(TRUE)
}

.max_subperiod_ann_m <- function(x, m, model, ...) {
  dots <- list(...)
  model <- tolower(model)

  if (model == "uniform") {
    omega <- dots$omega

    if (is.null(omega)) {
      stop("For model = 'uniform', omega must be supplied.", call. = FALSE)
    }

    return(max(0L, floor((omega - x) * m)))
  }

  Inf
}

.sum_mthly_annuity_terms <- function(x, i, m, model, ..., due = FALSE,
                                     temporary = NULL, k_max = 200000,
                                     tol = 1e-12) {
  x <- as.numeric(x)
  i <- as.numeric(i)

  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stop("x must be a single nonnegative finite number.", call. = FALSE)
  }

  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite value greater than -1.", call. = FALSE)
  }

  if (!is.null(temporary)) {
    temporary <- as.numeric(temporary)

    if (length(temporary) != 1L || !is.finite(temporary) || temporary < 0) {
      stop("temporary must be a single nonnegative finite number.", call. = FALSE)
    }
  }

  j <- (1 + i)^(1 / m) - 1
  vsub <- 1 / (1 + j)

  max_sub <- .max_subperiod_ann_m(x, m, model, ...)

  if (due) {
    if (!is.null(temporary)) {
      n_pay <- round(temporary * m)

      if (abs(temporary * m - n_pay) > 1e-10) {
        stop("n * m must be an integer for m-thly annuity functions.", call. = FALSE)
      }

      if (n_pay == 0L) return(0)
      idx <- 0:(n_pay - 1L)

    } else if (is.finite(max_sub)) {
      idx <- 0:max_sub
    } else {
      idx <- 0:k_max
    }
  } else {
    if (!is.null(temporary)) {
      n_pay <- round(temporary * m)

      if (abs(temporary * m - n_pay) > 1e-10) {
        stop("n * m must be an integer for m-thly annuity functions.", call. = FALSE)
      }

      if (n_pay == 0L) return(0)
      idx <- 1:n_pay

    } else if (is.finite(max_sub)) {
      if (max_sub <= 0L) return(0)
      idx <- 1:max_sub
    } else {
      idx <- 1:k_max
    }
  }

  times <- idx / m
  terms <- (1 / m) * (vsub^idx) * tpx(times, x = x, model = model, ...)

  terms[!is.finite(terms)] <- 0

  if (!is.null(temporary) || is.finite(max_sub)) {
    return(sum(terms))
  }

  acc <- 0
  small_run <- 0L

  for (jj in seq_along(terms)) {
    acc <- acc + terms[jj]

    if (terms[jj] < tol) {
      small_run <- small_run + 1L
      if (small_run >= 20L) break
    } else {
      small_run <- 0L
    }
  }

  acc
}

# -------------------------------------------------------------------------
# Whole life m-thly annuities
# -------------------------------------------------------------------------

#' Whole life m-thly annuity-immediate
#'
#' @rdname annuity_mthly
#' @export
ax_m <- function(x, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, .check_i_ann_m(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_x_ann_m(x)

  sapply(seq_along(x), function(j) {
    .sum_mthly_annuity_terms(
      x = x[j], i = i[j], m = m, model = model, ...,
      due = FALSE, temporary = NULL, k_max = k_max, tol = tol
    )
  })
}

#' Whole life m-thly annuity-due
#'
#' @rdname annuity_mthly
#' @export
adotx_m <- function(x, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, .check_i_ann_m(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_x_ann_m(x)

  sapply(seq_along(x), function(j) {
    .sum_mthly_annuity_terms(
      x = x[j], i = i[j], m = m, model = model, ...,
      due = TRUE, temporary = NULL, k_max = k_max, tol = tol
    )
  })
}

# -------------------------------------------------------------------------
# Temporary m-thly annuities
# -------------------------------------------------------------------------

#' Temporary m-thly annuity-immediate
#'
#' @rdname annuity_mthly
#' @export
axn_m <- function(x, n, m, i, model, ...) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, n, .check_i_ann_m(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_ann_m(x)
  .check_n_ann_m(n)

  sapply(seq_along(x), function(j) {
    .sum_mthly_annuity_terms(
      x = x[j], i = i[j], m = m, model = model, ...,
      due = FALSE, temporary = n[j]
    )
  })
}

#' Temporary m-thly annuity-due
#'
#' @rdname annuity_mthly
#' @export
adotxn_m <- function(x, n, m, i, model, ...) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, n, .check_i_ann_m(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_ann_m(x)
  .check_n_ann_m(n)

  sapply(seq_along(x), function(j) {
    .sum_mthly_annuity_terms(
      x = x[j], i = i[j], m = m, model = model, ...,
      due = TRUE, temporary = n[j]
    )
  })
}

# -------------------------------------------------------------------------
# Deferred m-thly annuities
# -------------------------------------------------------------------------

#' Deferred whole life m-thly annuity-immediate
#'
#' @rdname annuity_mthly
#' @export
nax_m <- function(x, n, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, n, .check_i_ann_m(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_ann_m(x)
  .check_n_ann_m(n)

  sapply(seq_along(x), function(j) {
    nEx(x[j], n[j], i[j], model = model, ...) *
      ax_m(
        x[j] + n[j],
        m = m,
        i = i[j],
        model = model,
        ...,
        k_max = k_max,
        tol = tol
      )
  })
}

#' Deferred whole life m-thly annuity-due
#'
#' @rdname annuity_mthly
#' @export
nadotx_m <- function(x, n, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, n, .check_i_ann_m(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_ann_m(x)
  .check_n_ann_m(n)

  sapply(seq_along(x), function(j) {
    nEx(x[j], n[j], i[j], model = model, ...) *
      adotx_m(
        x[j] + n[j],
        m = m,
        i = i[j],
        model = model,
        ...,
        k_max = k_max,
        tol = tol
      )
  })
}

# -------------------------------------------------------------------------
# Actuarial accumulated values
# -------------------------------------------------------------------------

#' Temporary m-thly annuity-immediate actuarial accumulated value
#'
#' @rdname annuity_mthly
#' @export
sxn_m <- function(x, n, m, i, model, ...) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, n, .check_i_ann_m(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_ann_m(x)
  .check_n_ann_m(n)

  sapply(seq_along(x), function(j) {
    ex <- nEx(x[j], n[j], i[j], model = model, ...)

    if (is.na(ex) || ex <= 0) return(Inf)

    axn_m(x[j], n[j], m = m, i = i[j], model = model, ...) / ex
  })
}

#' Temporary m-thly annuity-due actuarial accumulated value
#'
#' @rdname annuity_mthly
#' @export
sdotxn_m <- function(x, n, m, i, model, ...) {
  m <- .check_ann_m(m)

  args <- .recycle_ann_m_vectors(x, n, .check_i_ann_m(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_ann_m(x)
  .check_n_ann_m(n)

  sapply(seq_along(x), function(j) {
    ex <- nEx(x[j], n[j], i[j], model = model, ...)

    if (is.na(ex) || ex <= 0) return(Inf)

    adotxn_m(x[j], n[j], m = m, i = i[j], model = model, ...) / ex
  })
}
