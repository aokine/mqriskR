#' m-thly contingent annuity functions (Chapter 8)
#'
#' Chapter 8 functions for life annuities payable m-thly.
#'
#' These functions implement the exact m-thly formulas from Section 8.5
#' using the survival model functions already available in the package.
#'
#' In all cases, the annuity is one unit per year payable in m equal
#' installments, so each payment is of size \eqn{1/m}.
#'
#' @name annuity_mthly
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_i_m <- function(i, m) {
  i <- as.numeric(i)
  m <- as.numeric(m)

  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite number greater than -1.", call. = FALSE)
  }
  if (length(m) != 1L || !is.finite(m) || m <= 0 || abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }

  list(i = i, m = as.integer(round(m)))
}

#' @noRd
.check_x_scalar_m <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stop("x must be a single nonnegative finite number.", call. = FALSE)
  }
  x
}

#' @noRd
.check_n_scalar_m <- function(n) {
  n <- as.numeric(n)
  if (length(n) != 1L || !is.finite(n) || n < 0) {
    stop("n must be a single nonnegative finite number.", call. = FALSE)
  }
  n
}

#' @noRd
.recycle_x_n_m <- function(x, n) {
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
.max_subperiod_from_model <- function(x, m, model, ...) {
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

#' @noRd
.sum_mthly_annuity_terms <- function(x, i, m, model, ..., due = FALSE,
                                     temporary = NULL, k_max = 200000,
                                     tol = 1e-12) {
  x <- .check_x_scalar_m(x)
  chk <- .check_i_m(i, m)
  i <- chk$i
  m <- chk$m

  if (!is.null(temporary)) {
    temporary <- .check_n_scalar_m(temporary)
  }

  j <- (1 + i)^(1 / m) - 1
  vsub <- 1 / (1 + j)

  max_sub <- .max_subperiod_from_model(x, m, model, ...)

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
#' Computes the exact whole life m-thly annuity-immediate,
#' with annual payment rate 1 split into m equal payments of size 1/m.
#'
#' \deqn{a_x^{(m)} = \frac{1}{m}\sum_{t=1}^{\infty} v^{t/m}\,{}_{t/m}p_x}
#'
#' @name annuity_mthly_whole_immediate
#' @aliases ax_m
#' @param x Age.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param model Survival model name.
#' @param ... Additional parameters passed to the survival model.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return Numeric vector of annuity values.
#' @examples
#' ax_m(40, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
ax_m <- function(x, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  x <- as.numeric(x)

  sapply(x, function(xx) {
    .sum_mthly_annuity_terms(
      x = xx, i = i, m = m, model = model, ...,
      due = FALSE, temporary = NULL, k_max = k_max, tol = tol
    )
  })
}

#' Whole life m-thly annuity-due
#'
#' Computes the exact whole life m-thly annuity-due.
#'
#' \deqn{\ddot{a}_x^{(m)} = \frac{1}{m}\sum_{t=0}^{\infty} v^{t/m}\,{}_{t/m}p_x}
#'
#' @name annuity_mthly_whole_due
#' @aliases adotx_m
#' @inheritParams annuity_mthly_whole_immediate
#'
#' @return Numeric vector of annuity values.
#' @examples
#' adotx_m(40, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
adotx_m <- function(x, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  x <- as.numeric(x)

  sapply(x, function(xx) {
    .sum_mthly_annuity_terms(
      x = xx, i = i, m = m, model = model, ...,
      due = TRUE, temporary = NULL, k_max = k_max, tol = tol
    )
  })
}

# -------------------------------------------------------------------------
# Temporary m-thly annuities
# -------------------------------------------------------------------------

#' Temporary m-thly annuity-immediate
#'
#' Computes the exact temporary m-thly annuity-immediate.
#'
#' \deqn{a_{x:\angl{n}}^{(m)} = \frac{1}{m}\sum_{t=1}^{mn} v^{t/m}\,{}_{t/m}p_x}
#'
#' @name annuity_mthly_temp_immediate
#' @aliases axn_m
#' @param x Age.
#' @param n Term in years.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param model Survival model name.
#' @param ... Additional parameters passed to the survival model.
#'
#' @return Numeric vector of annuity values.
#' @examples
#' axn_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
axn_m <- function(x, n, m, i, model, ...) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    .sum_mthly_annuity_terms(
      x = x[j], i = i, m = m, model = model, ...,
      due = FALSE, temporary = n[j]
    )
  })
}

#' Temporary m-thly annuity-due
#'
#' Computes the exact temporary m-thly annuity-due.
#'
#' \deqn{\ddot{a}_{x:\angl{n}}^{(m)} = \frac{1}{m}\sum_{t=0}^{mn-1} v^{t/m}\,{}_{t/m}p_x}
#'
#' @name annuity_mthly_temp_due
#' @aliases adotxn_m
#' @inheritParams annuity_mthly_temp_immediate
#'
#' @return Numeric vector of annuity values.
#' @examples
#' adotxn_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
adotxn_m <- function(x, n, m, i, model, ...) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    .sum_mthly_annuity_terms(
      x = x[j], i = i, m = m, model = model, ...,
      due = TRUE, temporary = n[j]
    )
  })
}

# -------------------------------------------------------------------------
# Deferred m-thly annuities
# -------------------------------------------------------------------------

#' Deferred whole life m-thly annuity-immediate
#'
#' Computes the deferred whole life m-thly annuity-immediate using
#' \deqn{{}_nE_x \, a_{x+n}^{(m)}}
#'
#' @name annuity_mthly_deferred_immediate
#' @aliases nax_m
#' @param x Age.
#' @param n Deferral period in years.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param model Survival model name.
#' @param ... Additional parameters passed to the survival model.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return Numeric vector of annuity values.
#' @examples
#' nax_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
nax_m <- function(x, n, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nEx(x[j], n[j], i, model = model, ...) *
      ax_m(x[j] + n[j], m = m, i = i, model = model, ..., k_max = k_max, tol = tol)
  })
}

#' Deferred whole life m-thly annuity-due
#'
#' Computes the deferred whole life m-thly annuity-due using
#' \deqn{{}_nE_x \, \ddot{a}_{x+n}^{(m)}}
#'
#' @name annuity_mthly_deferred_due
#' @aliases nadotx_m
#' @inheritParams annuity_mthly_deferred_immediate
#'
#' @return Numeric vector of annuity values.
#' @examples
#' nadotx_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
nadotx_m <- function(x, n, m, i, model, ..., k_max = 200000, tol = 1e-12) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    nEx(x[j], n[j], i, model = model, ...) *
      adotx_m(x[j] + n[j], m = m, i = i, model = model, ..., k_max = k_max, tol = tol)
  })
}

# -------------------------------------------------------------------------
# Actuarial accumulated values
# -------------------------------------------------------------------------

#' Temporary m-thly annuity-immediate actuarial accumulated value
#'
#' Computes
#' \deqn{s_{x:\angl{n}}^{(m)} = a_{x:\angl{n}}^{(m)} / {}_nE_x}
#'
#' @name annuity_mthly_accum_immediate
#' @aliases sxn_m
#' @inheritParams annuity_mthly_temp_immediate
#'
#' @return Numeric vector of actuarial accumulated values.
#' @examples
#' sxn_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
sxn_m <- function(x, n, m, i, model, ...) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    ex <- nEx(x[j], n[j], i, model = model, ...)
    if (ex <= 0) return(Inf)
    axn_m(x[j], n[j], m = m, i = i, model = model, ...) / ex
  })
}

#' Temporary m-thly annuity-due actuarial accumulated value
#'
#' Computes
#' \deqn{\ddot{s}_{x:\angl{n}}^{(m)} = \ddot{a}_{x:\angl{n}}^{(m)} / {}_nE_x}
#'
#' @name annuity_mthly_accum_due
#' @aliases sdotxn_m
#' @inheritParams annuity_mthly_temp_immediate
#'
#' @return Numeric vector of actuarial accumulated values.
#' @examples
#' sdotxn_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
sdotxn_m <- function(x, n, m, i, model, ...) {
  chk <- .check_i_m(i, m)
  m <- chk$m
  xn <- .recycle_x_n_m(x, n)
  x <- xn$x
  n <- xn$n

  sapply(seq_along(x), function(j) {
    ex <- nEx(x[j], n[j], i, model = model, ...)
    if (ex <= 0) return(Inf)
    adotxn_m(x[j], n[j], m = m, i = i, model = model, ...) / ex
  })
}
