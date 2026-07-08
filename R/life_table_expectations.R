# life_table_expectations.R
#
# Expectation functions for discrete tabular life tables.

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_n_temp <- function(n) {
  n <- as.numeric(n)

  if (length(n) == 0L || any(!is.finite(n)) || any(n < 0)) {
    stop("n must contain finite nonnegative values.", call. = FALSE)
  }

  n
}

.recycle_expectation_args <- function(...) {
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

.check_age_values <- function(x) {
  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Curtate expectation of life
# -------------------------------------------------------------------------

#' Curtate expectation of life from a life table
#'
#' Computes the curtate expectation of life
#' \eqn{e_x = \sum_{k=1}^{\infty} {}_k p_x}
#' in the discrete tabular setting.
#'
#' @param tbl A \code{life_table} object.
#' @param x Numeric vector of integer ages.
#'
#' @return Numeric vector of curtate expectations \eqn{e_x}.
#' @export
ex_curtate_tab <- function(tbl, x) {
  .validate_life_table(tbl)

  args <- .recycle_expectation_args(x)
  x <- args[[1L]]
  .check_age_values(x)

  sapply(x, function(xx) {
    lx_x <- lx(tbl, xx)

    if (is.na(lx_x) || lx_x <= 0) {
      return(NA_real_)
    }

    max_k <- max(tbl$x, na.rm = TRUE) - xx

    if (max_k <= 0) {
      return(0)
    }

    ks <- seq_len(max_k)
    probs <- npx(tbl, x = xx, n = ks)

    sum(probs, na.rm = TRUE)
  })
}

# -------------------------------------------------------------------------
# Temporary curtate expectation of life
# -------------------------------------------------------------------------

#' Temporary curtate expectation of life from a life table
#'
#' Computes
#' \eqn{e_{x:\overline{n}|} = \sum_{k=1}^{n} {}_k p_x}
#' for integer \eqn{n} in the discrete tabular setting.
#'
#' @param tbl A \code{life_table} object.
#' @param x Numeric vector of integer ages.
#' @param n Numeric vector of nonnegative integers.
#'
#' @return Numeric vector of temporary curtate expectations.
#' @export
ex_temp_curtate_tab <- function(tbl, x, n) {
  .validate_life_table(tbl)

  args <- .recycle_expectation_args(x, n)
  x <- args[[1L]]
  n <- args[[2L]]

  .check_age_values(x)

  if (any(n < 0 | n != floor(n))) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]

    lx_x <- lx(tbl, xx)

    if (is.na(lx_x) || lx_x <= 0) {
      return(NA_real_)
    }

    if (nn == 0) {
      return(0)
    }

    ks <- seq_len(nn)
    probs <- npx(tbl, x = xx, n = ks)

    sum(probs, na.rm = TRUE)
  })
}

# -------------------------------------------------------------------------
# Temporary complete expectation of life
# -------------------------------------------------------------------------

#' Temporary complete expectation of life from a life table
#'
#' Computes
#' \eqn{\overset{\circ}{e}_{x:\overline{n}|} =
#' \int_0^n {}_t p_x \, dt}
#' using a within-year assumption:
#' \code{"udd"}, \code{"cf"}, or \code{"balducci"}.
#'
#' @param tbl A \code{life_table} object.
#' @param x Numeric vector of integer ages.
#' @param n Numeric vector of nonnegative durations.
#' @param assumption One of \code{"udd"}, \code{"cf"}, \code{"balducci"}.
#'
#' @return Numeric vector of temporary complete expectations.
#' @export
ex_temp_complete_tab <- function(tbl, x, n,
                                 assumption = c("udd", "cf", "balducci")) {
  .validate_life_table(tbl)

  assumption <- .match_assumption(assumption)

  args <- .recycle_expectation_args(x, n)
  x <- args[[1L]]
  n <- .check_n_temp(args[[2L]])

  .check_age_values(x)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]

    lx_x <- lx(tbl, xx)

    if (is.na(lx_x) || lx_x <= 0) {
      return(NA_real_)
    }

    if (nn == 0) {
      return(0)
    }

    whole <- floor(nn)
    frac <- nn - whole
    total <- 0

    if (whole > 0) {
      rs <- 0:(whole - 1L)
      surv_to_r <- npx(tbl, x = xx, n = rs)

      one_year_piece <- vapply(rs, function(r) {
        age_now <- xx + r
        lx_now <- lx(tbl, age_now)

        if (is.na(lx_now) || lx_now <= 0) {
          return(0)
        }

        integrate(
          function(s) {
            tpx_tab(tbl, x = age_now, t = s, assumption = assumption)
          },
          lower = 0,
          upper = 1,
          rel.tol = 1e-9
        )$value
      }, numeric(1))

      total <- total + sum(surv_to_r * one_year_piece, na.rm = TRUE)
    }

    if (frac > 0) {
      surv_to_whole <- npx(tbl, x = xx, n = whole)
      age_now <- xx + whole
      lx_now <- lx(tbl, age_now)

      if (!is.na(surv_to_whole) &&
          !is.na(lx_now) &&
          surv_to_whole > 0 &&
          lx_now > 0) {
        frac_piece <- integrate(
          function(s) {
            tpx_tab(tbl, x = age_now, t = s, assumption = assumption)
          },
          lower = 0,
          upper = frac,
          rel.tol = 1e-9
        )$value

        total <- total + surv_to_whole * frac_piece
      }
    }

    total
  })
}

# -------------------------------------------------------------------------
# Complete expectation of life
# -------------------------------------------------------------------------

#' Complete expectation of life from a life table
#'
#' Computes
#' \eqn{\overset{\circ}{e}_x = \int_0^\infty {}_t p_x \, dt}
#' using a within-year assumption:
#' \code{"udd"}, \code{"cf"}, or \code{"balducci"}.
#'
#' @param tbl A \code{life_table} object.
#' @param x Numeric vector of integer ages.
#' @param assumption One of \code{"udd"}, \code{"cf"}, \code{"balducci"}.
#'
#' @return Numeric vector of complete expectations \eqn{\overset{\circ}{e}_x}.
#' @export
ex_complete_tab <- function(tbl, x,
                            assumption = c("udd", "cf", "balducci")) {
  .validate_life_table(tbl)

  assumption <- .match_assumption(assumption)

  args <- .recycle_expectation_args(x)
  x <- args[[1L]]

  .check_age_values(x)

  sapply(x, function(xx) {
    lx_x <- lx(tbl, xx)

    if (is.na(lx_x) || lx_x <= 0) {
      return(NA_real_)
    }

    max_n <- max(tbl$x, na.rm = TRUE) - xx

    if (max_n <= 0) {
      return(0)
    }

    ex_temp_complete_tab(tbl, x = xx, n = max_n, assumption = assumption)
  })
}
