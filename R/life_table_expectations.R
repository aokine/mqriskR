# life_table_expectations.R
#
# Expectation functions for Chapter 6 of Models for Quantifying Risk.
# These functions work in the discrete tabular life-table setting.
#
# Functions:
#   ex_complete_tab()       complete expectation of life, \circ e_x
#   ex_curtate_tab()        curtate expectation of life, e_x
#   ex_temp_complete_tab()  temporary complete expectation of life, \circ e_{x:\angl{n}}
#   ex_temp_curtate_tab()   temporary curtate expectation of life, e_{x:\angl{n}}
#
# By default, complete expectations are computed under a fractional-age assumption:
#   - "udd"       : \circ e_x = e_x + 1/2 (within each year)
#   - "cf"        : constant force within each year
#   - "balducci"  : Balducci / hyperbolic within each year
#
# These functions assume the core life-table functions already exist:
#   life_table(), lx(), npx(), qx_tab(), tpx_tab()

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_n_temp <- function(n) {
  n <- as.numeric(n)
  if (any(!is.finite(n)) || any(n < 0)) {
    stop("n must be a finite nonnegative number.", call. = FALSE)
  }
  n
}

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

# -------------------------------------------------------------------------
# Curtate expectation of life: e_x
# -------------------------------------------------------------------------

#' Curtate expectation of life from a life table
#'
#' Computes the curtate expectation of life
#' \eqn{e_x = \sum_{k=1}^\infty {}_k p_x}
#' in the discrete tabular setting.
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of integer ages.
#'
#' @return Numeric vector of curtate expectations \eqn{e_x}.
#' @export
ex_curtate_tab <- function(tbl, x) {
  .validate_life_table(tbl)
  x <- as.numeric(x)

  sapply(x, function(xx) {
    lx_x <- lx(tbl, xx)
    if (is.na(lx_x) || lx_x <= 0) return(NA_real_)

    future_ages <- tbl$x[tbl$x > xx]
    if (length(future_ages) == 0L) return(0)

    sum(lx(tbl, future_ages) / lx_x, na.rm = TRUE)
  })
}

# -------------------------------------------------------------------------
# Temporary curtate expectation of life: e_{x:\angl{n}}
# -------------------------------------------------------------------------

#' Temporary curtate expectation of life from a life table
#'
#' Computes
#' \eqn{e_{x:\angl{n}} = \sum_{k=1}^{n} {}_k p_x}
#' for integer n in the discrete tabular setting.
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of integer ages.
#' @param n Numeric vector of nonnegative integers.
#'
#' @return Numeric vector of temporary curtate expectations.
#' @export
ex_temp_curtate_tab <- function(tbl, x, n) {
  .validate_life_table(tbl)
  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n

  if (any(n != floor(n))) {
    stop("For ex_temp_curtate_tab(), n must be a nonnegative integer.", call. = FALSE)
  }

  sapply(seq_along(x), function(i) {
    xx <- x[i]
    nn <- n[i]

    lx_x <- lx(tbl, xx)
    if (is.na(lx_x) || lx_x <= 0) return(NA_real_)
    if (nn == 0) return(0)

    ks <- seq_len(nn)
    probs <- npx(tbl, xx, ks)
    sum(probs, na.rm = TRUE)
  })
}

# -------------------------------------------------------------------------
# Temporary complete expectation of life: \circ e_{x:\angl{n}}
# -------------------------------------------------------------------------

#' Temporary complete expectation of life from a life table
#'
#' Computes
#' \eqn{\overset{\circ}{e}_{x:\angl{n}} = \int_0^n {}_t p_x \, dt}
#' using a within-year assumption:
#' \code{"udd"}, \code{"cf"}, or \code{"balducci"}.
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of integer ages.
#' @param n Numeric vector of nonnegative numbers.
#' @param assumption One of \code{"udd"}, \code{"cf"}, \code{"balducci"}.
#'
#' @return Numeric vector of temporary complete expectations.
#' @export
ex_temp_complete_tab <- function(tbl, x, n,
                                 assumption = c("udd", "cf", "balducci")) {
  .validate_life_table(tbl)
  assumption <- .match_assumption(assumption)

  xn <- .recycle_xn(x, n)
  x <- xn$x
  n <- xn$n
  n <- .check_n_temp(n)

  sapply(seq_along(x), function(i) {
    xx <- x[i]
    nn <- n[i]

    lx_x <- lx(tbl, xx)
    if (is.na(lx_x) || lx_x <= 0) return(NA_real_)
    if (nn == 0) return(0)

    whole <- floor(nn)
    frac  <- nn - whole

    total <- 0

    # Full years: \sum_{r=0}^{whole-1} {}_r p_x \int_0^1 {}_s p_{x+r} ds
    if (whole > 0) {
      rs <- 0:(whole - 1)
      surv_to_r <- npx(tbl, xx, rs)

      one_year_piece <- sapply(xx + rs, function(age_now) {
        integrate(
          function(s) tpx_tab(tbl, x = age_now, t = s, assumption = assumption),
          lower = 0, upper = 1, rel.tol = 1e-9
        )$value
      })

      total <- total + sum(surv_to_r * one_year_piece)
    }

    # Remaining fractional year
    if (frac > 0) {
      surv_to_whole <- npx(tbl, xx, whole)

      frac_piece <- integrate(
        function(s) tpx_tab(tbl, x = xx + whole, t = s, assumption = assumption),
        lower = 0, upper = frac, rel.tol = 1e-9
      )$value

      total <- total + surv_to_whole * frac_piece
    }

    total
  })
}

# -------------------------------------------------------------------------
# Complete expectation of life: \circ e_x
# -------------------------------------------------------------------------

#' Complete expectation of life from a life table
#'
#' Computes
#' \eqn{\overset{\circ}{e}_x = \int_0^\infty {}_t p_x \, dt}
#' using a within-year assumption:
#' \code{"udd"}, \code{"cf"}, or \code{"balducci"}.
#'
#' @param tbl A life_table object.
#' @param x Numeric vector of integer ages.
#' @param assumption One of \code{"udd"}, \code{"cf"}, \code{"balducci"}.
#'
#' @return Numeric vector of complete expectations \eqn{\overset{\circ}{e}_x}.
#' @export
ex_complete_tab <- function(tbl, x,
                            assumption = c("udd", "cf", "balducci")) {
  .validate_life_table(tbl)
  assumption <- .match_assumption(assumption)
  x <- as.numeric(x)

  sapply(x, function(xx) {
    lx_x <- lx(tbl, xx)
    if (is.na(lx_x) || lx_x <= 0) return(NA_real_)

    # Maximum duration available from age x
    max_n <- max(tbl$x, na.rm = TRUE) - xx
    if (max_n <= 0) return(0)

    ex_temp_complete_tab(tbl, x = xx, n = max_n, assumption = assumption)
  })
}
