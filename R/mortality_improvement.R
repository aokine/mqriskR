#' Mortality improvement projection functions (Chapter 8)
#'
#' Functions for Chapter 8 Section 8.8 mortality improvement projection.
#'
#' These functions project one-year death and survival probabilities forward
#' from a base year and evaluate annuity values under projected mortality.
#'
#' The standard projection used is
#' \deqn{q_x^{[Y]} = q_x^{[B]} (1 - AA_x)^{Y-B}}
#' where:
#' \itemize{
#'   \item \eqn{B} is the base year,
#'   \item \eqn{Y} is the projection year,
#'   \item \eqn{AA_x} is the mortality improvement factor at age \eqn{x}.
#' }
#'
#' @name mortality_improvement
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_scalar_year <- function(y, name) {
  y <- as.numeric(y)
  if (length(y) != 1L || !is.finite(y) || abs(y - round(y)) > 1e-10) {
    stop(paste0(name, " must be a single integer-like year."), call. = FALSE)
  }
  as.integer(round(y))
}

#' @noRd
.check_year_vector <- function(y, name) {
  y <- as.numeric(y)
  if (length(y) == 0L) {
    stop(paste0(name, " must have positive length."), call. = FALSE)
  }
  if (any(!is.finite(y)) || any(abs(y - round(y)) > 1e-10)) {
    stop(paste0(name, " must contain integer-like years."), call. = FALSE)
  }
  as.integer(round(y))
}

#' @noRd
.check_prob_vector <- function(x, name) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop(paste0(name, " must be finite."), call. = FALSE)
  }
  if (any(x < 0 | x > 1)) {
    stop(paste0(name, " must lie in [0, 1]."), call. = FALSE)
  }
  x
}

#' @noRd
.recycle_proj_inputs <- function(qx_base, AAx, proj_year) {
  qx_base <- as.numeric(qx_base)
  AAx <- as.numeric(AAx)
  proj_year <- as.numeric(proj_year)

  lens <- c(length(qx_base), length(AAx), length(proj_year))
  L <- max(lens)

  if (!all(lens %in% c(1L, L))) {
    stop("qx_base, AAx, and proj_year must have the same length, or length 1.", call. = FALSE)
  }

  if (length(qx_base) == 1L) qx_base <- rep(qx_base, L)
  if (length(AAx) == 1L) AAx <- rep(AAx, L)
  if (length(proj_year) == 1L) proj_year <- rep(proj_year, L)

  list(qx_base = qx_base, AAx = AAx, proj_year = proj_year)
}

# -------------------------------------------------------------------------
# One-year projected probabilities
# -------------------------------------------------------------------------

#' Project one-year death probability under mortality improvement
#'
#' Computes projected one-year death probability
#' \deqn{q_x^{[Y]} = q_x^{[B]} (1-AA_x)^{Y-B}}.
#'
#' @param qx_base Base-year one-year death probability \eqn{q_x^{[B]}}.
#' @param AAx Mortality improvement factor \eqn{AA_x}.
#' @param base_year Base year \eqn{B}.
#' @param proj_year Projection year \eqn{Y}. May be scalar or vector.
#'
#' @return Numeric vector of projected one-year death probabilities.
#' @export
qx_proj <- function(qx_base, AAx, base_year, proj_year) {
  qx_base <- .check_prob_vector(qx_base, "qx_base")
  AAx <- .check_prob_vector(AAx, "AAx")
  base_year <- .check_scalar_year(base_year, "base_year")
  proj_year <- .check_year_vector(proj_year, "proj_year")

  rr <- .recycle_proj_inputs(qx_base, AAx, proj_year)
  qx_base <- rr$qx_base
  AAx <- rr$AAx
  proj_year <- rr$proj_year

  qx_base * (1 - AAx)^(proj_year - base_year)
}

#' Project one-year survival probability under mortality improvement
#'
#' Computes projected one-year survival probability
#' \deqn{p_x^{[Y]} = 1 - q_x^{[Y]}}.
#'
#' @inheritParams qx_proj
#'
#' @return Numeric vector of projected one-year survival probabilities.
#' @export
px_proj <- function(qx_base, AAx, base_year, proj_year) {
  1 - qx_proj(
    qx_base = qx_base,
    AAx = AAx,
    base_year = base_year,
    proj_year = proj_year
  )
}

# -------------------------------------------------------------------------
# Multi-year survival under improvement
# -------------------------------------------------------------------------

#' Multi-year survival probability under mortality improvement
#'
#' Computes survival over \code{n} years starting at age \code{x0} in year
#' \code{issue_year}, using base-year mortality \code{qx_base_vec} and
#' improvement factors \code{AAx_vec}.
#'
#' The vectors should correspond to ages
#' \code{x0, x0+1, ..., x0+n-1}.
#'
#' The survival probability is
#' \deqn{\prod_{j=0}^{n-1} \left(1 - q_{x+j}^{[issue\_year+j]}\right)}.
#'
#' @param x0 Issue age.
#' @param n Number of years.
#' @param qx_base_vec Base-year one-year death probabilities for successive ages.
#' @param AAx_vec Mortality improvement factors for successive ages.
#' @param base_year Base year.
#' @param issue_year Issue year.
#'
#' @return Numeric scalar.
#' @export
tpx_improved <- function(x0, n, qx_base_vec, AAx_vec, base_year, issue_year) {
  x0 <- as.numeric(x0)
  n <- as.numeric(n)

  if (length(x0) != 1L || !is.finite(x0) || x0 < 0) {
    stop("x0 must be a single nonnegative number.", call. = FALSE)
  }
  if (length(n) != 1L || !is.finite(n) || n < 0 || abs(n - round(n)) > 1e-10) {
    stop("n must be a single nonnegative integer.", call. = FALSE)
  }

  n <- as.integer(round(n))
  base_year <- .check_scalar_year(base_year, "base_year")
  issue_year <- .check_scalar_year(issue_year, "issue_year")

  qx_base_vec <- .check_prob_vector(qx_base_vec, "qx_base_vec")
  AAx_vec <- .check_prob_vector(AAx_vec, "AAx_vec")

  if (n == 0L) return(1)

  if (length(qx_base_vec) < n || length(AAx_vec) < n) {
    stop("qx_base_vec and AAx_vec must each have length at least n.", call. = FALSE)
  }

  years <- issue_year + 0:(n - 1L)

  qproj <- qx_proj(
    qx_base = qx_base_vec[1:n],
    AAx = AAx_vec[1:n],
    base_year = base_year,
    proj_year = years
  )

  prod(1 - qproj)
}

# -------------------------------------------------------------------------
# Chapter 8 annuity functions under improvement
# -------------------------------------------------------------------------

#' Temporary annuity-immediate under mortality improvement
#'
#' Computes
#' \deqn{a_{x:\angl{n}} = \sum_{t=1}^n v^t \, {}_tp_x^{(\mathrm{improved})}}
#' using projected mortality rates.
#'
#' @param x0 Issue age.
#' @param n Term in years.
#' @param i Effective annual interest rate.
#' @param qx_base_vec Base-year one-year death probabilities for ages
#'   \code{x0, x0+1, ..., x0+n-1}.
#' @param AAx_vec Mortality improvement factors for ages
#'   \code{x0, x0+1, ..., x0+n-1}.
#' @param base_year Base year.
#' @param issue_year Issue year.
#'
#' @return Numeric scalar.
#' @export
axn_improved <- function(x0, n, i, qx_base_vec, AAx_vec, base_year, issue_year) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value > -1.", call. = FALSE)
  }

  n <- as.numeric(n)
  if (length(n) != 1L || !is.finite(n) || n < 0 || abs(n - round(n)) > 1e-10) {
    stop("n must be a single nonnegative integer.", call. = FALSE)
  }
  n <- as.integer(round(n))

  if (n == 0L) return(0)

  survs <- sapply(1:n, function(t) {
    tpx_improved(
      x0 = x0,
      n = t,
      qx_base_vec = qx_base_vec,
      AAx_vec = AAx_vec,
      base_year = base_year,
      issue_year = issue_year
    )
  })

  sum(discount(i, 1:n) * survs)
}

#' Deferred temporary annuity-immediate under mortality improvement
#'
#' Computes
#' \deqn{{}_{u\mid}a_{x:\angl{n}} = \sum_{t=u+1}^{u+n} v^t \, {}_tp_x^{(\mathrm{improved})}}
#'
#' @param x0 Issue age.
#' @param u Deferral period in years.
#' @param n Temporary payment period in years.
#' @param i Effective annual interest rate.
#' @param qx_base_vec Base-year one-year death probabilities for ages
#'   \code{x0, x0+1, ..., x0+u+n-1}.
#' @param AAx_vec Mortality improvement factors for ages
#'   \code{x0, x0+1, ..., x0+u+n-1}.
#' @param base_year Base year.
#' @param issue_year Issue year.
#'
#' @return Numeric scalar.
#' @export
naxn_improved <- function(x0, u, n, i, qx_base_vec, AAx_vec, base_year, issue_year) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value > -1.", call. = FALSE)
  }

  u <- as.numeric(u)
  n <- as.numeric(n)

  if (length(u) != 1L || !is.finite(u) || u < 0 || abs(u - round(u)) > 1e-10) {
    stop("u must be a single nonnegative integer.", call. = FALSE)
  }
  if (length(n) != 1L || !is.finite(n) || n < 0 || abs(n - round(n)) > 1e-10) {
    stop("n must be a single nonnegative integer.", call. = FALSE)
  }

  u <- as.integer(round(u))
  n <- as.integer(round(n))

  if (n == 0L) return(0)

  times <- (u + 1L):(u + n)

  survs <- sapply(times, function(t) {
    tpx_improved(
      x0 = x0,
      n = t,
      qx_base_vec = qx_base_vec,
      AAx_vec = AAx_vec,
      base_year = base_year,
      issue_year = issue_year
    )
  })

  sum(discount(i, times) * survs)
}

#' Whole life annuity-immediate under mortality improvement
#'
#' Computes a truncated whole life annuity-immediate under projected mortality.
#'
#' @param x0 Issue age.
#' @param i Effective annual interest rate.
#' @param qx_base_vec Base-year one-year death probabilities for successive ages.
#' @param AAx_vec Mortality improvement factors for successive ages.
#' @param base_year Base year.
#' @param issue_year Issue year.
#'
#' @return Numeric scalar.
#' @export
ax_improved <- function(x0, i, qx_base_vec, AAx_vec, base_year, issue_year) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value > -1.", call. = FALSE)
  }

  qx_base_vec <- .check_prob_vector(qx_base_vec, "qx_base_vec")
  AAx_vec <- .check_prob_vector(AAx_vec, "AAx_vec")

  L <- min(length(qx_base_vec), length(AAx_vec))
  if (L == 0L) return(0)

  survs <- sapply(1:L, function(t) {
    tpx_improved(
      x0 = x0,
      n = t,
      qx_base_vec = qx_base_vec[1:L],
      AAx_vec = AAx_vec[1:L],
      base_year = base_year,
      issue_year = issue_year
    )
  })

  sum(discount(i, 1:L) * survs)
}
