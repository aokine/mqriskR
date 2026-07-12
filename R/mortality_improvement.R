#' Mortality improvement projection functions
#'
#' Functions for projecting one-year death and survival probabilities under
#' mortality improvement and evaluating annuity values under projected
#' mortality.
#'
#' The standard projection used is
#' \deqn{q_x^{[Y]} = q_x^{[B]} (1 - AA_x)^{Y-B},}
#' where \eqn{B} is the base year, \eqn{Y} is the projection year, and
#' \eqn{AA_x} is the mortality improvement factor at age \eqn{x}.
#'
#' @name mortality_improvement_projection
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_scalar_year <- function(y, name) {
  y <- as.numeric(y)
  if (length(y) != 1L || !is.finite(y) || abs(y - round(y)) > 1e-10) {
    stop(paste0(name, " must be a single integer-like year."), call. = FALSE)
  }
  as.integer(round(y))
}

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

.check_prob_vector <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) == 0L) {
    stop(paste0(name, " must have positive length."), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(paste0(name, " must be finite."), call. = FALSE)
  }
  if (any(x < 0 | x > 1)) {
    stop(paste0(name, " must lie in [0, 1]."), call. = FALSE)
  }
  x
}

.recycle_proj_inputs <- function(qx_base, AAx, proj_year) {
  qx_base <- as.numeric(qx_base)
  AAx <- as.numeric(AAx)
  proj_year <- as.numeric(proj_year)

  lens <- c(length(qx_base), length(AAx), length(proj_year))
  L <- max(lens)

  if (!all(lens %in% c(1L, L))) {
    stop("qx_base, AAx, and proj_year must have compatible lengths.", call. = FALSE)
  }

  list(
    qx_base = rep_len(qx_base, L),
    AAx = rep_len(AAx, L),
    proj_year = rep_len(proj_year, L)
  )
}

.recycle_improved_annuity_inputs <- function(x0, n = NULL, u = NULL, i) {
  args <- list(x0 = as.numeric(x0), i = as.numeric(i))
  if (!is.null(n)) args$n <- as.numeric(n)
  if (!is.null(u)) args$u <- as.numeric(u)

  lens <- vapply(args, length, integer(1))
  L <- max(lens)

  if (any(lens == 0L)) {
    stop("Vectorized inputs must have positive length.", call. = FALSE)
  }
  if (!all(lens %in% c(1L, L))) {
    stop("Vectorized inputs must have compatible lengths.", call. = FALSE)
  }
  if (any(vapply(args, function(z) any(!is.finite(z)), logical(1)))) {
    stop("Vectorized inputs must be finite.", call. = FALSE)
  }

  lapply(args, rep_len, length.out = L)
}

.check_nonnegative_integer_scalar <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x < 0 || abs(x - round(x)) > 1e-10) {
    stop(paste0(name, " must be a single nonnegative integer."), call. = FALSE)
  }
  as.integer(round(x))
}

.check_improved_interest <- function(i) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value > -1.", call. = FALSE)
  }
  i
}

# -------------------------------------------------------------------------
# One-year projected probabilities
# -------------------------------------------------------------------------

#' Project one-year death probability under mortality improvement
#'
#' @param qx_base Base-year one-year death probability.
#' @param AAx Mortality improvement factor.
#' @param base_year Base year.
#' @param proj_year Projection year. May be scalar or vector.
#'
#' @return Numeric vector of projected one-year death probabilities.
#' @rdname mortality_improvement_projection
#' @export
qx_proj <- function(qx_base, AAx, base_year, proj_year) {
  qx_base <- .check_prob_vector(qx_base, "qx_base")
  AAx <- .check_prob_vector(AAx, "AAx")
  base_year <- .check_scalar_year(base_year, "base_year")
  proj_year <- .check_year_vector(proj_year, "proj_year")

  rr <- .recycle_proj_inputs(qx_base, AAx, proj_year)

  rr$qx_base * (1 - rr$AAx)^(rr$proj_year - base_year)
}

#' Project one-year survival probability under mortality improvement
#'
#' @rdname mortality_improvement_projection
#' @export
px_proj <- function(qx_base, AAx, base_year, proj_year) {
  1 - qx_proj(qx_base, AAx, base_year, proj_year)
}

# -------------------------------------------------------------------------
# Multi-year survival under improvement
# -------------------------------------------------------------------------

#' Multi-year survival probability under mortality improvement
#'
#' @param x0 Issue age.
#' @param n Number of years.
#' @param qx_base_vec Base-year death probabilities for successive ages.
#' @param AAx_vec Mortality improvement factors for successive ages.
#' @param issue_year Issue year.
#'
#' @rdname mortality_improvement_projection
#' @export
tpx_improved <- function(x0, n, qx_base_vec, AAx_vec, base_year, issue_year) {
  x0 <- as.numeric(x0)
  if (length(x0) != 1L || !is.finite(x0) || x0 < 0) {
    stop("x0 must be a single nonnegative number.", call. = FALSE)
  }

  n <- .check_nonnegative_integer_scalar(n, "n")
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
# Annuity functions under improvement
# -------------------------------------------------------------------------

#' Temporary annuity-immediate under mortality improvement
#'
#' @param i Effective annual interest rate.
#'
#' @rdname mortality_improvement_projection
#' @export
axn_improved <- function(x0, n, i, qx_base_vec, AAx_vec, base_year, issue_year) {
  args <- .recycle_improved_annuity_inputs(x0 = x0, n = n, i = i)
  x0 <- args$x0
  n <- args$n
  i <- args$i

  sapply(seq_along(x0), function(j) {
    ii <- .check_improved_interest(i[j])
    nn <- .check_nonnegative_integer_scalar(n[j], "n")
    if (nn == 0L) return(0)

    survs <- sapply(1:nn, function(t) {
      tpx_improved(x0[j], t, qx_base_vec, AAx_vec, base_year, issue_year)
    })

    sum(discount(ii, 1:nn) * survs)
  })
}

#' Deferred temporary annuity-immediate under mortality improvement
#'
#' @param u Deferral period in years.
#'
#' @rdname mortality_improvement_projection
#' @export
naxn_improved <- function(x0, u, n, i, qx_base_vec, AAx_vec, base_year, issue_year) {
  args <- .recycle_improved_annuity_inputs(x0 = x0, u = u, n = n, i = i)
  x0 <- args$x0
  u <- args$u
  n <- args$n
  i <- args$i

  sapply(seq_along(x0), function(j) {
    ii <- .check_improved_interest(i[j])
    uu <- .check_nonnegative_integer_scalar(u[j], "u")
    nn <- .check_nonnegative_integer_scalar(n[j], "n")
    if (nn == 0L) return(0)

    times <- (uu + 1L):(uu + nn)

    survs <- sapply(times, function(t) {
      tpx_improved(x0[j], t, qx_base_vec, AAx_vec, base_year, issue_year)
    })

    sum(discount(ii, times) * survs)
  })
}

#' Whole life annuity-immediate under mortality improvement
#'
#' @rdname mortality_improvement_projection
#' @export
ax_improved <- function(x0, i, qx_base_vec, AAx_vec, base_year, issue_year) {
  args <- .recycle_improved_annuity_inputs(x0 = x0, i = i)
  x0 <- args$x0
  i <- args$i

  qx_base_vec <- .check_prob_vector(qx_base_vec, "qx_base_vec")
  AAx_vec <- .check_prob_vector(AAx_vec, "AAx_vec")

  L <- min(length(qx_base_vec), length(AAx_vec))
  if (L == 0L) return(rep(0, length(x0)))

  sapply(seq_along(x0), function(j) {
    ii <- .check_improved_interest(i[j])

    survs <- sapply(1:L, function(t) {
      tpx_improved(
        x0 = x0[j],
        n = t,
        qx_base_vec = qx_base_vec[1:L],
        AAx_vec = AAx_vec[1:L],
        base_year = base_year,
        issue_year = issue_year
      )
    })

    sum(discount(ii, 1:L) * survs)
  })
}
