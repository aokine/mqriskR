#' Premium, loss, and expense functions
#'
#' Functions for net and gross premiums, present-value-of-loss moments,
#' continuous-payment premium rates, and premiums payable more frequently
#' than annually.
#'
#' The functions include:
#' \itemize{
#'   \item net annual premiums under the equivalence principle,
#'   \item limited-payment premiums,
#'   \item continuous-payment premium rates,
#'   \item fully continuous premium rates,
#'   \item true fractional premiums,
#'   \item present-value-of-loss means and variances,
#'   \item a basic gross premium formula for whole life insurance.
#' }
#'
#' The discrete premium functions may be evaluated from either a life table
#' supplied through \code{tbl} or a parametric survival model supplied through
#' \code{model}. Continuous-payment premium functions are evaluated through
#' the corresponding continuous insurance and annuity functions.
#'
#' Scalar inputs retain their existing behavior. Where mathematically
#' meaningful, numeric inputs may also be vectors. Inputs are evaluated
#' elementwise when they have a common length, and scalar inputs are recycled.
#'
#' @param x Age. May be scalar or vector.
#' @param n Term. May be scalar or vector of nonnegative integers.
#' @param t Premium-paying period. May be scalar or vector of nonnegative integers.
#' @param m Number of payments per year. Must be a positive integer scalar.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param P Premium amount or premium rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#' @param tol Numerical tolerance for functions that truncate infinite sums.
#' @param k_max Maximum summation horizon for functions that truncate infinite sums.
#' @param benefit Benefit amount. May be scalar or vector.
#' @param first_premium_pct First-year premium expense proportion. May be scalar or vector.
#' @param renewal_premium_pct Renewal premium expense proportion. May be scalar or vector.
#' @param first_policy_exp First-year fixed expense. May be scalar or vector.
#' @param renewal_policy_exp Renewal fixed expense after the first year. May be scalar or vector.
#' @param settlement_exp Settlement expense incurred at benefit payment. May be scalar or vector.
#'
#' @return Numeric vector.
#' @name premium_functions
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_premium_interest <- function(i) {
  i <- as.numeric(i)
  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.check_premium_frequency <- function(m) {
  m <- as.numeric(m)
  if (length(m) != 1L || !is.finite(m) || m <= 0 ||
      abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }
  as.integer(round(m))
}

#' @noRd
.check_premium_age <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L || any(!is.finite(x)) || any(x < 0)) {
    stop("x must contain nonnegative finite values.", call. = FALSE)
  }
  x
}

#' @noRd
.check_premium_integer <- function(z, name) {
  z <- as.numeric(z)
  if (length(z) == 0L || any(!is.finite(z)) || any(z < 0) ||
      any(abs(z - round(z)) > 1e-10)) {
    stop(sprintf("%s must contain nonnegative integer values.", name),
         call. = FALSE)
  }
  as.integer(round(z))
}

#' @noRd
.check_premium_numeric <- function(z, name) {
  z <- as.numeric(z)
  if (length(z) == 0L || any(!is.finite(z))) {
    stop(sprintf("%s must contain finite numeric values.", name),
         call. = FALSE)
  }
  z
}

#' @noRd
.recycle_premium_vectors <- function(...) {
  values <- list(...)
  lens <- vapply(values, length, integer(1))

  if (length(values) == 0L || any(lens == 0L)) {
    stop("Arguments must have positive length.", call. = FALSE)
  }

  target <- max(lens)
  if (!all(lens %in% c(1L, target))) {
    stop("Arguments must have compatible lengths.", call. = FALSE)
  }

  lapply(values, rep_len, length.out = target)
}

#' @noRd
.validate_premium_source <- function(tbl, model) {
  if (is.null(tbl) && is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.map_premium <- function(args, FUN) {
  n <- length(args[[1L]])
  vapply(seq_len(n), function(j) {
    current <- lapply(args, `[[`, j)
    do.call(FUN, current)
  }, numeric(1))
}

# -------------------------------------------------------------------------
# Annual net premiums
# -------------------------------------------------------------------------

#' Whole life annual net premium
#'
#' @rdname premium_functions
#' @export
Px <- function(x, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, i = i)

  .map_premium(args, function(x, i) {
    Ax(x, i = i, tbl = tbl, model = model, ...) /
      adotx(x, i = i, tbl = tbl, model = model, ...,
            k_max = 5000, tol = 1e-12)
  })
}

#' Term insurance annual net premium
#'
#' @rdname premium_functions
#' @export
Pxn1 <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Axn1(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, n, i = i, tbl = tbl, model = model, ...)
  })
}

#' Pure endowment annual net premium
#'
#' @rdname premium_functions
#' @export
PnEx <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    nEx(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, n, i = i, tbl = tbl, model = model, ...)
  })
}

#' Endowment insurance annual net premium
#'
#' @rdname premium_functions
#' @export
Pxn <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Axn(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, n, i = i, tbl = tbl, model = model, ...)
  })
}

#' Limited-payment whole life annual net premium
#'
#' @rdname premium_functions
#' @export
tPx <- function(x, t, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  t <- .check_premium_integer(t, "t")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, t = t, i = i)

  .map_premium(args, function(x, t, i) {
    Ax(x, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, t, i = i, tbl = tbl, model = model, ...)
  })
}

#' Limited-payment term insurance annual net premium
#'
#' @rdname premium_functions
#' @export
tPxn1 <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  t <- .check_premium_integer(t, "t")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, t = t, i = i)

  if (any(args$t > args$n)) {
    stop("For t-pay n-year contracts, t must satisfy t <= n.", call. = FALSE)
  }

  .map_premium(args, function(x, n, t, i) {
    Axn1(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, t, i = i, tbl = tbl, model = model, ...)
  })
}

#' Limited-payment pure endowment annual net premium
#'
#' @rdname premium_functions
#' @export
tPnEx <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  t <- .check_premium_integer(t, "t")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, t = t, i = i)

  if (any(args$t > args$n)) {
    stop("For t-pay n-year contracts, t must satisfy t <= n.", call. = FALSE)
  }

  .map_premium(args, function(x, n, t, i) {
    nEx(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, t, i = i, tbl = tbl, model = model, ...)
  })
}

#' Limited-payment endowment insurance annual net premium
#'
#' @rdname premium_functions
#' @export
tPxn <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  t <- .check_premium_integer(t, "t")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, t = t, i = i)

  if (any(args$t > args$n)) {
    stop("For t-pay n-year contracts, t must satisfy t <= n.", call. = FALSE)
  }

  .map_premium(args, function(x, n, t, i) {
    Axn(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, t, i = i, tbl = tbl, model = model, ...)
  })
}

#' Deferred insurance annual net premium
#'
#' @rdname premium_functions
#' @export
PnAx <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    nAx(x, n, i = i, tbl = tbl, model = model, ...) /
      adotx(x, i = i, tbl = tbl, model = model, ...,
            k_max = 5000, tol = 1e-12)
  })
}

#' Limited-payment deferred insurance annual net premium
#'
#' @rdname premium_functions
#' @export
tPnAx <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  t <- .check_premium_integer(t, "t")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, t = t, i = i)

  if (any(args$t > args$n)) {
    stop("For deferred insurance with limited premiums, t must satisfy t <= n.",
         call. = FALSE)
  }

  .map_premium(args, function(x, n, t, i) {
    nAx(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn(x, t, i = i, tbl = tbl, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# Continuous-payment premiums for discrete benefits
# -------------------------------------------------------------------------

#' Continuous-payment premium rate for discrete whole life insurance
#'
#' @rdname premium_functions
#' @export
Pbarx <- function(x, i, model, ..., tol = 1e-10) {
  x <- .check_premium_age(x)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, i = i)

  .map_premium(args, function(x, i) {
    Ax(x, i = i, model = model, ...) /
      abarx(x, i = i, model = model, ..., tol = tol)
  })
}

#' Continuous-payment premium rate for discrete term insurance
#'
#' @rdname premium_functions
#' @export
Pbarxn1 <- function(x, n, i, model, ...) {
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Axn1(x, n, i = i, model = model, ...) /
      abarxn(x, n, i = i, model = model, ...)
  })
}

#' Continuous-payment premium rate for discrete endowment insurance
#'
#' @rdname premium_functions
#' @export
Pbarxn <- function(x, n, i, model, ...) {
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Axn(x, n, i = i, model = model, ...) /
      abarxn(x, n, i = i, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# Fully continuous premium rates
# -------------------------------------------------------------------------

#' Fully continuous premium rate for whole life insurance
#'
#' @rdname premium_functions
#' @export
PbarAbarx <- function(x, i, model, ..., tol = 1e-10) {
  x <- .check_premium_age(x)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, i = i)

  .map_premium(args, function(x, i) {
    Abarx(x, i = i, model = model, ...) /
      abarx(x, i = i, model = model, ..., tol = tol)
  })
}

#' Fully continuous premium rate for term insurance
#'
#' @rdname premium_functions
#' @export
PbarAbarxn1 <- function(x, n, i, model, ...) {
  x <- .check_premium_age(x)
  n <- .check_premium_numeric(n, "n")
  if (any(n < 0)) stop("n must be nonnegative.", call. = FALSE)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Abarxn1(x, n, i = i, model = model, ...) /
      abarxn(x, n, i = i, model = model, ...)
  })
}

#' Fully continuous premium rate for endowment insurance
#'
#' @rdname premium_functions
#' @export
PbarAbarxn <- function(x, n, i, model, ...) {
  x <- .check_premium_age(x)
  n <- .check_premium_numeric(n, "n")
  if (any(n < 0)) stop("n must be nonnegative.", call. = FALSE)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Abarxn(x, n, i = i, model = model, ...) /
      abarxn(x, n, i = i, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# True fractional premiums
# -------------------------------------------------------------------------

#' True fractional whole life premium
#'
#' @rdname premium_functions
#' @export
Px_m <- function(x, m, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  m <- .check_premium_frequency(m)
  x <- .check_premium_age(x)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, i = i)

  .map_premium(args, function(x, i) {
    Ax(x, i = i, tbl = tbl, model = model, ...) /
      adotx_m(x, m = m, i = i, tbl = tbl, model = model, ...)
  })
}

#' True fractional term insurance premium
#'
#' @rdname premium_functions
#' @export
Pxn1_m <- function(x, n, m, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  m <- .check_premium_frequency(m)
  x <- .check_premium_age(x)
  n <- .check_premium_numeric(n, "n")
  if (any(n < 0)) stop("n must be nonnegative.", call. = FALSE)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Axn1(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn_m(x, n, m = m, i = i, tbl = tbl, model = model, ...)
  })
}

#' True fractional endowment insurance premium
#'
#' @rdname premium_functions
#' @export
Pxn_m <- function(x, n, m, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  m <- .check_premium_frequency(m)
  x <- .check_premium_age(x)
  n <- .check_premium_numeric(n, "n")
  if (any(n < 0)) stop("n must be nonnegative.", call. = FALSE)
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    Axn(x, n, i = i, tbl = tbl, model = model, ...) /
      adotxn_m(x, n, m = m, i = i, tbl = tbl, model = model, ...)
  })
}

#' True fractional deferred insurance premium
#'
#' @rdname premium_functions
#' @export
PnAx_m <- function(x, n, m, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  m <- .check_premium_frequency(m)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, i = i)

  .map_premium(args, function(x, n, i) {
    nAx(x, n, i = i, tbl = tbl, model = model, ...) /
      adotx_m(x, m = m, i = i, tbl = tbl, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# Present-value-of-loss means and variances
# -------------------------------------------------------------------------

#' Expected present value of loss for whole life insurance
#'
#' @rdname premium_functions
#' @export
EL0x <- function(x, P, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, P = P, i = i)

  .map_premium(args, function(x, P, i) {
    Ax(x, i = i, tbl = tbl, model = model, ...) -
      P * adotx(x, i = i, tbl = tbl, model = model, ...)
  })
}

#' Variance of present value of loss for whole life insurance
#'
#' @rdname premium_functions
#' @export
varL0x <- function(x, P, i, tbl = NULL, model = NULL, ...,
                   tol = 1e-12, k_max = 5000) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, P = P, i = i)

  .map_premium(args, function(x, P, i) {
    d <- i / (1 + i)
    (1 + P / d)^2 *
      var_Ax(x, i = i, tbl = tbl, model = model, ...,
             tol = tol, k_max = k_max)
  })
}

#' Expected present value of loss for term insurance
#'
#' @rdname premium_functions
#' @export
EL0xn1 <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, P = P, i = i)

  .map_premium(args, function(x, n, P, i) {
    Axn1(x, n, i = i, tbl = tbl, model = model, ...) -
      P * adotxn(x, n, i = i, tbl = tbl, model = model, ...)
  })
}

#' Variance of present value of loss for term insurance
#'
#' @rdname premium_functions
#' @export
varL0xn1 <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, P = P, i = i)

  .map_premium(args, function(x, n, P, i) {
    d <- i / (1 + i)
    (1 + P / d)^2 *
      var_Axn1(x, n, i = i, tbl = tbl, model = model, ...)
  })
}

#' Expected present value of loss for endowment insurance
#'
#' @rdname premium_functions
#' @export
EL0xn <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, P = P, i = i)

  .map_premium(args, function(x, n, P, i) {
    Axn(x, n, i = i, tbl = tbl, model = model, ...) -
      P * adotxn(x, n, i = i, tbl = tbl, model = model, ...)
  })
}

#' Variance of present value of loss for endowment insurance
#'
#' @rdname premium_functions
#' @export
varL0xn <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .validate_premium_source(tbl, model)
  x <- .check_premium_age(x)
  n <- .check_premium_integer(n, "n")
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, n = n, P = P, i = i)

  .map_premium(args, function(x, n, P, i) {
    d <- i / (1 + i)
    (1 + P / d)^2 *
      var_Axn(x, n, i = i, tbl = tbl, model = model, ...)
  })
}

#' Expected present value of loss for fully continuous whole life insurance
#'
#' @rdname premium_functions
#' @export
EL0barAbarx <- function(x, P, i, model, ..., tol = 1e-10) {
  x <- .check_premium_age(x)
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, P = P, i = i)

  .map_premium(args, function(x, P, i) {
    Abarx(x, i = i, model = model, ...) -
      P * abarx(x, i = i, model = model, ..., tol = tol)
  })
}

#' Variance of present value of loss for fully continuous whole life insurance
#'
#' @rdname premium_functions
#' @export
varL0barAbarx <- function(x, P, i, model, ...) {
  x <- .check_premium_age(x)
  P <- .check_premium_numeric(P, "P")
  i <- .check_premium_interest(i)
  args <- .recycle_premium_vectors(x = x, P = P, i = i)

  .map_premium(args, function(x, P, i) {
    delta <- log(1 + i)
    (1 + P / delta)^2 *
      var_Abarx(x, i = i, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# Gross premium with expenses: whole life, fully discrete
# -------------------------------------------------------------------------

#' Gross premium for whole life insurance with expenses
#'
#' @rdname premium_functions
#' @export
Gx <- function(x,
               i,
               benefit = 1,
               first_premium_pct = 0,
               renewal_premium_pct = 0,
               first_policy_exp = 0,
               renewal_policy_exp = 0,
               settlement_exp = 0,
               tbl = NULL,
               model = NULL,
               ...) {
  .validate_premium_source(tbl, model)

  x <- .check_premium_age(x)
  i <- .check_premium_interest(i)
  benefit <- .check_premium_numeric(benefit, "benefit")
  first_premium_pct <- .check_premium_numeric(
    first_premium_pct, "first_premium_pct"
  )
  renewal_premium_pct <- .check_premium_numeric(
    renewal_premium_pct, "renewal_premium_pct"
  )
  first_policy_exp <- .check_premium_numeric(
    first_policy_exp, "first_policy_exp"
  )
  renewal_policy_exp <- .check_premium_numeric(
    renewal_policy_exp, "renewal_policy_exp"
  )
  settlement_exp <- .check_premium_numeric(
    settlement_exp, "settlement_exp"
  )

  args <- .recycle_premium_vectors(
    x = x,
    i = i,
    benefit = benefit,
    first_premium_pct = first_premium_pct,
    renewal_premium_pct = renewal_premium_pct,
    first_policy_exp = first_policy_exp,
    renewal_policy_exp = renewal_policy_exp,
    settlement_exp = settlement_exp
  )

  .map_premium(
    args,
    function(x, i, benefit, first_premium_pct, renewal_premium_pct,
             first_policy_exp, renewal_policy_exp, settlement_exp) {
      Ax_val <- Ax(x, i = i, tbl = tbl, model = model, ...)
      ax_val <- ax(x, i = i, tbl = tbl, model = model, ...)
      adotx_val <- adotx(x, i = i, tbl = tbl, model = model, ...)

      numerator <- (benefit + settlement_exp) * Ax_val +
        first_policy_exp + renewal_policy_exp * ax_val

      denominator <- adotx_val -
        first_premium_pct - renewal_premium_pct * ax_val

      if (!is.finite(denominator) || denominator <= 0) {
        stop("Gross premium denominator must be positive.", call. = FALSE)
      }

      numerator / denominator
    }
  )
}

