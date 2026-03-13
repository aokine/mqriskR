#' Chapter 9 premium, loss, and expense functions
#'
#' Functions in this file implement Chapter 9 funding-plan formulas,
#' including:
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
#' Naming follows Chapter 9 notation as closely as possible:
#' \itemize{
#'   \item \code{Px()} = whole life annual premium
#'   \item \code{Pxn1()} = term insurance annual premium
#'   \item \code{PnEx()} = pure endowment annual premium
#'   \item \code{Pxn()} = endowment insurance annual premium
#'   \item \code{tPx()} = limited-payment whole life annual premium
#'   \item \code{tPxn1()} = limited-payment term insurance annual premium
#'   \item \code{tPnEx()} = limited-payment pure endowment annual premium
#'   \item \code{tPxn()} = limited-payment endowment insurance annual premium
#'   \item \code{PnAx()} = deferred insurance annual premium
#'   \item \code{tPnAx()} = limited-payment deferred insurance annual premium
#'   \item \code{Pbarx()} = continuous-payment premium for discrete whole life insurance
#'   \item \code{Pbarxn1()} = continuous-payment premium for discrete term insurance
#'   \item \code{Pbarxn()} = continuous-payment premium for discrete endowment insurance
#'   \item \code{PbarAbarx()} = fully continuous premium for continuous whole life insurance
#'   \item \code{PbarAbarxn1()} = fully continuous premium for continuous term insurance
#'   \item \code{PbarAbarxn()} = fully continuous premium for continuous endowment insurance
#'   \item \code{Px_m()} = true fractional whole life annual premium
#'   \item \code{Pxn1_m()} = true fractional term insurance annual premium
#'   \item \code{Pxn_m()} = true fractional endowment insurance annual premium
#'   \item \code{PnAx_m()} = true fractional deferred insurance annual premium
#' }
#'
#' The discrete premium functions can be evaluated from either a life table
#' via \code{tbl = ...} or from a parametric model via \code{model = ...}.
#'
#' The continuous-premium-rate functions use the continuous annuity functions
#' already in the package, so they are written for the parametric survival
#' model framework.
#'
#' @name premium_ch9
#' @aliases Px Pxn1 PnEx Pxn tPx tPxn1 tPnEx tPxn PnAx tPnAx
#' @aliases Pbarx Pbarxn1 Pbarxn PbarAbarx PbarAbarxn1 PbarAbarxn
#' @aliases Px_m Pxn1_m Pxn_m PnAx_m
#' @aliases EL0x varL0x EL0xn1 varL0xn1 EL0xn varL0xn
#' @aliases EL0barAbarx varL0barAbarx
#' @aliases Gx
#' @param x Age.
#' @param n Term.
#' @param t Premium-paying period.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param P Premium amount or premium rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#' @param tol Numerical tolerance for functions that truncate infinite sums.
#' @param k_max Maximum summation horizon for functions that truncate infinite sums.
#' @param benefit Benefit amount.
#' @param first_premium_pct First-year premium expense proportion.
#' @param renewal_premium_pct Renewal premium expense proportion.
#' @param first_policy_exp First-year fixed expense.
#' @param renewal_policy_exp Renewal fixed expense each year after the first.
#' @param settlement_exp Settlement expense incurred at benefit payment.
#' @return Numeric vector.
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_i_ch9 <- function(i) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite number greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.check_m_ch9 <- function(m) {
  m <- as.numeric(m)
  if (length(m) != 1L || !is.finite(m) || m <= 0 || abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }
  as.integer(round(m))
}

#' @noRd
.check_t_scalar_ch9 <- function(t, name = "t") {
  t <- as.numeric(t)
  if (length(t) != 1L || !is.finite(t) || t < 0 || abs(t - round(t)) > 1e-10) {
    stop(sprintf("%s must be a single nonnegative integer.", name), call. = FALSE)
  }
  as.integer(round(t))
}

#' @noRd
.recycle_xn_ch9 <- function(x, n) {
  x <- as.numeric(x)
  n <- as.numeric(n)

  if (length(x) == 0L || length(n) == 0L) {
    stop("x and n must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("x must contain nonnegative finite values.", call. = FALSE)
  }
  if (any(!is.finite(n)) || any(n < 0) || any(abs(n - round(n)) > 1e-10)) {
    stop("n must contain nonnegative integer values.", call. = FALSE)
  }

  if (length(x) == 1L && length(n) > 1L) x <- rep(x, length(n))
  if (length(n) == 1L && length(x) > 1L) n <- rep(n, length(x))

  if (length(x) != length(n)) {
    stop("x and n must have the same length, or one must have length 1.", call. = FALSE)
  }

  list(x = x, n = as.integer(round(n)))
}

#' @noRd
.recycle_xnt_ch9 <- function(x, n, t) {
  xn <- .recycle_xn_ch9(x, n)
  x <- xn$x
  n <- xn$n
  t <- as.numeric(t)

  if (length(t) == 0L || any(!is.finite(t)) || any(t < 0) || any(abs(t - round(t)) > 1e-10)) {
    stop("t must contain nonnegative integer values.", call. = FALSE)
  }

  if (length(t) == 1L && length(x) > 1L) t <- rep(t, length(x))
  if (length(x) == 1L && length(t) > 1L) {
    x <- rep(x, length(t))
    n <- rep(n, length(t))
  }

  if (length(x) != length(t)) {
    stop("x, n, and t must have compatible lengths.", call. = FALSE)
  }

  list(x = x, n = n, t = as.integer(round(t)))
}

# -------------------------------------------------------------------------
# Annual net premiums
# -------------------------------------------------------------------------

#' @rdname premium_ch9
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to survival-model functions.
#' @return Numeric vector.
#' @export
Px <- function(x, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  Ax(x, i = i, tbl = tbl, model = model, ...) /
    adotx(x, i = i, model = model, ..., k_max = 5000, tol = 1e-12)
}

#' @rdname premium_ch9
#' @param n Term.
#' @return Numeric vector.
#' @export
Pxn1 <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  Axn1(x, n, i = i, tbl = tbl, model = model, ...) /
    adotxn(x, n, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
PnEx <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  nEx(x, n, i = i, tbl = tbl, model = model, ...) /
    adotxn(x, n, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
Pxn <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  Axn(x, n, i = i, tbl = tbl, model = model, ...) /
    adotxn(x, n, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @param t Premium-paying period.
#' @return Numeric vector.
#' @export
tPx <- function(x, t, i, tbl = NULL, model = NULL, ...) {
  tt <- .check_t_scalar_ch9(t, "t")
  .check_i_ch9(i)

  Ax(x, i = i, tbl = tbl, model = model, ...) /
    adotxn(x, tt, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
tPxn1 <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  xnt <- .recycle_xnt_ch9(x, n, t)
  if (any(xnt$t > xnt$n)) {
    stop("For t-pay n-year contracts, t must satisfy t <= n.", call. = FALSE)
  }

  Axn1(xnt$x, xnt$n, i = i, tbl = tbl, model = model, ...) /
    adotxn(xnt$x, xnt$t, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
tPnEx <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  xnt <- .recycle_xnt_ch9(x, n, t)
  if (any(xnt$t > xnt$n)) {
    stop("For t-pay n-year contracts, t must satisfy t <= n.", call. = FALSE)
  }

  nEx(xnt$x, xnt$n, i = i, tbl = tbl, model = model, ...) /
    adotxn(xnt$x, xnt$t, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
tPxn <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  xnt <- .recycle_xnt_ch9(x, n, t)
  if (any(xnt$t > xnt$n)) {
    stop("For t-pay n-year contracts, t must satisfy t <= n.", call. = FALSE)
  }

  Axn(xnt$x, xnt$n, i = i, tbl = tbl, model = model, ...) /
    adotxn(xnt$x, xnt$t, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
PnAx <- function(x, n, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  nAx(x, n, i = i, tbl = tbl, model = model, ...) /
    adotx(x, i = i, model = model, ..., k_max = 5000, tol = 1e-12)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
tPnAx <- function(x, n, t, i, tbl = NULL, model = NULL, ...) {
  xnt <- .recycle_xnt_ch9(x, n, t)
  if (any(xnt$t > xnt$n)) {
    stop("For deferred insurance with limited premiums, t must satisfy t <= n.", call. = FALSE)
  }

  nAx(xnt$x, xnt$n, i = i, tbl = tbl, model = model, ...) /
    adotxn(xnt$x, xnt$t, i = i, model = model, ...)
}

# -------------------------------------------------------------------------
# Continuous payment premiums for discrete benefits
# -------------------------------------------------------------------------

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
Pbarx <- function(x, i, model, ..., tol = 1e-10) {
  .check_i_ch9(i)
  Ax(x, i = i, model = model, ...) /
    abarx(x, i = i, model = model, ..., tol = tol)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
Pbarxn1 <- function(x, n, i, model, ...) {
  .check_i_ch9(i)
  Axn1(x, n, i = i, model = model, ...) /
    abarxn(x, n, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
Pbarxn <- function(x, n, i, model, ...) {
  .check_i_ch9(i)
  Axn(x, n, i = i, model = model, ...) /
    abarxn(x, n, i = i, model = model, ...)
}

# -------------------------------------------------------------------------
# Fully continuous premium rates
# -------------------------------------------------------------------------

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
PbarAbarx <- function(x, i, model, ..., tol = 1e-10) {
  .check_i_ch9(i)
  Abarx(x, i = i, model = model, ...) /
    abarx(x, i = i, model = model, ..., tol = tol)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
PbarAbarxn1 <- function(x, n, i, model, ...) {
  .check_i_ch9(i)
  Abarxn1(x, n, i = i, model = model, ...) /
    abarxn(x, n, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
PbarAbarxn <- function(x, n, i, model, ...) {
  .check_i_ch9(i)
  Abarxn(x, n, i = i, model = model, ...) /
    abarxn(x, n, i = i, model = model, ...)
}

# -------------------------------------------------------------------------
# True fractional premiums
# -------------------------------------------------------------------------

#' @rdname premium_ch9
#' @param m Number of payments per year.
#' @return Numeric vector.
#' @export
Px_m <- function(x, m, i, tbl = NULL, model = NULL, ...) {
  m <- .check_m_ch9(m)
  Ax(x, i = i, tbl = tbl, model = model, ...) /
    adotx_m(x, m = m, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
Pxn1_m <- function(x, n, m, i, tbl = NULL, model = NULL, ...) {
  m <- .check_m_ch9(m)
  Axn1(x, n, i = i, tbl = tbl, model = model, ...) /
    adotxn_m(x, n, m = m, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
Pxn_m <- function(x, n, m, i, tbl = NULL, model = NULL, ...) {
  m <- .check_m_ch9(m)
  Axn(x, n, i = i, tbl = tbl, model = model, ...) /
    adotxn_m(x, n, m = m, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
PnAx_m <- function(x, n, m, i, tbl = NULL, model = NULL, ...) {
  m <- .check_m_ch9(m)
  nAx(x, n, i = i, tbl = tbl, model = model, ...) /
    adotx_m(x, m = m, i = i, model = model, ...)
}

# -------------------------------------------------------------------------
# Present-value-of-loss means and variances
# -------------------------------------------------------------------------

#' @rdname premium_ch9
#' @param P Premium amount or premium rate.
#' @return Numeric vector.
#' @export
EL0x <- function(x, P, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  Ax(x, i = i, tbl = tbl, model = model, ...) -
    P * adotx(x, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @return Numeric vector.
#' @export
varL0x <- function(x, P, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  d <- i / (1 + i)

  (1 + P / d)^2 *
    var_Ax(x, i = i, tbl = tbl, model = model, ..., tol = tol, k_max = k_max)
}

#' @rdname premium_ch9
#' @export
EL0xn1 <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  Axn1(x, n, i = i, tbl = tbl, model = model, ...) -
    P * adotxn(x, n, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @export
varL0xn1 <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  d <- i / (1 + i)

  (1 + P / d)^2 *
    var_Axn1(x, n, i = i, tbl = tbl, model = model, ...)
}

#' @rdname premium_ch9
#' @export
EL0xn <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  Axn(x, n, i = i, tbl = tbl, model = model, ...) -
    P * adotxn(x, n, i = i, model = model, ...)
}

#' @rdname premium_ch9
#' @export
varL0xn <- function(x, n, P, i, tbl = NULL, model = NULL, ...) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  d <- i / (1 + i)

  (1 + P / d)^2 *
    var_Axn(x, n, i = i, tbl = tbl, model = model, ...)
}

#' @rdname premium_ch9
#' @export
EL0barAbarx <- function(x, P, i, model, ..., tol = 1e-10) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  Abarx(x, i = i, model = model, ...) -
    P * abarx(x, i = i, model = model, ..., tol = tol)
}

#' @rdname premium_ch9
#' @export
varL0barAbarx <- function(x, P, i, model, ...) {
  .check_i_ch9(i)
  P <- as.numeric(P)
  delta <- log(1 + i)

  (1 + P / delta)^2 *
    var_Abarx(x, i = i, model = model, ...)
}

# -------------------------------------------------------------------------
# Gross premium with expenses: whole life, fully discrete
# -------------------------------------------------------------------------

#' @rdname premium_ch9
#' @param benefit Benefit amount.
#' @param first_premium_pct First-year premium expense proportion.
#' @param renewal_premium_pct Renewal premium expense proportion.
#' @param first_policy_exp First-year fixed expense.
#' @param renewal_policy_exp Renewal fixed expense each year after the first.
#' @param settlement_exp Settlement expense incurred at benefit payment.
#' @return Numeric vector.
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
  .check_i_ch9(i)

  Ax_val <- Ax(x, i = i, tbl = tbl, model = model, ...)
  ax_val <- ax(x, i = i, model = model, ...)
  adotx_val <- adotx(x, i = i, model = model, ...)

  num <- (benefit + settlement_exp) * Ax_val +
    first_policy_exp +
    renewal_policy_exp * ax_val

  den <- adotx_val -
    first_premium_pct -
    renewal_premium_pct * ax_val

  if (any(den <= 0)) {
    stop("Gross premium denominator must be positive.", call. = FALSE)
  }

  num / den
}
