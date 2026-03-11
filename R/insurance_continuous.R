#' Continuous insurance models (Chapter 7)
#'
#' Continuous contingent payment / insurance functions matching Chapter 7 notation.
#'
#' These functions handle:
#' \itemize{
#'   \item continuous whole life insurance: \eqn{\bar{A}_x},
#'   \item continuous term insurance: \eqn{\bar{A}_{x:\angl{n}}^{1}},
#'   \item continuous deferred insurance: \eqn{{}_{n\mid}\bar{A}_x},
#'   \item continuous endowment insurance: \eqn{\bar{A}_{x:\angl{n}}},
#'   \item second moments and variances.
#' }
#'
#' These functions are evaluated from the parametric survival-model framework
#' developed earlier in the package.
#'
#' @name insurance_continuous
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

#' @noRd
.check_cont_inputs <- function(x, i) {
  x <- as.numeric(x)
  i <- as.numeric(i)

  if (any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be >= 0.", call. = FALSE)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value > -1.", call. = FALSE)
  }

  list(x = x, i = i)
}

#' @noRd
.recycle_x_n_cont <- function(x, n) {
  x <- as.numeric(x)
  n <- as.numeric(n)

  if (length(x) == 0L || length(n) == 0L) {
    stop("x and n must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x))) stop("x must be finite.", call. = FALSE)
  if (any(!is.finite(n))) stop("n must be finite.", call. = FALSE)
  if (any(x < 0)) stop("x must be >= 0.", call. = FALSE)
  if (any(n < 0)) stop("n must be >= 0.", call. = FALSE)

  if (length(x) == 1L && length(n) > 1L) x <- rep(x, length(n))
  if (length(n) == 1L && length(x) > 1L) n <- rep(n, length(x))
  if (length(x) != length(n)) {
    stop("x and n must have compatible lengths.", call. = FALSE)
  }

  list(x = x, n = n)
}

#' @noRd
.integrate_cont <- function(f, lower, upper) {
  stats::integrate(
    f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-9,
    subdivisions = 1000L
  )$value
}

#' @noRd
.get_upper_cont <- function(x, n = NULL, model, ...) {
  if (tolower(model) == "uniform") {
    dots <- list(...)
    if (is.null(dots$omega)) {
      stop("For model = 'uniform', omega must be supplied.", call. = FALSE)
    }
    upper <- dots$omega - x
    if (!is.null(n)) upper <- pmin(upper, n)
    return(pmax(upper, 0))
  }

  if (is.null(n)) return(Inf)
  n
}

#' @noRd
.pure_endowment_cont <- function(x, n, i, model, ...) {
  discount(i, n) * tpx(n, x = x, model = model, ...)
}

# ------------------------------------------------------------
# First moments
# ------------------------------------------------------------

#' Continuous whole life insurance APV
#'
#' Computes
#' \eqn{\bar{A}_x = \int_0^\infty v^t \, {}_t p_x \mu_{x+t}\,dt}.
#'
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Abarx <- function(x, i, model, ...) {
  chk <- .check_cont_inputs(x, i)
  x <- chk$x

  delta <- interest_convert(i = i)$delta

  sapply(x, function(xx) {
    upper <- .get_upper_cont(xx, model = model, ...)

    if (upper <= 0) return(0)

    f <- function(t) {
      val <- exp(-delta * t) *
        tpx(t, x = xx, model = model, ...) *
        hazard0(xx + t, model = model, ...)
      val[!is.finite(val)] <- 0
      val
    }

    .integrate_cont(f, 0, upper)
  })
}

#' Continuous term insurance APV
#'
#' Computes
#' \eqn{\bar{A}_{x:\angl{n}}^{1} = \int_0^n v^t \, {}_t p_x \mu_{x+t}\,dt}.
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Abarxn1 <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_cont(x, n)
  x <- xn$x
  n <- xn$n

  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single value > -1.", call. = FALSE)
  }

  delta <- interest_convert(i = i)$delta

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]

    upper <- .get_upper_cont(xx, n = nn, model = model, ...)
    if (upper <= 0) return(0)

    f <- function(t) {
      val <- exp(-delta * t) *
        tpx(t, x = xx, model = model, ...) *
        hazard0(xx + t, model = model, ...)
      val[!is.finite(val)] <- 0
      val
    }

    .integrate_cont(f, 0, upper)
  })
}

#' Continuous deferred insurance APV
#'
#' Computes
#' \eqn{{}_{n\mid}\bar{A}_x = v^n\,{}_np_x\,\bar{A}_{x+n}}.
#'
#' @param x Age.
#' @param n Deferral period.
#' @param i Effective annual interest rate.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
nAbarx <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_cont(x, n)
  x <- xn$x
  n <- xn$n

  .pure_endowment_cont(x, n, i, model = model, ...) *
    Abarx(x + n, i, model = model, ...)
}

#' Continuous endowment insurance APV
#'
#' Computes
#' \eqn{\bar{A}_{x:\angl{n}} = \bar{A}_{x:\angl{n}}^{1} + v^n\,{}_np_x}.
#'
#' @param x Age.
#' @param n Term.
#' @param i Effective annual interest rate.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Abarxn <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_cont(x, n)
  x <- xn$x
  n <- xn$n

  Abarxn1(x, n, i, model = model, ...) +
    .pure_endowment_cont(x, n, i, model = model, ...)
}

# ------------------------------------------------------------
# Second moments
# ------------------------------------------------------------

#' Second moment of continuous whole life insurance PV
#'
#' Computes \eqn{{}^{2}\bar{A}_x} by evaluating \eqn{\bar{A}_x}
#' at doubled force.
#'
#' @inheritParams Abarx
#' @return Numeric vector of second moments.
#' @export
A2barx <- function(x, i, model, ...) {
  Abarx(x, i = double_force_i(i), model = model, ...)
}

#' Second moment of continuous term insurance PV
#'
#' Computes \eqn{{}^{2}\bar{A}_{x:\angl{n}}^{1}} by evaluating
#' \eqn{\bar{A}_{x:\angl{n}}^{1}} at doubled force.
#'
#' @inheritParams Abarxn1
#' @return Numeric vector of second moments.
#' @export
A2barxn1 <- function(x, n, i, model, ...) {
  Abarxn1(x, n, i = double_force_i(i), model = model, ...)
}

#' Second moment of continuous deferred insurance PV
#'
#' Computes \eqn{{}^{2}{}_{n\mid}\bar{A}_x} by evaluating
#' \eqn{{}_{n\mid}\bar{A}_x} at doubled force.
#'
#' @inheritParams nAbarx
#' @return Numeric vector of second moments.
#' @export
A2nAbarx <- function(x, n, i, model, ...) {
  nAbarx(x, n, i = double_force_i(i), model = model, ...)
}

#' Second moment of continuous endowment insurance PV
#'
#' Computes \eqn{{}^{2}\bar{A}_{x:\angl{n}}}.
#'
#' @inheritParams Abarxn
#' @return Numeric vector of second moments.
#' @export
A2barxn <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_cont(x, n)
  x <- xn$x
  n <- xn$n

  A2barxn1(x, n, i, model = model, ...) +
    .pure_endowment_cont(x, n, double_force_i(i), model = model, ...)
}

# ------------------------------------------------------------
# Variances
# ------------------------------------------------------------

#' Variance of continuous whole life insurance PV
#'
#' Computes \eqn{\mathrm{Var}(\bar{Z}_x) = {}^{2}\bar{A}_x - \bar{A}_x^2}.
#'
#' @inheritParams Abarx
#' @return Numeric vector of variances.
#' @export
var_Abarx <- function(x, i, model, ...) {
  A2barx(x, i, model = model, ...) -
    Abarx(x, i, model = model, ...)^2
}

#' Variance of continuous term insurance PV
#'
#' Computes
#' \eqn{\mathrm{Var}(\bar{Z}_{x:\angl{n}}^{1}) = {}^{2}\bar{A}_{x:\angl{n}}^{1} - (\bar{A}_{x:\angl{n}}^{1})^2}.
#'
#' @inheritParams Abarxn1
#' @return Numeric vector of variances.
#' @export
var_Abarxn1 <- function(x, n, i, model, ...) {
  A2barxn1(x, n, i, model = model, ...) -
    Abarxn1(x, n, i, model = model, ...)^2
}

#' Variance of continuous deferred insurance PV
#'
#' @inheritParams nAbarx
#' @return Numeric vector of variances.
#' @export
var_nAbarx <- function(x, n, i, model, ...) {
  A2nAbarx(x, n, i, model = model, ...) -
    nAbarx(x, n, i, model = model, ...)^2
}

#' Variance of continuous endowment insurance PV
#'
#' @inheritParams Abarxn
#' @return Numeric vector of variances.
#' @export
var_Abarxn <- function(x, n, i, model, ...) {
  A2barxn(x, n, i, model = model, ...) -
    Abarxn(x, n, i, model = model, ...)^2
}
