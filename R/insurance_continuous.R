#' Continuous insurance models
#'
#' Continuous contingent payment and insurance functions.
#'
#' @name insurance_continuous
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

.check_i_cont <- function(i) {
  i <- as.numeric(i)

  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }

  i
}

.recycle_cont_vectors <- function(...) {
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


.recycle_x_n_cont <- function(x, n) {
  args <- .recycle_cont_vectors(x, n)

  x <- args[[1L]]
  n <- args[[2L]]

  .check_x_cont(x)
  .check_n_cont(n)

  list(x = x, n = n)
}


.check_x_cont <- function(x) {
  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  invisible(TRUE)
}

.check_n_cont <- function(n) {
  if (any(n < 0)) {
    stop("n must be nonnegative.", call. = FALSE)
  }

  invisible(TRUE)
}

.integrate_cont <- function(f, lower, upper) {
  stats::integrate(
    f,
    lower = lower,
    upper = upper,
    rel.tol = 1e-9,
    subdivisions = 1000L
  )$value
}

.get_upper_cont <- function(x, n = NULL, model, ...) {
  if (tolower(model) == "uniform") {
    dots <- list(...)

    if (is.null(dots$omega)) {
      stop("For model = 'uniform', omega must be supplied.", call. = FALSE)
    }

    upper <- dots$omega - x

    if (!is.null(n)) {
      upper <- pmin(upper, n)
    }

    return(pmax(upper, 0))
  }

  if (is.null(n)) {
    return(Inf)
  }

  n
}

.pure_endowment_cont <- function(x, n, i, model, ...) {
  discount(i, n) * tpx(n, x = x, model = model, ...)
}

# ------------------------------------------------------------
# First moments
# ------------------------------------------------------------

#' Continuous whole life insurance APV
#'
#' Computes
#' \eqn{\bar{A}_x = \int_0^\infty v^t {}_t p_x \mu_{x+t}\,dt}.
#'
#' @param x Age.
#' @param i Effective annual interest rate.
#' @param model Parametric survival model name.
#' @param ... Additional model parameters passed to survival-model functions.
#'
#' @return Numeric vector of APVs.
#' @export
Abarx <- function(x, i, model, ...) {
  args <- .recycle_cont_vectors(x, .check_i_cont(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_x_cont(x)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    delta <- log(1 + i[j])

    upper <- .get_upper_cont(xx, model = model, ...)

    if (upper <= 0) {
      return(0)
    }

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
#' \eqn{\bar{A}_{x:\overline{n}|}^{1} =
#' \int_0^n v^t {}_t p_x \mu_{x+t}\,dt}.
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
  args <- .recycle_cont_vectors(x, n, .check_i_cont(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_cont(x)
  .check_n_cont(n)

  sapply(seq_along(x), function(j) {
    xx <- x[j]
    nn <- n[j]
    delta <- log(1 + i[j])

    upper <- .get_upper_cont(xx, n = nn, model = model, ...)

    if (upper <= 0) {
      return(0)
    }

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
#' \eqn{{}_{n\mid}\bar{A}_x = v^n {}_n p_x \bar{A}_{x+n}}.
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
  args <- .recycle_cont_vectors(x, n, .check_i_cont(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_x_cont(x)
  .check_n_cont(n)

  .pure_endowment_cont(x, n, i, model = model, ...) *
    Abarx(x + n, i, model = model, ...)
}

#' Continuous endowment insurance APV
#'
#' Computes
#' \eqn{\bar{A}_{x:\overline{n}|} =
#' \bar{A}_{x:\overline{n}|}^{1} + v^n {}_n p_x}.
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
  args <- .recycle_cont_vectors(x, n, .check_i_cont(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

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
#' Computes \eqn{{}^{2}\bar{A}_{x:\overline{n}|}^{1}} by evaluating
#' \eqn{\bar{A}_{x:\overline{n}|}^{1}} at doubled force.
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
#' Computes \eqn{{}^{2}\bar{A}_{x:\overline{n}|}}.
#'
#' @inheritParams Abarxn
#' @return Numeric vector of second moments.
#' @export
A2barxn <- function(x, n, i, model, ...) {
  args <- .recycle_cont_vectors(x, n, .check_i_cont(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  A2barxn1(x, n, i, model = model, ...) +
    .pure_endowment_cont(x, n, double_force_i(i), model = model, ...)
}

# ------------------------------------------------------------
# Variances
# ------------------------------------------------------------

#' Variance of continuous whole life insurance PV
#'
#' Computes \eqn{\mathrm{Var}(\bar{Z}_x) =
#' {}^{2}\bar{A}_x - \bar{A}_x^2}.
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
#' \eqn{\mathrm{Var}(\bar{Z}_{x:\overline{n}|}^{1}) =
#' {}^{2}\bar{A}_{x:\overline{n}|}^{1}
#' - (\bar{A}_{x:\overline{n}|}^{1})^2}.
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
