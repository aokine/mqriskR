#' Annuity approximations (Chapter 8)
#'
#' Chapter 8 approximation formulas for m-thly and continuous life annuities.
#'
#' This file implements:
#' \itemize{
#'   \item UDD approximations for m-thly annuities,
#'   \item UDD approximations for continuous annuities,
#'   \item Woolhouse 2-term approximations,
#'   \item Woolhouse 3-term approximations.
#' }
#'
#' @name annuity_approximations
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_i_scalar_ann <- function(i) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite number greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.check_m_scalar_ann <- function(m) {
  m <- as.numeric(m)
  if (length(m) != 1L || !is.finite(m) || m <= 0 || abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }
  as.integer(round(m))
}

#' @noRd
.recycle_x_n_ann <- function(x, n) {
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
.alpha_m <- function(i, m) {
  i <- .check_i_scalar_ann(i)
  m <- .check_m_scalar_ann(m)

  ic <- interest_convert(i = i, m = m)
  d <- ic$d
  im <- ic$im
  dm <- m * ((im / m) / (1 + im / m))

  (i * d) / (im * dm)
}

#' @noRd
.beta_m <- function(i, m) {
  i <- .check_i_scalar_ann(i)
  m <- .check_m_scalar_ann(m)

  ic <- interest_convert(i = i, m = m)
  im <- ic$im
  dm <- m * ((im / m) / (1 + im / m))

  (i - im) / (im * dm)
}

#' @noRd
.gamma_m <- function(i, m) {
  i <- .check_i_scalar_ann(i)
  m <- .check_m_scalar_ann(m)

  ic <- interest_convert(i = i, m = m)
  d <- ic$d
  im <- ic$im
  dm <- m * ((im / m) / (1 + im / m))

  (dm - d) / (im * dm)
}
#' @noRd
.delta_from_i <- function(i) {
  interest_convert(i = .check_i_scalar_ann(i))$delta
}

#' @noRd
.mu_x <- function(x, model, ...) {
  hazard0(x, model = model, ...)
}

# -------------------------------------------------------------------------
# Shared documentation: UDD approximations
# -------------------------------------------------------------------------

#' UDD annuity approximations
#'
#' UDD-based approximations for Chapter 8 annuity functions.
#'
#' These functions implement the standard Uniform Distribution of Deaths
#' approximations linking annual, m-thly, and continuous annuity values.
#'
#' The exported functions documented on this page are:
#' \itemize{
#'   \item \code{ax_m_udd()}
#'   \item \code{axn_m_udd()}
#'   \item \code{nax_m_udd()}
#'   \item \code{adotx_m_udd()}
#'   \item \code{adotxn_m_udd()}
#'   \item \code{nadotx_m_udd()}
#'   \item \code{sxn_m_udd()}
#'   \item \code{sdotxn_m_udd()}
#'   \item \code{abarx_udd()}
#'   \item \code{abarxn_udd()}
#'   \item \code{nabarx_udd()}
#' }
#'
#' @name annuity_approximations_udd
#' @aliases ax_m_udd axn_m_udd nax_m_udd
#' @aliases adotx_m_udd adotxn_m_udd nadotx_m_udd
#' @aliases sxn_m_udd sdotxn_m_udd
#' @aliases abarx_udd abarxn_udd nabarx_udd
NULL

#' UDD approximation to m-thly whole life annuity-due
#'
#' Computes
#' \deqn{\ddot{a}_x^{(m)} \approx \alpha(m)\ddot{a}_x - \beta(m).}
#'
#' @param x Age.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @rdname annuity_approximations_udd
#' @export
adotx_m_udd <- function(x, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  x <- as.numeric(x)

  .alpha_m(i, m) * adotx(x, i, model = model, ...) - .beta_m(i, m)
}

#' UDD approximation to m-thly temporary annuity-due
#'
#' Computes
#' \deqn{\ddot{a}_{x:\angl{n}}^{(m)} \approx
#' \alpha(m)\ddot{a}_{x:\angl{n}} - \beta(m)(1-{}_nE_x).}
#'
#' @param n Term.
#' @rdname annuity_approximations_udd
#' @export
adotxn_m_udd <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  x <- xn$x
  n <- xn$n

  .alpha_m(i, m) * adotxn(x, n, i, model = model, ...) -
    .beta_m(i, m) * (1 - nEx(x, n, i, model = model, ...))
}

#' UDD approximation to m-thly deferred whole life annuity-due
#'
#' Computes
#' \deqn{{}_{n\mid}\ddot{a}_x^{(m)} \approx
#' \alpha(m)\,{}_{n\mid}\ddot{a}_x - \beta(m)\,{}_nE_x.}
#'
#' @rdname annuity_approximations_udd
#' @export
nadotx_m_udd <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  x <- xn$x
  n <- xn$n

  .alpha_m(i, m) * nadotx(x, n, i, model = model, ...) -
    .beta_m(i, m) * nEx(x, n, i, model = model, ...)
}

#' UDD approximation to m-thly whole life annuity-immediate
#'
#' Computes
#' \deqn{a_x^{(m)} \approx \alpha(m)a_x + \gamma(m).}
#'
#' @rdname annuity_approximations_udd
#' @export
ax_m_udd <- function(x, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  x <- as.numeric(x)

  .alpha_m(i, m) * ax(x, i, model = model, ...) + .gamma_m(i, m)
}

#' UDD approximation to m-thly temporary annuity-immediate
#'
#' Computes
#' \deqn{a_{x:\angl{n}}^{(m)} \approx
#' \alpha(m)a_{x:\angl{n}} + \gamma(m)(1-{}_nE_x).}
#'
#' @rdname annuity_approximations_udd
#' @export
axn_m_udd <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  x <- xn$x
  n <- xn$n

  .alpha_m(i, m) * axn(x, n, i, model = model, ...) +
    .gamma_m(i, m) * (1 - nEx(x, n, i, model = model, ...))
}

#' UDD approximation to m-thly deferred whole life annuity-immediate
#'
#' Computes
#' \deqn{{}_{n\mid}a_x^{(m)} \approx
#' \alpha(m)\,{}_{n\mid}a_x + \gamma(m)\,{}_nE_x.}
#'
#' @rdname annuity_approximations_udd
#' @export
nax_m_udd <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  x <- xn$x
  n <- xn$n

  .alpha_m(i, m) * nax(x, n, i, model = model, ...) +
    .gamma_m(i, m) * nEx(x, n, i, model = model, ...)
}

#' UDD approximation to m-thly temporary annuity accumulated value
#'
#' Computes
#' \deqn{\ddot{s}_{x:\angl{n}}^{(m)} \approx
#' \alpha(m)\ddot{s}_{x:\angl{n}} - \beta(m)\left(\frac{1}{{}_nE_x}-1\right).}
#'
#' @rdname annuity_approximations_udd
#' @export
sdotxn_m_udd <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  x <- xn$x
  n <- xn$n

  ex <- nEx(x, n, i, model = model, ...)
  .alpha_m(i, m) * sdotxn(x, n, i, model = model, ...) -
    .beta_m(i, m) * (1 / ex - 1)
}

#' UDD approximation to m-thly temporary annuity accumulated value
#'
#' Computes
#' \deqn{s_{x:\angl{n}}^{(m)} \approx
#' \alpha(m)s_{x:\angl{n}} + \gamma(m)\left(\frac{1}{{}_nE_x}-1\right).}
#'
#' @rdname annuity_approximations_udd
#' @export
sxn_m_udd <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  x <- xn$x
  n <- xn$n

  ex <- nEx(x, n, i, model = model, ...)
  .alpha_m(i, m) * sxn(x, n, i, model = model, ...) +
    .gamma_m(i, m) * (1 / ex - 1)
}

#' UDD approximation to continuous whole life annuity
#'
#' Computes
#' \deqn{\bar{a}_x \approx \frac{id}{\delta^2}\ddot{a}_x - \frac{i-\delta}{\delta^2}.}
#'
#' @rdname annuity_approximations_udd
#' @export
abarx_udd <- function(x, i, model, ...) {
  i <- .check_i_scalar_ann(i)
  delta <- .delta_from_i(i)
  d <- interest_convert(i = i)$d

  ((i * d) / delta^2) * adotx(x, i, model = model, ...) -
    (i - delta) / delta^2
}

#' UDD approximation to continuous temporary annuity
#'
#' Uses the identity
#' \deqn{\bar{a}_{x:\angl{n}} \approx \frac{1-\bar{A}_{x:\angl{n}}}{\delta}}
#' together with the package's existing Chapter 7 UDD insurance approximation
#' for \eqn{\bar{A}_{x:\angl{n}}}.
#'
#' Note that this function relies on the already-existing \code{Abarxn_udd()}
#' implementation in the package, so extra survival-model arguments are not used.
#'
#' @rdname annuity_approximations_udd
#' @export
abarxn_udd <- function(x, n, i, model, ...) {
  i <- .check_i_scalar_ann(i)
  xn <- .recycle_x_n_ann(x, n)
  delta <- .delta_from_i(i)

  # existing package function does not accept model/...
  (1 - Abarxn_udd(xn$x, xn$n, i)) / delta
}

#' UDD approximation to continuous deferred whole life annuity
#'
#' Computes
#' \deqn{{}_{n\mid}\bar{a}_x \approx {}_nE_x \, \bar{a}_{x+n}.}
#'
#' @rdname annuity_approximations_udd
#' @export
nabarx_udd <- function(x, n, i, model, ...) {
  xn <- .recycle_x_n_ann(x, n)
  nEx(xn$x, xn$n, i, model = model, ...) *
    abarx_udd(xn$x + xn$n, i, model = model, ...)
}

# -------------------------------------------------------------------------
# Shared documentation: Woolhouse 2-term approximations
# -------------------------------------------------------------------------

#' Woolhouse 2-term annuity approximations
#'
#' Woolhouse 2-term approximations for Chapter 8 annuity functions.
#'
#' @name annuity_approximations_woolhouse2
#' @aliases ax_m_woolhouse2 adotx_m_woolhouse2
#' @aliases axn_m_woolhouse2 adotxn_m_woolhouse2
#' @aliases nax_m_woolhouse2 nadotx_m_woolhouse2
#' @aliases sxn_m_woolhouse2 sdotxn_m_woolhouse2
#' @aliases abarx_woolhouse2
NULL

#' Woolhouse 2-term approximation to m-thly whole life annuity-immediate
#'
#' @param x Age.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @rdname annuity_approximations_woolhouse2
#' @export
ax_m_woolhouse2 <- function(x, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  ax(x, i, model = model, ...) + (m - 1) / (2 * m)
}

#' Woolhouse 2-term approximation to m-thly whole life annuity-due
#'
#' @rdname annuity_approximations_woolhouse2
#' @export
adotx_m_woolhouse2 <- function(x, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  adotx(x, i, model = model, ...) - (m - 1) / (2 * m)
}

#' Woolhouse 2-term approximation to m-thly deferred annuity-immediate
#'
#' @param n Term.
#' @rdname annuity_approximations_woolhouse2
#' @export
nax_m_woolhouse2 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)

  nax(xn$x, xn$n, i, model = model, ...) +
    (m - 1) / (2 * m) * nEx(xn$x, xn$n, i, model = model, ...)
}

#' Woolhouse 2-term approximation to m-thly deferred annuity-due
#'
#' @rdname annuity_approximations_woolhouse2
#' @export
nadotx_m_woolhouse2 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)

  nadotx(xn$x, xn$n, i, model = model, ...) -
    (m - 1) / (2 * m) * nEx(xn$x, xn$n, i, model = model, ...)
}

#' Woolhouse 2-term approximation to m-thly temporary annuity-immediate
#'
#' @rdname annuity_approximations_woolhouse2
#' @export
axn_m_woolhouse2 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)

  axn(xn$x, xn$n, i, model = model, ...) +
    (m - 1) / (2 * m) * (1 - nEx(xn$x, xn$n, i, model = model, ...))
}

#' Woolhouse 2-term approximation to m-thly temporary annuity-due
#'
#' @rdname annuity_approximations_woolhouse2
#' @export
adotxn_m_woolhouse2 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)

  adotxn(xn$x, xn$n, i, model = model, ...) -
    (m - 1) / (2 * m) * (1 - nEx(xn$x, xn$n, i, model = model, ...))
}

#' Woolhouse 2-term approximation to m-thly temporary accumulated value
#'
#' @rdname annuity_approximations_woolhouse2
#' @export
sxn_m_woolhouse2 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  ex <- nEx(xn$x, xn$n, i, model = model, ...)

  sxn(xn$x, xn$n, i, model = model, ...) +
    (m - 1) / (2 * m) * (1 / ex - 1)
}

#' Woolhouse 2-term approximation to m-thly temporary accumulated value
#'
#' @rdname annuity_approximations_woolhouse2
#' @export
sdotxn_m_woolhouse2 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  xn <- .recycle_x_n_ann(x, n)
  ex <- nEx(xn$x, xn$n, i, model = model, ...)

  sdotxn(xn$x, xn$n, i, model = model, ...) -
    (m - 1) / (2 * m) * (1 / ex - 1)
}

#' Woolhouse 2-term approximation to continuous whole life annuity
#'
#' @rdname annuity_approximations_woolhouse2
#' @export
abarx_woolhouse2 <- function(x, i, model, ...) {
  adotx(x, i, model = model, ...) - 0.5
}

# -------------------------------------------------------------------------
# Shared documentation: Woolhouse 3-term approximations
# -------------------------------------------------------------------------

#' Woolhouse 3-term annuity approximations
#'
#' Woolhouse 3-term approximations for Chapter 8 annuity functions.
#'
#' @name annuity_approximations_woolhouse3
#' @aliases ax_m_woolhouse3 adotx_m_woolhouse3
#' @aliases axn_m_woolhouse3 adotxn_m_woolhouse3
#' @aliases nax_m_woolhouse3 nadotx_m_woolhouse3
#' @aliases abarx_woolhouse3
NULL

#' Woolhouse 3-term approximation to m-thly whole life annuity-immediate
#'
#' @param x Age.
#' @param m Number of payments per year.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @rdname annuity_approximations_woolhouse3
#' @export
ax_m_woolhouse3 <- function(x, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  delta <- .delta_from_i(i)
  x <- as.numeric(x)

  ax(x, i, model = model, ...) +
    (m - 1) / (2 * m) -
    ((m^2 - 1) / (12 * m^2)) * (.mu_x(x, model = model, ...) + delta)
}

#' Woolhouse 3-term approximation to m-thly whole life annuity-due
#'
#' @rdname annuity_approximations_woolhouse3
#' @export
adotx_m_woolhouse3 <- function(x, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  delta <- .delta_from_i(i)
  x <- as.numeric(x)

  adotx(x, i, model = model, ...) -
    (m - 1) / (2 * m) -
    ((m^2 - 1) / (12 * m^2)) * (.mu_x(x, model = model, ...) + delta)
}

#' Woolhouse 3-term approximation to m-thly deferred annuity-immediate
#'
#' @param n Term.
#' @rdname annuity_approximations_woolhouse3
#' @export
nax_m_woolhouse3 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  delta <- .delta_from_i(i)
  xn <- .recycle_x_n_ann(x, n)
  ex <- nEx(xn$x, xn$n, i, model = model, ...)

  nax(xn$x, xn$n, i, model = model, ...) +
    ex * ((m - 1) / (2 * m) -
            ((m^2 - 1) / (12 * m^2)) * (.mu_x(xn$x + xn$n, model = model, ...) + delta))
}

#' Woolhouse 3-term approximation to m-thly deferred annuity-due
#'
#' @rdname annuity_approximations_woolhouse3
#' @export
nadotx_m_woolhouse3 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  delta <- .delta_from_i(i)
  xn <- .recycle_x_n_ann(x, n)
  ex <- nEx(xn$x, xn$n, i, model = model, ...)

  nadotx(xn$x, xn$n, i, model = model, ...) -
    ex * ((m - 1) / (2 * m) +
            ((m^2 - 1) / (12 * m^2)) * (.mu_x(xn$x + xn$n, model = model, ...) + delta))
}

#' Woolhouse 3-term approximation to m-thly temporary annuity-immediate
#'
#' @rdname annuity_approximations_woolhouse3
#' @export
axn_m_woolhouse3 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  delta <- .delta_from_i(i)
  xn <- .recycle_x_n_ann(x, n)
  ex <- nEx(xn$x, xn$n, i, model = model, ...)

  axn(xn$x, xn$n, i, model = model, ...) +
    (m - 1) / (2 * m) * (1 - ex) -
    ((m^2 - 1) / (12 * m^2)) *
    (.mu_x(xn$x, model = model, ...) -
       ex * .mu_x(xn$x + xn$n, model = model, ...) +
       delta * (1 - ex))
}

#' Woolhouse 3-term approximation to m-thly temporary annuity-due
#'
#' @rdname annuity_approximations_woolhouse3
#' @export
adotxn_m_woolhouse3 <- function(x, n, m, i, model, ...) {
  m <- .check_m_scalar_ann(m)
  delta <- .delta_from_i(i)
  xn <- .recycle_x_n_ann(x, n)
  ex <- nEx(xn$x, xn$n, i, model = model, ...)

  adotxn(xn$x, xn$n, i, model = model, ...) -
    (m - 1) / (2 * m) * (1 - ex) -
    ((m^2 - 1) / (12 * m^2)) *
    (.mu_x(xn$x, model = model, ...) -
       ex * .mu_x(xn$x + xn$n, model = model, ...) +
       delta * (1 - ex))
}

#' Woolhouse 3-term approximation to continuous whole life annuity
#'
#' @rdname annuity_approximations_woolhouse3
#' @export
abarx_woolhouse3 <- function(x, i, model, ...) {
  delta <- .delta_from_i(i)
  x <- as.numeric(x)

  adotx(x, i, model = model, ...) -
    0.5 -
    (.mu_x(x, model = model, ...) + delta) / 12
}
