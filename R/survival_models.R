#' Survival models (Chapter 5): core survival functions + parametric models
#'
#' This file provides actuarial survival-model utilities for Chapter 5.
#'
#' ## Age-at-failure random variable (T0)
#' - `S0(t, ...)` : survival function  S0(t) = Pr(T0 > t)
#' - `F0(t, ...)` : CDF              F0(t) = Pr(T0 <= t) = 1 - S0(t)
#' - `f0(t, ...)` : density          f0(t) = d/dt F0(t)
#' - `hazard0(t, ...)` : hazard      lambda0(t) = f0(t) / S0(t)
#' - `cumhaz0(t, ...)` : cum hazard  Lambda0(t) = \eqn{\int_0^t \lambda_0(y) dy}
#'
#' ## Conditional (future lifetime) random variable (Tx), given alive at age x
#' - `tpx(t, x, ...)` :  \eqn{{}_t p_x = S0(x+t)/S0(x)}
#' - `tqx(t, x, ...)` :  \eqn{{}_t q_x = 1 - {}_t p_x}
#' - `fx(t, x, ...)`  :  \eqn{f_x(t) = f0(x+t)/S0(x)}
#'
#' ## Expectations
#' - `ex_complete(x, ...)` : \eqn{\overset{\circ}{e}_x = \int_0^\infty {}_t p_x\,dt}
#' - `ex_curtate(x, ...)`  : \eqn{e_x = E[K_x] = \sum_{k\ge 1} {}_k p_x} (truncated sum)
#'
#' ## Models supported
#' - `"uniform"`     (de Moivre)      parameter: `omega`
#' - `"exponential"` (constant force) parameter: `lambda`
#' - `"gompertz"`                  parameters: `B`, `c`
#' - `"makeham"`                   parameters: `A`, `B`, `c`
#' - `"weibull"`                   parameters: `shape`, `scale`
#'
#' All functions are vectorized over `t` (and over `x` where applicable).
#'
#' @name survival_models
#' @keywords internal
NULL

# -----------------------------
# Internal helpers (not exported)
# -----------------------------

#' @noRd
.check_model <- function(model) {
  model <- tolower(model)
  ok <- model %in% c("uniform", "exponential", "gompertz", "makeham", "weibull")
  if (!ok) stop("Unknown model='", model, "'. Use one of: uniform, exponential, gompertz, makeham, weibull.")
  model
}

#' @noRd
.check_params <- function(model, ...) {
  p <- list(...)
  model <- .check_model(model)

  need <- switch(
    model,
    "uniform"     = c("omega"),
    "exponential" = c("lambda"),
    "gompertz"    = c("B", "c"),
    "makeham"     = c("A", "B", "c"),
    "weibull"     = c("shape", "scale")
  )

  miss <- setdiff(need, names(p))
  if (length(miss) > 0) {
    stop("Missing parameter(s) for model='", model, "': ", paste(miss, collapse = ", "))
  }

  if (model == "uniform") {
    if (!is.numeric(p$omega) || length(p$omega) != 1 || !is.finite(p$omega) || p$omega <= 0) {
      stop("For model='uniform', omega must be a single positive number.")
    }
  }

  if (model == "exponential") {
    if (!is.numeric(p$lambda) || length(p$lambda) != 1 || !is.finite(p$lambda) || p$lambda <= 0) {
      stop("For model='exponential', lambda must be a single positive number.")
    }
  }

  if (model == "gompertz") {
    if (!is.numeric(p$B) || length(p$B) != 1 || !is.finite(p$B) || p$B <= 0) stop("For gompertz, B must be > 0.")
    if (!is.numeric(p$c) || length(p$c) != 1 || !is.finite(p$c) || p$c <= 1) stop("For gompertz, c must be > 1.")
  }

  if (model == "makeham") {
    if (!is.numeric(p$A) || length(p$A) != 1 || !is.finite(p$A)) stop("For makeham, A must be finite.")
    if (!is.numeric(p$B) || length(p$B) != 1 || !is.finite(p$B) || p$B <= 0) stop("For makeham, B must be > 0.")
    if (!is.numeric(p$c) || length(p$c) != 1 || !is.finite(p$c) || p$c <= 1) stop("For makeham, c must be > 1.")
    if (p$A <= -p$B) stop("For makeham, need A > -B (so force of mortality can be positive at t=0).")
  }

  if (model == "weibull") {
    if (!is.numeric(p$shape) || length(p$shape) != 1 || !is.finite(p$shape) || p$shape <= 0) stop("For weibull, shape must be > 0.")
    if (!is.numeric(p$scale) || length(p$scale) != 1 || !is.finite(p$scale) || p$scale <= 0) stop("For weibull, scale must be > 0.")
  }

  invisible(TRUE)
}

#' @noRd
.recycle_tx <- function(t, x) {
  t <- as.numeric(t)
  x <- as.numeric(x)

  if (length(t) == 0L || length(x) == 0L) stop("t and x must have positive length.")
  if (length(t) == 1L && length(x) > 1L) t <- rep(t, length(x))
  if (length(x) == 1L && length(t) > 1L) x <- rep(x, length(t))

  if (length(t) != length(x)) {
    stop("t and x must have the same length, or one of them must have length 1.")
  }

  list(t = t, x = x)
}

#' @noRd
.find_u_max <- function(x, model, ..., tol = 1e-10, u_max_cap = 5000) {
  u <- 10
  tpx_u <- tpx(u, x, model = model, ...)
  if (is.na(tpx_u) || tpx_u < tol) return(u)

  while (u < u_max_cap) {
    u <- u * 2
    tpx_u <- tpx(u, x, model = model, ...)
    if (is.na(tpx_u) || tpx_u < tol) return(u)
  }
  u_max_cap
}

# -----------------------------
# Exported functions
# -----------------------------

#' Survival function for age-at-failure T0
#'
#' Computes the survival distribution function \eqn{S_0(t)=Pr(T_0>t)}.
#'
#' Supported models (Chapter 5): uniform (de Moivre), exponential,
#' Gompertz, Makeham, Weibull.
#'
#' @param t Numeric vector of times (\eqn{t \ge 0}).
#' @param model One of \code{"uniform"}, \code{"exponential"}, \code{"gompertz"},
#'   \code{"makeham"}, \code{"weibull"}.
#' @param ... Model parameters:
#'   \itemize{
#'     \item uniform: \code{omega}
#'     \item exponential: \code{lambda}
#'     \item gompertz: \code{B}, \code{c}
#'     \item makeham: \code{A}, \code{B}, \code{c}
#'     \item weibull: \code{shape}, \code{scale}
#'   }
#' @return Numeric vector of survival probabilities in \eqn{[0,1]}.
#' @export
S0 <- function(t, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)
  t <- as.numeric(t)
  if (any(!is.finite(t))) stop("t must be finite.")
  if (any(t < 0)) stop("t must be >= 0.")

  p <- list(...)
  out <- switch(
    model,
    "uniform" = {
      omega <- p$omega
      ifelse(t <= omega, pmax(omega - t, 0) / omega, 0)
    },
    "exponential" = {
      lambda <- p$lambda
      exp(-lambda * t)
    },
    "gompertz" = {
      B <- p$B; c <- p$c
      exp((B / log(c)) * (1 - c^t))
    },
    "makeham" = {
      A <- p$A; B <- p$B; c <- p$c
      exp((B / log(c)) * (1 - c^t) - A * t)
    },
    "weibull" = {
      shape <- p$shape; scale <- p$scale
      exp(- (t / scale)^shape)
    }
  )
  pmin(pmax(out, 0), 1)
}

#' Distribution functions for age-at-failure T0
#'
#' Convenience functions for the CDF and density of \eqn{T_0}.
#'
#' \itemize{
#'   \item \code{F0(t)} computes \eqn{F_0(t)=Pr(T_0 \le t)=1-S_0(t)}
#'   \item \code{f0(t)} computes \eqn{f_0(t)=\frac{d}{dt}F_0(t)}
#' }
#'
#' @name dist0
#' @aliases F0 f0
#' @inheritParams S0
#' @return Numeric vector. For \code{F0}: CDF values in \eqn{[0,1]}. For \code{f0}: density values (\eqn{\ge 0}).
NULL

#' @rdname dist0
#' @export
F0 <- function(t, model, ...) {
  1 - S0(t, model = model, ...)
}

#' @rdname dist0
#' @export
f0 <- function(t, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)
  t <- as.numeric(t)
  if (any(!is.finite(t))) stop("t must be finite.")
  if (any(t < 0)) stop("t must be >= 0.")

  p <- list(...)
  out <- switch(
    model,
    "uniform" = {
      omega <- p$omega
      ifelse(t > 0 & t < omega, 1 / omega, 0)
    },
    "exponential" = {
      lambda <- p$lambda
      lambda * exp(-lambda * t)
    },
    "gompertz" = {
      B <- p$B; c <- p$c
      mu <- B * c^t
      mu * S0(t, model = model, ...)
    },
    "makeham" = {
      A <- p$A; B <- p$B; c <- p$c
      mu <- A + B * c^t
      mu * S0(t, model = model, ...)
    },
    "weibull" = {
      shape <- p$shape; scale <- p$scale
      (shape / scale) * (t / scale)^(shape - 1) * exp(- (t / scale)^shape)
    }
  )
  pmax(out, 0)
}

#' Hazard / force for age-at-failure T0
#'
#' Computes \eqn{\lambda_0(t)=f_0(t)/S_0(t)}.
#'
#' @inheritParams S0
#' @return Numeric vector of hazard/force values.
#' @export
hazard0 <- function(t, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)
  t <- as.numeric(t)
  if (any(!is.finite(t))) stop("t must be finite.")
  if (any(t < 0)) stop("t must be >= 0.")

  p <- list(...)
  out <- switch(
    model,
    "uniform" = {
      omega <- p$omega
      ifelse(t < omega, 1 / (omega - t), Inf)
    },
    "exponential" = {
      lambda <- p$lambda
      rep(lambda, length(t))
    },
    "gompertz" = {
      B <- p$B; c <- p$c
      B * c^t
    },
    "makeham" = {
      A <- p$A; B <- p$B; c <- p$c
      A + B * c^t
    },
    "weibull" = {
      shape <- p$shape; scale <- p$scale
      (shape / scale) * (t / scale)^(shape - 1)
    }
  )
  out
}

#' Cumulative hazard for age-at-failure T0
#'
#' Computes \eqn{\Lambda_0(t)=\int_0^t \lambda_0(y)\,dy}.
#'
#' @inheritParams S0
#' @return Numeric vector of cumulative hazard values.
#' @export
cumhaz0 <- function(t, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)
  t <- as.numeric(t)
  if (any(!is.finite(t))) stop("t must be finite.")
  if (any(t < 0)) stop("t must be >= 0.")

  p <- list(...)
  out <- switch(
    model,
    "uniform" = {
      omega <- p$omega
      ifelse(t < omega, log(omega / (omega - t)), Inf)
    },
    "exponential" = {
      lambda <- p$lambda
      lambda * t
    },
    "gompertz" = {
      B <- p$B; c <- p$c
      (B / log(c)) * (c^t - 1)
    },
    "makeham" = {
      A <- p$A; B <- p$B; c <- p$c
      A * t + (B / log(c)) * (c^t - 1)
    },
    "weibull" = {
      shape <- p$shape; scale <- p$scale
      (t / scale)^shape
    }
  )
  out
}

#' Conditional survival probability for Tx
#'
#' Computes \eqn{{}_t p_x = S_0(x+t)/S_0(x)}.
#'
#' Vectorization rule:
#' - If \code{t} and \code{x} are the same length, values are computed elementwise.
#' - If one of \code{t} or \code{x} has length 1, it is recycled to match the other.
#'
#' @param t Numeric vector of durations (\eqn{t \ge 0}).
#' @param x Numeric vector of ages (\eqn{x \ge 0}).
#' @inheritParams S0
#' @return Numeric vector in \eqn{[0,1]}.
#' @export
tpx <- function(t, x, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)

  tx <- .recycle_tx(t, x)
  t <- tx$t
  x <- tx$x

  if (any(!is.finite(t))) stop("t must be finite.")
  if (any(!is.finite(x))) stop("x must be finite.")
  if (any(t < 0)) stop("t must be >= 0.")
  if (any(x < 0)) stop("x must be >= 0.")

  Sx  <- S0(x, model = model, ...)
  Sxt <- S0(x + t, model = model, ...)

  out <- ifelse(Sx <= 0, 0, Sxt / Sx)
  pmin(pmax(out, 0), 1)
}

#' Conditional failure probability for Tx
#'
#' Computes \eqn{{}_t q_x = 1 - {}_t p_x}.
#'
#' @inheritParams tpx
#' @return Numeric vector in \eqn{[0,1]}.
#' @export
tqx <- function(t, x, model, ...) {
  1 - tpx(t, x, model = model, ...)
}

#' Conditional density for Tx
#'
#' Computes \eqn{f_x(t)=f_0(x+t)/S_0(x)}.
#'
#' Vectorization rule:
#' - If \code{t} and \code{x} are the same length, values are computed elementwise.
#' - If one of \code{t} or \code{x} has length 1, it is recycled to match the other.
#'
#' @inheritParams tpx
#' @return Numeric vector of densities (\eqn{\ge 0}).
#' @export
fx <- function(t, x, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)

  tx <- .recycle_tx(t, x)
  t <- tx$t
  x <- tx$x

  if (any(!is.finite(t))) stop("t must be finite.")
  if (any(!is.finite(x))) stop("x must be finite.")
  if (any(t < 0)) stop("t must be >= 0.")
  if (any(x < 0)) stop("x must be >= 0.")

  Sx <- S0(x, model = model, ...)
  num <- f0(x + t, model = model, ...)
  out <- ifelse(Sx <= 0, 0, num / Sx)
  pmax(out, 0)
}

#' Complete expectation of life at age x
#'
#' Computes \eqn{\overset{\circ}{e}_x=\int_0^\infty {}_t p_x \, dt}.
#' Uses closed form for uniform and exponential; numeric integration otherwise.
#'
#' @param x Numeric vector of ages (\eqn{x \ge 0}).
#' @inheritParams S0
#' @param tol Tolerance used to choose a finite integration bound (numeric).
#' @return Numeric vector of complete expectations.
#' @export
ex_complete <- function(x, model, ..., tol = 1e-10) {
  model <- .check_model(model)
  .check_params(model, ...)
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("x must be finite.")
  if (any(x < 0)) stop("x must be >= 0.")

  p <- list(...)

  if (model == "uniform") {
    omega <- p$omega
    return(ifelse(x >= omega, 0, (omega - x) / 2))
  }

  if (model == "exponential") {
    lambda <- p$lambda
    return(rep(1 / lambda, length(x)))
  }

  sapply(x, function(xx) {
    if (S0(xx, model = model, ...) <= 0) return(0)

    u_max <- .find_u_max(xx, model = model, ..., tol = tol)
    integrand <- function(u) tpx(u, xx, model = model, ...)

    stats::integrate(integrand, lower = 0, upper = u_max, rel.tol = 1e-8)$value
  })
}

#' Curtate expectation of life at age x
#'
#' Computes \eqn{e_x = E[K_x] = \sum_{k=1}^\infty {}_k p_x} with truncation.
#'
#' @param x Numeric vector of ages (\eqn{x \ge 0}).
#' @inheritParams S0
#' @param k_max Maximum integer duration to sum to.
#' @param tol Stop early if the summand is < tol for several steps.
#' @return Numeric vector of curtate expectations.
#' @export
ex_curtate <- function(x, model, ..., k_max = 5000, tol = 1e-12) {
  model <- .check_model(model)
  .check_params(model, ...)
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("x must be finite.")
  if (any(x < 0)) stop("x must be >= 0.")

  sapply(x, function(xx) {
    if (S0(xx, model = model, ...) <= 0) return(0)

    acc <- 0
    small_run <- 0
    for (k in 1:k_max) {
      term <- tpx(k, xx, model = model, ...)
      acc <- acc + term

      if (term < tol) {
        small_run <- small_run + 1
        if (small_run >= 10) break
      } else {
        small_run <- 0
      }
    }
    acc
  })
}
