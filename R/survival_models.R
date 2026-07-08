#' Survival models: core survival functions and parametric models
#'
#' This file provides actuarial survival-model utilities.
#'
#' @name survival_models
#' @keywords internal
NULL

# -----------------------------
# Internal helpers
# -----------------------------

.check_model <- function(model) {
  if (length(model) != 1L || !is.character(model) || is.na(model)) {
    stop("model must be a single character string.", call. = FALSE)
  }

  model <- tolower(model)
  ok <- model %in% c("uniform", "exponential", "gompertz", "makeham", "weibull")

  if (!ok) {
    stop(
      "Unknown model='", model,
      "'. Use one of: uniform, exponential, gompertz, makeham, weibull.",
      call. = FALSE
    )
  }

  model
}

.check_params <- function(model, ...) {
  p <- list(...)
  model <- .check_model(model)

  need <- switch(
    model,
    "uniform" = "omega",
    "exponential" = "lambda",
    "gompertz" = c("B", "c"),
    "makeham" = c("A", "B", "c"),
    "weibull" = c("shape", "scale")
  )

  miss <- setdiff(need, names(p))
  if (length(miss) > 0L) {
    stop(
      "Missing parameter(s) for model='", model, "': ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  if (model == "uniform") {
    if (!is.numeric(p$omega) || length(p$omega) != 1L ||
        !is.finite(p$omega) || p$omega <= 0) {
      stop("For model='uniform', omega must be a single positive number.",
           call. = FALSE)
    }
  }

  if (model == "exponential") {
    if (!is.numeric(p$lambda) || length(p$lambda) != 1L ||
        !is.finite(p$lambda) || p$lambda <= 0) {
      stop("For model='exponential', lambda must be a single positive number.",
           call. = FALSE)
    }
  }

  if (model == "gompertz") {
    if (!is.numeric(p$B) || length(p$B) != 1L ||
        !is.finite(p$B) || p$B <= 0) {
      stop("For model='gompertz', B must be a single positive number.",
           call. = FALSE)
    }
    if (!is.numeric(p$c) || length(p$c) != 1L ||
        !is.finite(p$c) || p$c <= 1) {
      stop("For model='gompertz', c must be a single number greater than 1.",
           call. = FALSE)
    }
  }

  if (model == "makeham") {
    if (!is.numeric(p$A) || length(p$A) != 1L || !is.finite(p$A)) {
      stop("For model='makeham', A must be a single finite number.",
           call. = FALSE)
    }
    if (!is.numeric(p$B) || length(p$B) != 1L ||
        !is.finite(p$B) || p$B <= 0) {
      stop("For model='makeham', B must be a single positive number.",
           call. = FALSE)
    }
    if (!is.numeric(p$c) || length(p$c) != 1L ||
        !is.finite(p$c) || p$c <= 1) {
      stop("For model='makeham', c must be a single number greater than 1.",
           call. = FALSE)
    }
    if (p$A <= -p$B) {
      stop("For model='makeham', A must be greater than -B.",
           call. = FALSE)
    }
  }

  if (model == "weibull") {
    if (!is.numeric(p$shape) || length(p$shape) != 1L ||
        !is.finite(p$shape) || p$shape <= 0) {
      stop("For model='weibull', shape must be a single positive number.",
           call. = FALSE)
    }
    if (!is.numeric(p$scale) || length(p$scale) != 1L ||
        !is.finite(p$scale) || p$scale <= 0) {
      stop("For model='weibull', scale must be a single positive number.",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

.recycle_tx <- function(t, x) {
  t <- as.numeric(t)
  x <- as.numeric(x)

  if (length(t) == 0L || length(x) == 0L) {
    stop("t and x must have positive length.", call. = FALSE)
  }

  if (length(t) == 1L && length(x) > 1L) t <- rep(t, length(x))
  if (length(x) == 1L && length(t) > 1L) x <- rep(x, length(t))

  if (length(t) != length(x)) {
    stop("t and x must have the same length, or one must have length 1.",
         call. = FALSE)
  }

  list(t = t, x = x)
}

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

#' Survival function for age-at-failure
#'
#' Computes \eqn{S_0(t)=Pr(T_0>t)}.
#'
#' @param t Numeric vector of times.
#' @param model One of \code{"uniform"}, \code{"exponential"},
#'   \code{"gompertz"}, \code{"makeham"}, or \code{"weibull"}.
#' @param ... Model parameters.
#'
#' @return Numeric vector of survival probabilities.
#' @export
S0 <- function(t, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)

  t <- as.numeric(t)

  if (length(t) == 0L || any(!is.finite(t))) {
    stop("t must contain finite numeric values.", call. = FALSE)
  }
  if (any(t < 0)) {
    stop("t must be nonnegative.", call. = FALSE)
  }

  p <- list(...)

  out <- switch(
    model,
    "uniform" = {
      omega <- p$omega
      ifelse(t <= omega, pmax(omega - t, 0) / omega, 0)
    },
    "exponential" = {
      exp(-p$lambda * t)
    },
    "gompertz" = {
      exp((p$B / log(p$c)) * (1 - p$c^t))
    },
    "makeham" = {
      exp((p$B / log(p$c)) * (1 - p$c^t) - p$A * t)
    },
    "weibull" = {
      exp(-(t / p$scale)^p$shape)
    }
  )

  pmin(pmax(out, 0), 1)
}

#' Distribution functions for age-at-failure
#'
#' @name dist0
#' @aliases F0 f0
#' @inheritParams S0
#' @return Numeric vector.
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

  if (length(t) == 0L || any(!is.finite(t))) {
    stop("t must contain finite numeric values.", call. = FALSE)
  }
  if (any(t < 0)) {
    stop("t must be nonnegative.", call. = FALSE)
  }

  p <- list(...)

  out <- switch(
    model,
    "uniform" = {
      omega <- p$omega
      ifelse(t >= 0 & t < omega, 1 / omega, 0)
    },
    "exponential" = {
      p$lambda * exp(-p$lambda * t)
    },
    "gompertz" = {
      mu <- p$B * p$c^t
      mu * S0(t, model = model, ...)
    },
    "makeham" = {
      mu <- p$A + p$B * p$c^t
      mu * S0(t, model = model, ...)
    },
    "weibull" = {
      (p$shape / p$scale) *
        (t / p$scale)^(p$shape - 1) *
        exp(-(t / p$scale)^p$shape)
    }
  )

  pmax(out, 0)
}

#' Hazard or force of mortality for age-at-failure
#'
#' Computes \eqn{\lambda_0(t)}.
#'
#' @inheritParams S0
#' @return Numeric vector.
#' @export
hazard0 <- function(t, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)

  t <- as.numeric(t)

  if (length(t) == 0L || any(!is.finite(t))) {
    stop("t must contain finite numeric values.", call. = FALSE)
  }
  if (any(t < 0)) {
    stop("t must be nonnegative.", call. = FALSE)
  }

  p <- list(...)

  switch(
    model,
    "uniform" = {
      ifelse(t < p$omega, 1 / (p$omega - t), Inf)
    },
    "exponential" = {
      rep(p$lambda, length(t))
    },
    "gompertz" = {
      p$B * p$c^t
    },
    "makeham" = {
      p$A + p$B * p$c^t
    },
    "weibull" = {
      (p$shape / p$scale) * (t / p$scale)^(p$shape - 1)
    }
  )
}

#' Cumulative hazard for age-at-failure
#'
#' Computes \eqn{\Lambda_0(t)}.
#'
#' @inheritParams S0
#' @return Numeric vector.
#' @export
cumhaz0 <- function(t, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)

  t <- as.numeric(t)

  if (length(t) == 0L || any(!is.finite(t))) {
    stop("t must contain finite numeric values.", call. = FALSE)
  }
  if (any(t < 0)) {
    stop("t must be nonnegative.", call. = FALSE)
  }

  p <- list(...)

  switch(
    model,
    "uniform" = {
      ifelse(t < p$omega, log(p$omega / (p$omega - t)), Inf)
    },
    "exponential" = {
      p$lambda * t
    },
    "gompertz" = {
      (p$B / log(p$c)) * (p$c^t - 1)
    },
    "makeham" = {
      p$A * t + (p$B / log(p$c)) * (p$c^t - 1)
    },
    "weibull" = {
      (t / p$scale)^p$shape
    }
  )
}

#' Conditional survival probability
#'
#' Computes \eqn{{}_t p_x = S_0(x+t)/S_0(x)}.
#'
#' @param t Numeric vector of durations.
#' @param x Numeric vector of ages.
#' @inheritParams S0
#'
#' @return Numeric vector of survival probabilities.
#' @export
tpx <- function(t, x, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)

  tx <- .recycle_tx(t, x)
  t <- tx$t
  x <- tx$x

  if (any(!is.finite(t)) || any(!is.finite(x))) {
    stop("t and x must contain finite numeric values.", call. = FALSE)
  }
  if (any(t < 0) || any(x < 0)) {
    stop("t and x must be nonnegative.", call. = FALSE)
  }

  Sx <- S0(x, model = model, ...)
  Sxt <- S0(x + t, model = model, ...)

  out <- ifelse(Sx <= 0, 0, Sxt / Sx)

  pmin(pmax(out, 0), 1)
}

#' Conditional failure probability
#'
#' Computes \eqn{{}_t q_x = 1 - {}_t p_x}.
#'
#' @inheritParams tpx
#' @return Numeric vector of failure probabilities.
#' @export
tqx <- function(t, x, model, ...) {
  1 - tpx(t = t, x = x, model = model, ...)
}

#' Conditional density
#'
#' Computes \eqn{f_x(t)=f_0(x+t)/S_0(x)}.
#'
#' @inheritParams tpx
#' @return Numeric vector of density values.
#' @export
fx <- function(t, x, model, ...) {
  model <- .check_model(model)
  .check_params(model, ...)

  tx <- .recycle_tx(t, x)
  t <- tx$t
  x <- tx$x

  if (any(!is.finite(t)) || any(!is.finite(x))) {
    stop("t and x must contain finite numeric values.", call. = FALSE)
  }
  if (any(t < 0) || any(x < 0)) {
    stop("t and x must be nonnegative.", call. = FALSE)
  }

  Sx <- S0(x, model = model, ...)
  num <- f0(x + t, model = model, ...)

  out <- ifelse(Sx <= 0, 0, num / Sx)

  pmax(out, 0)
}

#' Complete expectation of life
#'
#' Computes \eqn{\overset{\circ}{e}_x=\int_0^\infty {}_t p_x\,dt}.
#'
#' @param x Numeric vector of ages.
#' @inheritParams S0
#' @param tol Tolerance used to choose a finite integration bound.
#'
#' @return Numeric vector of complete expectations.
#' @export
ex_complete <- function(x, model, ..., tol = 1e-10) {
  model <- .check_model(model)
  .check_params(model, ...)

  x <- as.numeric(x)

  if (length(x) == 0L || any(!is.finite(x))) {
    stop("x must contain finite numeric values.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  p <- list(...)

  if (model == "uniform") {
    return(ifelse(x >= p$omega, 0, (p$omega - x) / 2))
  }

  if (model == "exponential") {
    return(rep(1 / p$lambda, length(x)))
  }

  sapply(x, function(xx) {
    if (S0(xx, model = model, ...) <= 0) return(0)

    u_max <- .find_u_max(xx, model = model, ..., tol = tol)
    integrand <- function(u) tpx(u, xx, model = model, ...)

    stats::integrate(integrand, lower = 0, upper = u_max,
                     rel.tol = 1e-8)$value
  })
}

#' Curtate expectation of life
#'
#' Computes \eqn{e_x = E[K_x] = \sum_{k=1}^{\infty} {}_k p_x}.
#'
#' @param x Numeric vector of ages.
#' @inheritParams S0
#' @param k_max Maximum integer duration to sum to.
#' @param tol Stop early if the summand is smaller than this tolerance for
#'   several consecutive steps.
#'
#' @return Numeric vector of curtate expectations.
#' @export
ex_curtate <- function(x, model, ..., k_max = 5000, tol = 1e-12) {
  model <- .check_model(model)
  .check_params(model, ...)

  x <- as.numeric(x)

  if (length(x) == 0L || any(!is.finite(x))) {
    stop("x must contain finite numeric values.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  sapply(x, function(xx) {
    if (S0(xx, model = model, ...) <= 0) return(0)

    acc <- 0
    small_run <- 0

    for (k in seq_len(k_max)) {
      term <- tpx(k, xx, model = model, ...)
      acc <- acc + term

      if (term < tol) {
        small_run <- small_run + 1L
        if (small_run >= 10L) break
      } else {
        small_run <- 0L
      }
    }

    acc
  })
}
