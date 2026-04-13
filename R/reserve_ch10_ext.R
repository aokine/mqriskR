#' Extended reserve functions for Chapter 10
#'
#' Additional Chapter 10 reserve functions:
#' retrospective reserves, deferred insurance reserves, annuity reserves,
#' continuous-time gain/loss helpers, and Euler/Thiele approximation helpers.
#'
#' @name reserve_ch10_ext
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_prob_ch10 <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) == 0L || any(!is.finite(x)) || any(x < 0) || any(x > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }
  x
}

#' @noRd
.check_interest_ch10 <- function(i, name = "i") {
  i <- as.numeric(i)
  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop(name, " must contain finite values greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.common_length_ch10 <- function(...) {
  lens <- vapply(list(...), length, integer(1))
  max(lens)
}

#' @noRd
.recycle_common_ch10 <- function(...) {
  xs <- list(...)
  n <- .common_length_ch10(...)
  out <- vector("list", length(xs))

  for (j in seq_along(xs)) {
    if (!length(xs[[j]]) %in% c(1L, n)) {
      stop("Inputs must have common length or length 1.", call. = FALSE)
    }
    out[[j]] <- rep(xs[[j]], length.out = n)
  }

  out
}

# -------------------------------------------------------------------------
# Retrospective reserves
# -------------------------------------------------------------------------

#' Whole life net level premium reserve by retrospective method
#'
#' Computes the Chapter 10 retrospective reserve
#' \eqn{{}_tV_x = P_x \ddot{s}_{x:\overline{t}|} - {}_tk_x}.
#'
#' @param x Issue age.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tVx_ret(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVx_ret <- function(x, t, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")

  xt <- .recycle2_ch10(x, t, "x", "t")
  x <- xt$a
  t <- xt$b

  sapply(seq_along(x), function(j) {
    Px(x = x[j], i = i, model = model, ...) *
      sdotxn(x = x[j], n = t[j], i = i, model = model, ...) -
      Axn1(x = x[j], n = t[j], i = i, model = model, ...) /
      nEx(x = x[j], n = t[j], i = i, model = model, ...)
  })
}

#' Endowment insurance reserve by retrospective method
#'
#' Computes the Chapter 10 retrospective reserve
#' \eqn{{}_tV_{x:\overline{n}|} = P_{x:\overline{n}|} \ddot{s}_{x:\overline{t}|} - {}_tk_x}
#' for \eqn{t \le n}.
#'
#' @param x Issue age.
#' @param n Term in years.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tVxn_ret(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVxn_ret <- function(x, n, t, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  .check_t_le_n(t, n, allow_equal = TRUE)

  sapply(seq_along(x), function(j) {
    if (t[j] == n[j]) return(1)

    Pxn(x = x[j], n = n[j], i = i, model = model, ...) *
      sdotxn(x = x[j], n = t[j], i = i, model = model, ...) -
      Axn1(x = x[j], n = t[j], i = i, model = model, ...) /
      nEx(x = x[j], n = t[j], i = i, model = model, ...)
  })
}

#' Term insurance reserve by retrospective method
#'
#' Computes the retrospective term reserve for \eqn{t \le n}.
#'
#' @inheritParams tVxn_ret
#' @return Numeric vector.
#' @examples
#' tVxn1_ret(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVxn1_ret <- function(x, n, t, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  .check_t_le_n(t, n, allow_equal = TRUE)

  sapply(seq_along(x), function(j) {
    if (t[j] == n[j]) return(0)

    Pxn1(x = x[j], n = n[j], i = i, model = model, ...) *
      sdotxn(x = x[j], n = t[j], i = i, model = model, ...) -
      Axn1(x = x[j], n = t[j], i = i, model = model, ...) /
      nEx(x = x[j], n = t[j], i = i, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# Deferred insurance reserves
# -------------------------------------------------------------------------

#' Deferred insurance reserve functions
#'
#' Chapter 10 reserve functions for deferred insurance contracts.
#'
#' \code{tVnAx()} computes the reserve for an \eqn{n}-year deferred insurance
#' funded over the deferred period.
#'
#' \code{htVnAx()} computes the reserve when premiums are limited to the first
#' \eqn{h} years, with \eqn{h \le n}.
#'
#' @name reserve_deferred_insurance_ch10
#' @aliases tVnAx htVnAx
NULL

#' @rdname reserve_deferred_insurance_ch10
#' @param x Issue age.
#' @param n Deferral period.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tVnAx(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVnAx <- function(x, n, t, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  sapply(seq_along(x), function(j) {
    if (t[j] < n[j]) {
      nAx(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...) -
        PnAx(x = x[j], n = n[j], i = i, model = model, ...) *
        adotxn(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...)
    } else {
      Ax(x = x[j] + t[j], i = i, model = model, ...)
    }
  })
}

#' @rdname reserve_deferred_insurance_ch10
#' @param h Premium-paying period.
#' @inheritParams tVnAx
#'
#' @return Numeric vector.
#' @examples
#' htVnAx(40, n = 20, h = 10, t = 5, i = 0.05, model = "uniform", omega = 100)
#' @export
htVnAx <- function(x, n, h, t, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  h <- .check_nonneg_integerish(h, "h")
  t <- .check_nonneg_integerish(t, "t")

  xnht <- .recycle_common_ch10(x, n, h, t)
  x <- xnht[[1]]
  n <- xnht[[2]]
  h <- xnht[[3]]
  t <- xnht[[4]]

  if (any(h > n)) {
    stop("h must satisfy h <= n for this function.", call. = FALSE)
  }

  sapply(seq_along(x), function(j) {
    if (t[j] < h[j]) {
      nAx(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...) -
        tPnAx(x = x[j], n = n[j], t = h[j], i = i, model = model, ...) *
        adotxn(x = x[j] + t[j], n = h[j] - t[j], i = i, model = model, ...)
    } else if (t[j] < n[j]) {
      nAx(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...)
    } else {
      Ax(x = x[j] + t[j], i = i, model = model, ...)
    }
  })
}

# -------------------------------------------------------------------------
# Annuity reserves
# -------------------------------------------------------------------------

#' Deferred annuity-due premium
#'
#' Computes
#' \eqn{P({}_{n|}\ddot{a}_x) = {}_{n|}\ddot{a}_x / \ddot{a}_{x:\overline{n}|}}.
#'
#' @param x Issue age.
#' @param n Deferral period.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' PnAdotx(40, n = 20, i = 0.05, model = "uniform", omega = 100)
#' @export
PnAdotx <- function(x, n, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")

  xn <- .recycle2_ch10(x, n, "x", "n")
  x <- xn$a
  n <- xn$b

  sapply(seq_along(x), function(j) {
    nadotx(x = x[j], n = n[j], i = i, model = model, ...) /
      adotxn(x = x[j], n = n[j], i = i, model = model, ...)
  })
}

#' Deferred annuity-due reserve
#'
#' Computes the Chapter 10 reserve
#' \eqn{{}_tV({}_{n|}\ddot{a}_x)}
#' for \eqn{t < n}.
#'
#' @param x Issue age.
#' @param n Deferral period.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tVnAdotx(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVnAdotx <- function(x, n, t, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  if (any(t >= n)) {
    stop("t must satisfy t < n for deferred annuity reserve.", call. = FALSE)
  }

  sapply(seq_along(x), function(j) {
    nadotx(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...) -
      PnAdotx(x = x[j], n = n[j], i = i, model = model, ...) *
      adotxn(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...)
  })
}

#' Deferred annuity-immediate premium
#'
#' Computes
#' \eqn{P({}_{n|}a_x) = {}_{n|}a_x / \ddot{a}_{x:\overline{n}|}}.
#'
#' @inheritParams PnAdotx
#' @return Numeric vector.
#' @examples
#' Pnax(40, n = 20, i = 0.05, model = "uniform", omega = 100)
#' @export
Pnax <- function(x, n, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")

  xn <- .recycle2_ch10(x, n, "x", "n")
  x <- xn$a
  n <- xn$b

  sapply(seq_along(x), function(j) {
    nax(x = x[j], n = n[j], i = i, model = model, ...) /
      adotxn(x = x[j], n = n[j], i = i, model = model, ...)
  })
}

#' Deferred annuity-immediate reserve
#'
#' Computes the reserve for an n-year deferred annuity-immediate
#' for \eqn{t < n}.
#'
#' @inheritParams tVnAdotx
#' @return Numeric vector.
#' @examples
#' tVnax(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVnax <- function(x, n, t, i, model, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  if (any(t >= n)) {
    stop("t must satisfy t < n for deferred annuity reserve.", call. = FALSE)
  }

  sapply(seq_along(x), function(j) {
    nax(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...) -
      Pnax(x = x[j], n = n[j], i = i, model = model, ...) *
      adotxn(x = x[j] + t[j], n = n[j] - t[j], i = i, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# Continuous-time gain / loss helpers
# -------------------------------------------------------------------------

#' Total gain for a continuous-style one-step recursion
#'
#' @param Vt Reserve at time t.
#' @param Vt1 Reserve at time t+h.
#' @param P Premium rate.
#' @param delta_actual Actual force of interest.
#' @param p_actual Actual survival probability over the step.
#' @param benefit Benefit paid at start of step. Default 0.
#' @param h Step length. Default 1.
#'
#' @return Numeric vector.
#' @examples
#' GT_cont(Vt = 10, Vt1 = 11, P = 1, delta_actual = 0.05, p_actual = 0.99)
#' @export
GT_cont <- function(Vt, Vt1, P, delta_actual, p_actual, benefit = 0, h = 1) {
  vals <- .recycle_common_ch10(
    .check_nonneg_numeric(Vt, "Vt"),
    .check_nonneg_numeric(Vt1, "Vt1"),
    .check_nonneg_numeric(P, "P"),
    .check_nonneg_numeric(delta_actual, "delta_actual"),
    .check_prob_ch10(p_actual, "p_actual"),
    .check_nonneg_numeric(benefit, "benefit"),
    .check_nonneg_numeric(h, "h")
  )

  Vt <- vals[[1]]
  Vt1 <- vals[[2]]
  P <- vals[[3]]
  delta_actual <- vals[[4]]
  p_actual <- vals[[5]]
  benefit <- vals[[6]]
  h <- vals[[7]]

  (Vt - benefit + P * h) * exp(delta_actual * h) - p_actual * Vt1
}

#' Mortality gain helper for continuous-style recursion
#'
#' @param delta_assumed Assumed force of interest.
#' @inheritParams GT_cont
#'
#' @return Numeric vector.
#' @examples
#' GM_cont(Vt = 10, Vt1 = 11, P = 1, delta_assumed = 0.05, p_actual = 0.99)
#' @export
GM_cont <- function(Vt, Vt1, P, delta_assumed, p_actual, benefit = 0, h = 1) {
  GT_cont(
    Vt = Vt, Vt1 = Vt1, P = P,
    delta_actual = delta_assumed,
    p_actual = p_actual,
    benefit = benefit, h = h
  )
}

#' Interest gain helper for continuous-style recursion
#'
#' @param delta_actual Actual force of interest.
#' @param p_assumed Assumed survival probability over the step.
#' @inheritParams GT_cont
#'
#' @return Numeric vector.
#' @examples
#' GI_cont(Vt = 10, Vt1 = 11, P = 1, delta_actual = 0.05, p_assumed = 0.99)
#' @export
GI_cont <- function(Vt, Vt1, P, delta_actual, p_assumed, benefit = 0, h = 1) {
  GT_cont(
    Vt = Vt, Vt1 = Vt1, P = P,
    delta_actual = delta_actual,
    p_actual = p_assumed,
    benefit = benefit, h = h
  )
}

# -------------------------------------------------------------------------
# Euler / Thiele helpers
# -------------------------------------------------------------------------

#' One backward Euler-style Thiele step
#'
#' Approximates the reserve at time \eqn{t} from a known reserve at
#' time \eqn{t+h}.
#'
#' @param V_next Reserve at time t+h.
#' @param P Premium rate.
#' @param delta Force of interest.
#' @param mu Force of mortality at time t.
#' @param benefit Benefit amount. Defaults to 1.
#' @param h Step size.
#'
#' @return Numeric vector.
#' @examples
#' thiele_backward_step(V_next = 1000, P = 26.96, delta = 0.058, mu = 0.002, benefit = 1000, h = 1)
#' @export
thiele_backward_step <- function(V_next, P, delta, mu, benefit = 1, h = 1) {
  vals <- .recycle_common_ch10(
    .check_nonneg_numeric(V_next, "V_next"),
    .check_nonneg_numeric(P, "P"),
    .check_nonneg_numeric(delta, "delta"),
    .check_nonneg_numeric(mu, "mu"),
    .check_nonneg_numeric(benefit, "benefit"),
    .check_nonneg_numeric(h, "h")
  )

  V_next <- vals[[1]]
  P <- vals[[2]]
  delta <- vals[[3]]
  mu <- vals[[4]]
  benefit <- vals[[5]]
  h <- vals[[6]]

  (V_next - h * P + h * mu * benefit) / (1 + h * (delta + mu))
}

#' Reserve derivative from Thiele's equation
#'
#' Computes
#' \eqn{dV/dt = P + \delta V - \mu(B - V)}.
#'
#' @param V Reserve at time t.
#' @param P Premium rate.
#' @param delta Force of interest.
#' @param mu Force of mortality.
#' @param benefit Benefit amount. Defaults to 1.
#'
#' @return Numeric vector.
#' @examples
#' thiele_dVdt(V = 900, P = 25, delta = 0.05, mu = 0.002, benefit = 1000)
#' @export
thiele_dVdt <- function(V, P, delta, mu, benefit = 1) {
  vals <- .recycle_common_ch10(
    .check_nonneg_numeric(V, "V"),
    .check_nonneg_numeric(P, "P"),
    .check_nonneg_numeric(delta, "delta"),
    .check_nonneg_numeric(mu, "mu"),
    .check_nonneg_numeric(benefit, "benefit")
  )

  V <- vals[[1]]
  P <- vals[[2]]
  delta <- vals[[3]]
  mu <- vals[[4]]
  benefit <- vals[[5]]

  P + delta * V - mu * (benefit - V)
}

#' Backward Euler reserve path from maturity
#'
#' Starting from a terminal reserve value at time T, computes reserves
#' backward on a grid using the backward Euler-style Thiele step.
#'
#' @param times Vector of times in increasing order.
#' @param V_terminal Reserve at the final time.
#' @param P Premium rate, scalar or vector of length length(times)-1.
#' @param delta Force of interest, scalar or vector of length length(times)-1.
#' @param mu Force of mortality, scalar or vector of length length(times)-1.
#' @param benefit Benefit amount, scalar or vector of length length(times)-1.
#'
#' @return Numeric vector of reserve values on the grid.
#' @examples
#' times <- seq(19, 20, by = 0.25)
#' thiele_backward_path(times, V_terminal = 1000, P = 26.96, delta = 0.058, mu = 0.002, benefit = 1000)
#' @export
thiele_backward_path <- function(times, V_terminal, P, delta, mu, benefit = 1) {
  times <- as.numeric(times)
  if (length(times) < 2L || any(!is.finite(times))) {
    stop("times must be a finite numeric vector of length at least 2.", call. = FALSE)
  }
  if (is.unsorted(times, strictly = TRUE)) {
    stop("times must be strictly increasing.", call. = FALSE)
  }

  n_steps <- length(times) - 1L

  vals <- .recycle_common_ch10(
    .check_nonneg_numeric(P, "P"),
    .check_nonneg_numeric(delta, "delta"),
    .check_nonneg_numeric(mu, "mu"),
    .check_nonneg_numeric(benefit, "benefit")
  )

  P <- vals[[1]]
  delta <- vals[[2]]
  mu <- vals[[3]]
  benefit <- vals[[4]]

  if (!length(P) %in% c(1L, n_steps)) {
    stop("P must have length 1 or length(length(times) - 1).", call. = FALSE)
  }
  if (!length(delta) %in% c(1L, n_steps)) {
    stop("delta must have length 1 or length(length(times) - 1).", call. = FALSE)
  }
  if (!length(mu) %in% c(1L, n_steps)) {
    stop("mu must have length 1 or length(length(times) - 1).", call. = FALSE)
  }
  if (!length(benefit) %in% c(1L, n_steps)) {
    stop("benefit must have length 1 or length(length(times) - 1).", call. = FALSE)
  }

  P <- rep(P, length.out = n_steps)
  delta <- rep(delta, length.out = n_steps)
  mu <- rep(mu, length.out = n_steps)
  benefit <- rep(benefit, length.out = n_steps)

  out <- numeric(length(times))
  out[length(times)] <- as.numeric(V_terminal)

  for (k in seq(n_steps, 1L, by = -1L)) {
    h <- times[k + 1L] - times[k]
    out[k] <- thiele_backward_step(
      V_next = out[k + 1L],
      P = P[k],
      delta = delta[k],
      mu = mu[k],
      benefit = benefit[k],
      h = h
    )
  }

  out
}
