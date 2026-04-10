# =========================================================
# Chapter 14 helpers for mqriskR
# Multiple-Decrement Models (Applications)
# =========================================================

#' Discrete multiple-decrement insurance APV \eqn{A_{x}^{(j)}}
#'
#' Computes the actuarial present value of a benefit payable at the end of the
#' year of decrement if decrement occurs by Cause \eqn{j}, matching Equation
#' (14.3b) in Chapter 14.
#'
#' The function evaluates
#' \deqn{
#' A_{x}^{(j)} = \sum_{k=0}^{n-1} v^{k+1} {}_{k}p_{x}^{(\tau)} q_{x+k}^{(j)}
#' }
#' with an optional benefit amount multiplier.
#'
#' @param qj Numeric vector of conditional probabilities
#'   \eqn{q_{x+k}^{(j)}} for Cause \eqn{j}.
#' @param ptau Numeric vector of survival probabilities
#'   \eqn{{}_{k}p_{x}^{(\tau)}} of remaining in force to duration \eqn{k}.
#' @param i Effective annual interest rate.
#' @param benefit Benefit amount payable on decrement by Cause \eqn{j}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' q1 <- c(.02, .02, .02, .02, .02)
#' q2 <- c(.03, .04, .05, .06, .00)
#' q3 <- c(.00, .00, .00, .00, .98)
#' qtau <- q1 + q2 + q3
#'
#' ptau <- numeric(length(qtau))
#' ptau[1] <- 1
#' for (k in 2:length(qtau)) {
#'   ptau[k] <- prod(1 - qtau[1:(k - 1)])
#' }
#'
#' Axj_md(qj = q1, ptau = ptau, i = 0.06, benefit = 1000)
#'
#' @export
Axj_md <- function(qj, ptau, i, benefit = 1) {
  if (!is.numeric(qj) || !is.numeric(ptau)) {
    stop("qj and ptau must be numeric vectors.", call. = FALSE)
  }
  if (length(qj) != length(ptau)) {
    stop("qj and ptau must have the same length.", call. = FALSE)
  }
  if (!is.numeric(i) || length(i) != 1 || !is.finite(i)) {
    stop("i must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(benefit) || length(benefit) != 1 || !is.finite(benefit)) {
    stop("benefit must be a finite numeric scalar.", call. = FALSE)
  }

  v <- 1 / (1 + i)
  k <- seq_along(qj)
  sum(benefit * v^k * ptau * qj)
}

#' Continuous multiple-decrement insurance APV \eqn{\overline{A}_{x}^{(j)}}
#'
#' Computes the actuarial present value of a benefit payable at the moment of
#' decrement by Cause \eqn{j}, matching Equation (14.4) in Chapter 14.
#'
#' The integral is evaluated numerically by the trapezoidal rule:
#' \deqn{
#' \overline{A}_{x}^{(j)} = \int_0^T v^t {}_{t}p_{x}^{(\tau)} \mu_{x+t}^{(j)} dt
#' }
#'
#' @param t Numeric vector of time points.
#' @param ptau Numeric vector of values \eqn{{}_{t}p_{x}^{(\tau)}}.
#' @param muj Numeric vector of values \eqn{\mu_{x+t}^{(j)}}.
#' @param delta Force of interest.
#' @param benefit Benefit amount payable on decrement by Cause \eqn{j}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' t <- seq(0, 20, by = 0.01)
#' ptau <- exp(-0.012 * t)
#' mu_ac <- rep(0.002, length(t))
#' Abarxj_md(t, ptau, mu_ac, delta = 0.05, benefit = 2000)
#'
#' @export
Abarxj_md <- function(t, ptau, muj, delta, benefit = 1) {
  if (!(is.numeric(t) && is.numeric(ptau) && is.numeric(muj))) {
    stop("t, ptau, and muj must be numeric vectors.", call. = FALSE)
  }
  if (!(length(t) == length(ptau) && length(t) == length(muj))) {
    stop("t, ptau, and muj must have the same length.", call. = FALSE)
  }
  if (length(t) < 2) {
    stop("Need at least two time points.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta)) {
    stop("delta must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(benefit) || length(benefit) != 1 || !is.finite(benefit)) {
    stop("benefit must be a finite numeric scalar.", call. = FALSE)
  }

  vt <- exp(-delta * t)
  y <- benefit * vt * ptau * muj
  sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
}

#' Projected asset share path \eqn{{}_{k}AS}
#'
#' Computes projected asset shares recursively using Equations (14.5b) and
#' (14.6b) of Chapter 14, with optional support for a survival benefit payable
#' at the end of year \eqn{k}.
#'
#' For policy year \eqn{k},
#' \deqn{
#' [{}_{k-1}AS + G(1-r_k) - e_k](1+i)
#' = b_k^{(1)} q_{x+k-1}^{(1)}
#' + b_k^{(2)} q_{x+k-1}^{(2)}
#' + p_{x+k-1}^{(\tau)} \left(b_k^{(3)} + {}_{k}AS\right)
#' }
#' so that
#' \deqn{
#' {}_{k}AS =
#' \frac{[{}_{k-1}AS + G(1-r_k) - e_k](1+i)
#' - b_k^{(1)} q_{x+k-1}^{(1)}
#' - b_k^{(2)} q_{x+k-1}^{(2)}}{p_{x+k-1}^{(\tau)}} - b_k^{(3)}
#' }
#'
#' @param AS0 Initial asset share \eqn{{}_{0}AS}.
#' @param G Level annual premium.
#' @param r Numeric vector of percent-of-premium expense factors.
#' @param e Numeric vector of fixed contract expenses.
#' @param b1 Numeric vector of Cause 1 benefit amounts.
#' @param b2 Numeric vector of Cause 2 benefit amounts.
#' @param q1 Numeric vector of Cause 1 decrement probabilities.
#' @param q2 Numeric vector of Cause 2 decrement probabilities.
#' @param p_tau Numeric vector of in-force probabilities.
#' @param i Effective annual interest rate.
#' @param b3 Optional numeric vector of survival benefit amounts payable at the
#'   end of year \eqn{k} conditional on survival through year \eqn{k}. Defaults
#'   to a zero vector.
#'
#' @return A data frame with columns \code{k} and \code{AS}.
#'
#' @export
AS_path <- function(AS0, G, r, e, b1, b2, q1, q2, p_tau, i, b3 = NULL) {
  n <- length(r)
  lens <- c(length(e), length(b1), length(b2), length(q1), length(q2), length(p_tau))
  if (any(lens != n)) {
    stop("All yearly vectors must have the same length.", call. = FALSE)
  }

  if (is.null(b3)) {
    b3 <- rep(0, n)
  }
  if (length(b3) != n) {
    stop("b3 must have the same length as the other yearly vectors.", call. = FALSE)
  }

  AS <- numeric(n + 1)
  AS[1] <- AS0

  for (k in seq_len(n)) {
    AS[k + 1] <- ((AS[k] + G * (1 - r[k]) - e[k]) * (1 + i) -
                    b1[k] * q1[k] - b2[k] * q2[k]) / p_tau[k] - b3[k]
  }

  data.frame(
    k = 0:n,
    AS = AS
  )
}

#' Euler approximation for \eqn{{}_{t}p_{x}^{00}} and \eqn{{}_{t}p_{x}^{01}}
#'
#' Computes the Euler approximations in the disability model allowing for
#' recovery, as in Equations (14.20) and (14.21).
#'
#' The model uses three states:
#' \itemize{
#'   \item State 0: healthy
#'   \item State 1: disabled
#'   \item State 2: deceased
#' }
#'
#' @param h Step size.
#' @param n Final time.
#' @param mu01 Function of time returning \eqn{\mu_{x+t}^{01}}.
#' @param mu02 Function of time returning \eqn{\mu_{x+t}^{02}}.
#' @param mu10 Function of time returning \eqn{\mu_{x+t}^{10}}.
#' @param mu12 Function of time returning \eqn{\mu_{x+t}^{12}}.
#' @param p00_0 Initial value of \eqn{{}_{0}p_{x}^{00}}.
#' @param p01_0 Initial value of \eqn{{}_{0}p_{x}^{01}}.
#'
#' @return A data frame with columns \code{t}, \code{tp00}, \code{tp01},
#'   and \code{tp02}.
#'
#' @examples
#' mu01 <- function(t) 0.10 * t + 0.20
#' mu02 <- function(t) 0.20
#' mu10 <- function(t) 0.50
#' mu12 <- function(t) 0.125 * t + 0.20
#'
#' tp00_tp01_euler(
#'   h = 0.10, n = 2.0,
#'   mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
#' )
#'
#' @export
tp00_tp01_euler <- function(h, n, mu01, mu02, mu10, mu12,
                            p00_0 = 1, p01_0 = 0) {
  t_grid <- seq(0, n, by = h)
  m <- length(t_grid)

  p00 <- numeric(m)
  p01 <- numeric(m)
  p02 <- numeric(m)

  p00[1] <- p00_0
  p01[1] <- p01_0
  p02[1] <- 1 - p00_0 - p01_0

  for (idx in 1:(m - 1)) {
    t <- t_grid[idx]

    dp00 <- p01[idx] * mu10(t) - p00[idx] * (mu01(t) + mu02(t))
    dp01 <- p00[idx] * mu01(t) - p01[idx] * (mu10(t) + mu12(t))

    p00[idx + 1] <- p00[idx] + h * dp00
    p01[idx + 1] <- p01[idx] + h * dp01
    p02[idx + 1] <- 1 - p00[idx + 1] - p01[idx + 1]
  }

  data.frame(
    t = t_grid,
    tp00 = p00,
    tp01 = p01,
    tp02 = p02
  )
}

#' Continuous premium approximation \eqn{\overline{P}} by trapezoidal rule
#'
#' Approximates the annual continuous premium in the disability model allowing
#' for recovery, as in Example 14.18.
#'
#' The numerator is
#' \deqn{
#' \int v^t \left[{}_{t}p_{x}^{00}\mu_{x+t}^{02}B^{02}
#' + {}_{t}p_{x}^{01}\mu_{x+t}^{12}B^{12}
#' + {}_{t}p_{x}^{01}R \right] dt
#' }
#' and the denominator is
#' \deqn{
#' \int v^t {}_{t}p_{x}^{00} dt
#' }
#'
#' @param t Numeric vector of time points.
#' @param tp00 Numeric vector of values \eqn{{}_{t}p_{x}^{00}}.
#' @param tp01 Numeric vector of values \eqn{{}_{t}p_{x}^{01}}.
#' @param delta Force of interest.
#' @param mu02 Function of time returning \eqn{\mu_{x+t}^{02}}.
#' @param mu12 Function of time returning \eqn{\mu_{x+t}^{12}}.
#' @param B02 Benefit payable on death while healthy.
#' @param B12 Benefit payable on death while disabled.
#' @param R Continuous income rate while disabled.
#'
#' @return A numeric scalar.
#'
#' @examples
#' mu01 <- function(t) 0.10 * t + 0.20
#' mu02 <- function(t) 0.20
#' mu10 <- function(t) 0.50
#' mu12 <- function(t) 0.125 * t + 0.20
#'
#' ex1410 <- tp00_tp01_euler(
#'   h = 0.10, n = 2.0,
#'   mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
#' )
#'
#' Pbar_trapz_ms(
#'   t = ex1410$t,
#'   tp00 = ex1410$tp00,
#'   tp01 = ex1410$tp01,
#'   delta = 0.04,
#'   mu02 = mu02,
#'   mu12 = mu12,
#'   B02 = 1000,
#'   B12 = 1000,
#'   R = 1000
#' )
#'
#' @export
Pbar_trapz_ms <- function(t, tp00, tp01, delta,
                          mu02, mu12,
                          B02 = 1, B12 = 1, R = 0) {
  if (!(length(t) == length(tp00) && length(t) == length(tp01))) {
    stop("t, tp00, and tp01 must have the same length.", call. = FALSE)
  }

  vt <- exp(-delta * t)

  num_y <- vt * (tp00 * mu02(t) * B02 + tp01 * mu12(t) * B12 + tp01 * R)
  den_y <- vt * tp00

  num <- sum(diff(t) * (head(num_y, -1) + tail(num_y, -1)) / 2)
  den <- sum(diff(t) * (head(den_y, -1) + tail(den_y, -1)) / 2)

  num / den
}

#' Reserve derivatives for the disability model with recovery
#'
#' Computes the right-hand sides of the coupled Thiele differential equations
#' in Equations (14.25) and (14.26) for the healthy-life reserve
#' \eqn{{}_{t}\overline{V}^{(0)}} and the disabled-life reserve
#' \eqn{{}_{t}\overline{V}^{(1)}}.
#'
#' The equations are
#' \deqn{
#' \frac{d}{dt}{}_{t}\overline{V}^{(0)}
#' = \overline{P} + \delta {}_{t}\overline{V}^{(0)}
#' - \mu_{x+t}^{02}(B-{}_{t}\overline{V}^{(0)})
#' - \mu_{x+t}^{01}({}_{t}\overline{V}^{(1)}-{}_{t}\overline{V}^{(0)})
#' }
#' and
#' \deqn{
#' \frac{d}{dt}{}_{t}\overline{V}^{(1)}
#' = \delta {}_{t}\overline{V}^{(1)} - R
#' - \mu_{x+t}^{12}(B-{}_{t}\overline{V}^{(1)})
#' - \mu_{x+t}^{10}({}_{t}\overline{V}^{(0)}-{}_{t}\overline{V}^{(1)})
#' }
#'
#' @param t Time.
#' @param V0 Value of \eqn{{}_{t}\overline{V}^{(0)}}.
#' @param V1 Value of \eqn{{}_{t}\overline{V}^{(1)}}.
#' @param delta Force of interest.
#' @param Pbar Continuous premium rate.
#' @param B Death benefit.
#' @param R Continuous disability income rate.
#' @param mu01 Function of time returning \eqn{\mu_{x+t}^{01}}.
#' @param mu02 Function of time returning \eqn{\mu_{x+t}^{02}}.
#' @param mu10 Function of time returning \eqn{\mu_{x+t}^{10}}.
#' @param mu12 Function of time returning \eqn{\mu_{x+t}^{12}}.
#'
#' @return A named numeric vector with components \code{dV0} and \code{dV1}.
#'
#' @examples
#' mu01 <- function(t) 0.10 * t + 0.20
#' mu02 <- function(t) 0.20
#' mu10 <- function(t) 0.50
#' mu12 <- function(t) 0.125 * t + 0.20
#'
#' thiele_dVdt_01(
#'   t = 2.0, V0 = 0, V1 = 0,
#'   delta = 0.04, Pbar = 446.95,
#'   B = 1000, R = 1000,
#'   mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
#' )
#'
#' @export
thiele_dVdt_01 <- function(t, V0, V1, delta, Pbar, B, R,
                           mu01, mu02, mu10, mu12) {
  if (!is.numeric(t) || length(t) != 1 || !is.finite(t)) {
    stop("t must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(V0) || length(V0) != 1 || !is.finite(V0)) {
    stop("V0 must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(V1) || length(V1) != 1 || !is.finite(V1)) {
    stop("V1 must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta)) {
    stop("delta must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(Pbar) || length(Pbar) != 1 || !is.finite(Pbar)) {
    stop("Pbar must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1 || !is.finite(B)) {
    stop("B must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1 || !is.finite(R)) {
    stop("R must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.function(mu01) || !is.function(mu02) ||
      !is.function(mu10) || !is.function(mu12)) {
    stop("mu01, mu02, mu10, and mu12 must be functions.", call. = FALSE)
  }

  dV0 <- Pbar + delta * V0 -
    mu02(t) * (B - V0) -
    mu01(t) * (V1 - V0)

  dV1 <- delta * V1 - R -
    mu12(t) * (B - V1) -
    mu10(t) * (V0 - V1)

  c(dV0 = dV0, dV1 = dV1)
}

#' Backward reserve path for the disability model with recovery
#'
#' Computes the backward Euler reserve path for the healthy-life reserve
#' \eqn{{}_{t}\overline{V}^{(0)}} and disabled-life reserve
#' \eqn{{}_{t}\overline{V}^{(1)}} using Equations (14.27) and (14.28).
#'
#' @param h Step size.
#' @param n Final time.
#' @param delta Force of interest.
#' @param Pbar Continuous premium rate.
#' @param B Death benefit.
#' @param R Continuous disability income rate.
#' @param mu01 Function of time returning \eqn{\mu_{x+t}^{01}}.
#' @param mu02 Function of time returning \eqn{\mu_{x+t}^{02}}.
#' @param mu10 Function of time returning \eqn{\mu_{x+t}^{10}}.
#' @param mu12 Function of time returning \eqn{\mu_{x+t}^{12}}.
#' @param V0_n Terminal value of \eqn{{}_{n}\overline{V}^{(0)}}.
#' @param V1_n Terminal value of \eqn{{}_{n}\overline{V}^{(1)}}.
#'
#' @return A data frame with columns \code{t}, \code{tV0}, and \code{tV1}.
#'
#' @examples
#' mu01 <- function(t) 0.10 * t + 0.20
#' mu02 <- function(t) 0.20
#' mu10 <- function(t) 0.50
#' mu12 <- function(t) 0.125 * t + 0.20
#'
#' thiele_path_01(
#'   h = 0.10, n = 2.0, delta = 0.04, Pbar = 446.95,
#'   B = 1000, R = 1000,
#'   mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12
#' )
#'
#' @export
thiele_path_01 <- function(h, n, delta, Pbar, B, R,
                           mu01, mu02, mu10, mu12,
                           V0_n = 0, V1_n = 0) {
  if (!is.numeric(h) || length(h) != 1 || !is.finite(h) || h <= 0) {
    stop("h must be a positive finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n < 0) {
    stop("n must be a nonnegative finite numeric scalar.", call. = FALSE)
  }

  t_grid <- seq(n, 0, by = -h)
  if (tail(t_grid, 1) != 0) {
    t_grid <- c(t_grid, 0)
  }
  m <- length(t_grid)

  V0 <- numeric(m)
  V1 <- numeric(m)
  V0[1] <- V0_n
  V1[1] <- V1_n

  for (idx in 1:(m - 1)) {
    t <- t_grid[idx]
    h_step <- t_grid[idx] - t_grid[idx + 1]

    dV <- thiele_dVdt_01(
      t = t,
      V0 = V0[idx],
      V1 = V1[idx],
      delta = delta,
      Pbar = Pbar,
      B = B,
      R = R,
      mu01 = mu01,
      mu02 = mu02,
      mu10 = mu10,
      mu12 = mu12
    )

    V0[idx + 1] <- V0[idx] - h_step * dV["dV0"]
    V1[idx + 1] <- V1[idx] - h_step * dV["dV1"]
  }

  out <- data.frame(
    t = rev(t_grid),
    tV0 = rev(V0),
    tV1 = rev(V1)
  )
  rownames(out) <- NULL
  out
}

#' n-step transition probability for a discrete-time Markov chain
#'
#' Computes the \eqn{(i,j)} entry of \eqn{P^n}, useful for Chapter 14 examples
#' involving discrete-time multi-state models such as CCRC and risk-class models.
#'
#' @param P Transition probability matrix.
#' @param n Nonnegative integer number of steps.
#' @param i Starting state index.
#' @param j Ending state index.
#'
#' @return A numeric scalar.
#'
#' @examples
#' P <- matrix(
#'   c(0.94, 0.03, 0.02, 0.01,
#'     0.50, 0.30, 0.18, 0.02,
#'     0.00, 0.00, 0.93, 0.07,
#'     0.00, 0.00, 0.00, 1.00),
#'   nrow = 4, byrow = TRUE
#' )
#'
#' markov_nstep_prob(P, n = 3, i = 1, j = 1)
#' markov_nstep_prob(P, n = 3, i = 1, j = 3)
#'
#' @export
markov_nstep_prob <- function(P, n, i, j) {
  if (!is.matrix(P)) stop("P must be a matrix.", call. = FALSE)
  if (n < 0 || n != as.integer(n)) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }
  if (i < 1 || i > nrow(P) || j < 1 || j > ncol(P)) {
    stop("i and j must be valid state indices.", call. = FALSE)
  }

  Pn <- diag(nrow(P))
  if (n > 0) {
    for (k in seq_len(n)) {
      Pn <- Pn %*% P
    }
  }

  Pn[i, j]
}

#' Gain or loss in a multiple-decrement model
#'
#' Computes the gain or loss expression from Section 14.6.
#'
#' With within-year decrement probabilities, the function evaluates
#' \deqn{
#' [{}_{t}V^G + G(1-r) - e](1+i)
#' - \left[(b^{(1)}+s^{(1)})q^{(1)} + (b^{(2)}+s^{(2)})q^{(2)}
#' + p^{(\tau)} {}_{t+1}V^G \right]
#' }
#'
#' If \code{year_end_cause2 = TRUE}, the Cause 2 decrement is treated as
#' occurring only at year end, matching Equation (14.30).
#'
#' @param Vt Gross premium reserve at time \eqn{t}.
#' @param G Gross premium for the year.
#' @param r Percent-of-premium expense factor.
#' @param e Fixed expense at the beginning of the year.
#' @param i Earned interest rate.
#' @param b1 Cause 1 benefit.
#' @param b2 Cause 2 benefit.
#' @param s1 Claim settlement expense for Cause 1.
#' @param s2 Claim settlement expense for Cause 2.
#' @param q1 Cause 1 decrement probability.
#' @param q2 Cause 2 decrement probability.
#' @param Vt1 Gross premium reserve at time \eqn{t+1}.
#' @param year_end_cause2 Logical; if \code{TRUE}, use the year-end Cause 2 form.
#' @param q1prime Single-decrement Cause 1 probability for the year-end Cause 2 case.
#' @param q2prime Single-decrement Cause 2 probability for the year-end Cause 2 case.
#'
#' @return A numeric scalar.
#'
#' @examples
#' gain_loss_md(
#'   Vt = 115.00, G = 16, r = 0, e = 3, i = 0.06,
#'   b1 = 1000, b2 = 110, s1 = 0, s2 = 0,
#'   q1 = 0.01, q2 = 0.10, Vt1 = 128.83
#' )
#'
#' @export
gain_loss_md <- function(Vt, G, r, e, i,
                         b1, b2, s1 = 0, s2 = 0,
                         q1, q2, Vt1,
                         year_end_cause2 = FALSE,
                         q1prime = NULL, q2prime = NULL) {
  lhs <- (Vt + G * (1 - r) - e) * (1 + i)

  rhs <- if (!year_end_cause2) {
    (b1 + s1) * q1 + (b2 + s2) * q2 + (1 - q1 - q2) * Vt1
  } else {
    if (is.null(q1prime) || is.null(q2prime)) {
      stop("q1prime and q2prime must be provided when year_end_cause2 = TRUE.",
           call. = FALSE)
    }
    (b1 + s1) * q1prime +
      (b2 + s2) * (1 - q1prime) * q2prime +
      (1 - q1prime) * (1 - q2prime) * Vt1
  }

  lhs - rhs
}

#' General projected asset share path (multiple decrements)
#'
#' @param AS0 Initial asset share.
#' @param G Premium.
#' @param r Expense percentages.
#' @param e Fixed expenses.
#' @param b_mat Matrix of benefits (rows = years, cols = causes).
#' @param q_mat Matrix of decrement probabilities (same shape as b_mat).
#' @param p_tau In-force probabilities.
#' @param i Interest rate.
#' @param b_surv Optional survival benefits.
#'
#' @export
AS_path_md <- function(AS0, G, r, e, b_mat, q_mat, p_tau, i, b_surv = NULL) {

  if (!is.matrix(b_mat) || !is.matrix(q_mat)) {
    stop("b_mat and q_mat must be matrices.")
  }
  if (!all(dim(b_mat) == dim(q_mat))) {
    stop("b_mat and q_mat must have same dimensions.")
  }

  n <- nrow(b_mat)

  if (length(r) != n || length(e) != n || length(p_tau) != n) {
    stop("All yearly vectors must have length n.")
  }

  if (is.null(b_surv)) {
    b_surv <- rep(0, n)
  }

  AS <- numeric(n + 1)
  AS[1] <- AS0

  for (k in seq_len(n)) {

    claim_term <- sum(b_mat[k, ] * q_mat[k, ])

    AS[k + 1] <- (
      (AS[k] + G * (1 - r[k]) - e[k]) * (1 + i)
      - claim_term
    ) / p_tau[k] - b_surv[k]
  }

  data.frame(k = 0:n, AS = AS)
}
