# =========================================================
# Chapter 17 helpers for mqriskR
# Profit Analysis
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.validate_numeric_ch17 <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite numeric vector.", name), call. = FALSE)
  }
  invisible(x)
}

.validate_nonneg_numeric_ch17 <- function(x, name) {
  .validate_numeric_ch17(x, name)
  if (any(x < 0)) {
    stop(sprintf("%s must be nonnegative.", name), call. = FALSE)
  }
  invisible(x)
}

.validate_prob_vector_ch17 <- function(x, name) {
  .validate_numeric_ch17(x, name)
  if (any(x < 0 | x > 1)) {
    stop(sprintf("%s must contain values in [0, 1].", name), call. = FALSE)
  }
  invisible(x)
}

.recycle_ch17 <- function(x, n, name) {
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(sprintf("%s must have length 1 or %d.", name, n), call. = FALSE)
  }
  x
}

# -------------------------------------------------------------------------
# Profit vector
# -------------------------------------------------------------------------

#' Profit vector for a discrete profit-analysis model
#'
#' Computes the Chapter 17 profit vector
#' \deqn{
#' \mathbf{Pr} = (Pr_0, Pr_1, \dots, Pr_n)
#' }
#' where \eqn{Pr_0} is the negative pre-contract expense and the yearly
#' expected profit values are calculated from the general discrete expression
#' in Equation (17.1).
#'
#' This implementation allows for two decrements, typically death and
#' withdrawal/surrender.
#'
#' @param V Vector of gross premium reserves \eqn{{}_tV^G} of length
#'   \eqn{n+1}, including the issue-time reserve and the terminal reserve.
#' @param G Gross premium vector for policy years 1 through \eqn{n}.
#' @param i Interest-rate vector for policy years 1 through \eqn{n}.
#' @param r Percent-of-premium expense vector.
#' @param e Fixed expense vector.
#' @param q1 First decrement probabilities, typically death.
#' @param q2 Second decrement probabilities, typically surrender or lapse.
#'   Defaults to 0.
#' @param b1 Benefit vector for decrement 1.
#' @param b2 Benefit vector for decrement 2. Defaults to 0.
#' @param s1 Settlement-expense vector for decrement 1. Defaults to 0.
#' @param s2 Settlement-expense vector for decrement 2. Defaults to 0.
#' @param p_tau Optional vector of in-force probabilities
#'   \eqn{p_{x+t}^{(\tau)}}. If omitted, it is computed as
#'   \eqn{1-q^{(1)}-q^{(2)}}.
#' @param pre_contract_expense Positive pre-contract expense amount. The
#'   returned first element is \eqn{Pr_0 = -\text{pre\_contract\_expense}}.
#'
#' @return Numeric vector of length \eqn{n+1}.
#'
#' @examples
#' V <- c(0, 5.66, 6.17, 0)
#' qx <- c(0.00142, 0.00153, 0.00166)
#' Pr_vector_disc(
#'   V = V, G = 95, i = 0.06, r = 0.05, e = 10,
#'   q1 = qx, b1 = 50000, pre_contract_expense = 15
#' )
#'
#' @export
Pr_vector_disc <- function(V, G, i, r = 0, e = 0,
                           q1, q2 = 0,
                           b1, b2 = 0,
                           s1 = 0, s2 = 0,
                           p_tau = NULL,
                           pre_contract_expense = 0) {
  .validate_numeric_ch17(V, "V")
  .validate_nonneg_numeric_ch17(G, "G")
  .validate_nonneg_numeric_ch17(i, "i")
  .validate_nonneg_numeric_ch17(r, "r")
  .validate_nonneg_numeric_ch17(e, "e")
  .validate_prob_vector_ch17(q1, "q1")
  .validate_prob_vector_ch17(q2, "q2")
  .validate_nonneg_numeric_ch17(b1, "b1")
  .validate_nonneg_numeric_ch17(b2, "b2")
  .validate_nonneg_numeric_ch17(s1, "s1")
  .validate_nonneg_numeric_ch17(s2, "s2")
  .validate_nonneg_numeric_ch17(pre_contract_expense, "pre_contract_expense")

  if (length(V) < 2L) {
    stop("V must have length at least 2.", call. = FALSE)
  }

  n <- length(V) - 1L

  G  <- .recycle_ch17(G,  n, "G")
  i  <- .recycle_ch17(i,  n, "i")
  r  <- .recycle_ch17(r,  n, "r")
  e  <- .recycle_ch17(e,  n, "e")
  q1 <- .recycle_ch17(q1, n, "q1")
  q2 <- .recycle_ch17(q2, n, "q2")
  b1 <- .recycle_ch17(b1, n, "b1")
  b2 <- .recycle_ch17(b2, n, "b2")
  s1 <- .recycle_ch17(s1, n, "s1")
  s2 <- .recycle_ch17(s2, n, "s2")

  if (is.null(p_tau)) {
    p_tau <- 1 - q1 - q2
  } else {
    .validate_prob_vector_ch17(p_tau, "p_tau")
    p_tau <- .recycle_ch17(p_tau, n, "p_tau")
  }

  if (any(p_tau < 0)) {
    stop("Some values of p_tau are negative.", call. = FALSE)
  }

  out <- numeric(n + 1L)
  out[1L] <- -pre_contract_expense

  for (k in seq_len(n)) {
    out[k + 1L] <-
      (V[k] + G[k] * (1 - r[k]) - e[k]) * (1 + i[k]) -
      ((b1[k] + s1[k]) * q1[k] +
         (b2[k] + s2[k]) * q2[k] +
         V[k + 1L] * p_tau[k])
  }

  names(out) <- paste0("Pr", 0:n)
  out
}

# -------------------------------------------------------------------------
# Profit signature
# -------------------------------------------------------------------------

#' Profit signature from a profit vector
#'
#' Converts a Chapter 17 profit vector
#' \eqn{\mathbf{Pr}=(Pr_0,\dots,Pr_n)}
#' into the corresponding profit signature
#' \eqn{\mathbf{\Pi}=(\Pi_0,\dots,\Pi_n)}
#' using Equation (17.3).
#'
#' @param Pr Profit vector of length \eqn{n+1}.
#' @param p_tau One-year in-force probabilities. This may have length
#'   \eqn{n-1} or \eqn{n}. If length \eqn{n}, the final entry is ignored for
#'   the profit-signature calculation.
#'
#' @return Numeric vector of length \eqn{n+1}.
#'
#' @examples
#' Pr <- c(-15.00, 8.42, 8.40, 8.61)
#' Pi_signature(Pr, p_tau = c(0.99858, 0.99847, 0.99834))
#'
#' @export
Pi_signature <- function(Pr, p_tau) {
  .validate_numeric_ch17(Pr, "Pr")
  .validate_prob_vector_ch17(p_tau, "p_tau")

  if (length(Pr) < 2L) {
    stop("Pr must have length at least 2.", call. = FALSE)
  }

  n <- length(Pr) - 1L

  if (!(length(p_tau) %in% c(n - 1L, n))) {
    stop("p_tau must have length n - 1 or n, where length(Pr) = n + 1.", call. = FALSE)
  }

  p_use <- p_tau[seq_len(n - 1L)]
  start_probs <- c(1, cumprod(p_use))

  out <- c(Pr[1L], Pr[-1L] * start_probs)
  names(out) <- paste0("Pi", 0:n)
  out
}

# -------------------------------------------------------------------------
# NPV, IRR, payback, profit margin
# -------------------------------------------------------------------------

#' Net present value of a profit signature
#'
#' Computes the Chapter 17 net present value:
#' \deqn{
#' NPV = \sum_{t=0}^{n}\frac{\Pi_t}{(1+r)^t}.
#' }
#'
#' @param Pi Profit signature vector.
#' @param r Risk discount rate.
#'
#' @return Numeric scalar.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' NPV_profit(Pi, r = 0.10)
#'
#' @export
NPV_profit <- function(Pi, r) {
  .validate_numeric_ch17(Pi, "Pi")
  .validate_nonneg_numeric_ch17(r, "r")

  if (length(r) != 1L) {
    stop("r must be a numeric scalar.", call. = FALSE)
  }

  sum(Pi / (1 + r)^(seq_along(Pi) - 1L))
}

#' Partial net present values
#'
#' Computes the sequence
#' \eqn{NPV(0), NPV(1), \dots, NPV(n)}
#' of partial net present values from a profit signature.
#'
#' @param Pi Profit signature vector.
#' @param r Risk discount rate.
#'
#' @return Numeric vector.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' NPV_partial(Pi, r = 0.10)
#'
#' @export
NPV_partial <- function(Pi, r) {
  .validate_numeric_ch17(Pi, "Pi")
  .validate_nonneg_numeric_ch17(r, "r")

  if (length(r) != 1L) {
    stop("r must be a numeric scalar.", call. = FALSE)
  }

  vals <- Pi / (1 + r)^(seq_along(Pi) - 1L)
  out <- cumsum(vals)
  names(out) <- paste0("NPV(", 0:(length(Pi) - 1L), ")")
  out
}

#' Discounted payback period
#'
#' Returns the first duration \eqn{t} for which the partial net present value
#' \eqn{NPV(t)} is nonnegative.
#'
#' @param Pi Profit signature vector.
#' @param r Risk discount rate.
#'
#' @return Integer scalar, or \code{NA_integer_} if the payback period is not
#'   reached.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' discounted_payback_period(Pi, r = 0.10)
#'
#' @export
discounted_payback_period <- function(Pi, r) {
  partial <- NPV_partial(Pi = Pi, r = r)
  idx <- which(partial >= 0)

  if (length(idx) == 0L) {
    return(NA_integer_)
  }

  as.integer(idx[1L] - 1L)
}

#' Internal rate of return of a profit signature
#'
#' Computes the internal rate of return (IRR), defined as the rate \eqn{r} for
#' which the net present value is zero.
#'
#' @param Pi Profit signature vector.
#' @param interval Numeric vector of length 2 giving the search interval for
#'   \code{uniroot()}.
#' @param tol Tolerance passed to \code{uniroot()}.
#'
#' @return Numeric scalar.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' IRR_profit(Pi)
#'
#' @export
IRR_profit <- function(Pi, interval = c(0, 1), tol = .Machine$double.eps^0.5) {
  .validate_numeric_ch17(Pi, "Pi")
  .validate_numeric_ch17(interval, "interval")
  .validate_nonneg_numeric_ch17(tol, "tol")

  if (length(interval) != 2L || interval[1L] <= -1 || interval[1L] >= interval[2L]) {
    stop("interval must be length 2 with -1 < interval[1] < interval[2].", call. = FALSE)
  }

  f <- function(x) NPV_profit(Pi = Pi, r = x)

  if (f(interval[1L]) * f(interval[2L]) > 0) {
    stop("NPV does not change sign on the supplied interval.", call. = FALSE)
  }

  uniroot(f, interval = interval, tol = tol)$root
}

#' APV of gross premiums under a risk discount rate
#'
#' Computes the actuarial present value of gross premiums:
#' \deqn{
#' APV_{GP} = \sum_{t=0}^{n-1} \frac{G_{t+1} \cdot {}_tp_x^{(\tau)}}{(1+r)^t}.
#' }
#'
#' @param G Gross premium vector for policy years 1 through \eqn{n}.
#' @param r Risk discount rate.
#' @param p_tau One-year in-force probabilities. This may have length
#'   \eqn{n-1} or \eqn{n}. If length \eqn{n}, the final entry is ignored.
#'
#' @return Numeric scalar.
#'
#' @examples
#' APV_gross_premiums(G = rep(95, 3),r = 0.10,p_tau = c(0.99858, 0.99847, 0.99834))
#'
#' @export
APV_gross_premiums <- function(G, r, p_tau) {
  .validate_nonneg_numeric_ch17(G, "G")
  .validate_nonneg_numeric_ch17(r, "r")
  .validate_prob_vector_ch17(p_tau, "p_tau")

  if (length(r) != 1L) {
    stop("r must be a numeric scalar.", call. = FALSE)
  }

  n <- length(G)

  if (!(length(p_tau) %in% c(n - 1L, n))) {
    stop("p_tau must have length n - 1 or n, where n = length(G).", call. = FALSE)
  }

  p_use <- if (n == 1L) numeric(0) else p_tau[seq_len(n - 1L)]
  start_probs <- c(1, cumprod(p_use))

  sum(G * start_probs / (1 + r)^(0:(n - 1L)))
}

#' Profit margin
#'
#' Computes the Chapter 17 profit margin:
#' \deqn{
#' \text{Profit Margin} = \frac{NPV}{APV_{GP}}.
#' }
#'
#' @param NPV Net present value of profits.
#' @param APV_GP Actuarial present value of gross premiums.
#'
#' @return Numeric scalar.
#'
#' @examples
#' profit_margin(6.03, 259.52)
#'
#' @export
profit_margin <- function(NPV, APV_GP) {
  .validate_numeric_ch17(NPV, "NPV")
  .validate_nonneg_numeric_ch17(APV_GP, "APV_GP")

  if (length(NPV) != 1L || length(APV_GP) != 1L) {
    stop("NPV and APV_GP must be numeric scalars.", call. = FALSE)
  }
  if (APV_GP == 0) {
    stop("APV_GP must be positive.", call. = FALSE)
  }

  NPV / APV_GP
}

# -------------------------------------------------------------------------
# Zeroized reserves
# -------------------------------------------------------------------------

#' Zeroized reserves for a discrete death-only contract
#'
#' Computes the zeroized reserve sequence by backward recursion, setting
#' negative reserves equal to zero.
#'
#' For a death-only contract with no settlement expense and no second
#' decrement, the recursion sets
#' \deqn{
#' Pr_{t+1} = ({}_tV^Z + G_{t+1}(1-r_{t+1}) - e_{t+1})(1+i_{t+1})
#' - [Bq_{x+t} + {}_{t+1}V^Z p_{x+t}]
#' }
#' equal to zero, solving backward for \eqn{{}_tV^Z}.
#'
#' @param qx Mortality vector.
#' @param i Interest-rate vector.
#' @param G Gross premium vector.
#' @param benefit Death-benefit vector.
#' @param r Percent-of-premium expense vector.
#' @param e Fixed-expense vector.
#' @param V_terminal Terminal reserve. Defaults to 0.
#' @param floor_zero Logical; if \code{TRUE}, negative reserves are reset to 0.
#'
#' @return Numeric vector of zeroized reserves of length \eqn{n+1}.
#'
#' @examples
#' V_zeroized(
#'   qx = c(.015, .017, .019, .021, .024),
#'   i = 0.06,
#'   G = 19279,
#'   benefit = 1000000,
#'   e = 240
#' )
#'
#' @export
V_zeroized <- function(qx, i, G, benefit, r = 0, e = 0,
                       V_terminal = 0, floor_zero = TRUE) {
  .validate_prob_vector_ch17(qx, "qx")
  .validate_nonneg_numeric_ch17(i, "i")
  .validate_nonneg_numeric_ch17(G, "G")
  .validate_nonneg_numeric_ch17(benefit, "benefit")
  .validate_nonneg_numeric_ch17(r, "r")
  .validate_nonneg_numeric_ch17(e, "e")
  .validate_nonneg_numeric_ch17(V_terminal, "V_terminal")

  n <- length(qx)

  i <- .recycle_ch17(i, n, "i")
  G <- .recycle_ch17(G, n, "G")
  benefit <- .recycle_ch17(benefit, n, "benefit")
  r <- .recycle_ch17(r, n, "r")
  e <- .recycle_ch17(e, n, "e")

  V <- numeric(n + 1L)
  V[n + 1L] <- V_terminal

  for (k in n:1L) {
    p <- 1 - qx[k]
    net_prem <- G[k] * (1 - r[k]) - e[k]

    Vk <- (benefit[k] * qx[k] + V[k + 1L] * p) / (1 + i[k]) - net_prem

    if (floor_zero) {
      Vk <- max(Vk, 0)
    }

    V[k] <- Vk
  }

  names(V) <- paste0("V", 0:n)
  V
}
