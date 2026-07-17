# =========================================================
# Profit analysis functions
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.profit_check_numeric <- function(x, name, allow_empty = FALSE) {
  if (!is.numeric(x)) {
    stop(name, " must be numeric.", call. = FALSE)
  }
  if (!allow_empty && length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain finite values.", call. = FALSE)
  }
  as.numeric(x)
}

#' @noRd
.profit_check_nonnegative <- function(x, name, allow_empty = FALSE) {
  x <- .profit_check_numeric(x, name, allow_empty = allow_empty)
  if (any(x < 0)) {
    stop(name, " must contain nonnegative values.", call. = FALSE)
  }
  x
}

#' @noRd
.profit_check_probability <- function(x, name, allow_empty = FALSE) {
  x <- .profit_check_numeric(x, name, allow_empty = allow_empty)
  if (any(x < 0 | x > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }
  x
}

#' @noRd
.profit_check_rate <- function(x, name) {
  x <- .profit_check_numeric(x, name)
  if (any(x <= -1)) {
    stop(name, " must contain values greater than -1.", call. = FALSE)
  }
  x
}

#' @noRd
.profit_check_scalar <- function(x, name, nonnegative = FALSE, positive = FALSE) {
  if (nonnegative || positive) {
    x <- .profit_check_nonnegative(x, name)
  } else {
    x <- .profit_check_numeric(x, name)
  }
  if (length(x) != 1L) {
    stop(name, " must be a numeric scalar.", call. = FALSE)
  }
  if (positive && x <= 0) {
    stop(name, " must be positive.", call. = FALSE)
  }
  x
}

#' @noRd
.profit_check_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

#' @noRd
.profit_recycle_to <- function(x, n, name) {
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(name, " must have length 1 or ", n, ".", call. = FALSE)
  }
  x
}

#' @noRd
.profit_recycle_common <- function(..., .names = NULL) {
  values <- list(...)
  lengths <- vapply(values, length, integer(1))
  if (any(lengths == 0L)) {
    stop("Inputs must have positive length.", call. = FALSE)
  }
  common_length <- max(lengths)
  if (is.null(.names)) {
    .names <- paste0("Argument ", seq_along(values))
  }
  if (length(.names) != length(values)) {
    stop(".names must have the same length as the supplied arguments.", call. = FALSE)
  }
  for (j in seq_along(values)) {
    if (!lengths[j] %in% c(1L, common_length)) {
      stop(
        .names[j],
        " must have length 1 or the common length ",
        common_length,
        ".",
        call. = FALSE
      )
    }
    values[[j]] <- rep(values[[j]], length.out = common_length)
  }
  values
}

#' @noRd
.profit_start_probabilities <- function(p_tau, n, name = "p_tau") {
  if (n < 1L) {
    stop("n must be positive.", call. = FALSE)
  }
  if (n == 1L) {
    if (!length(p_tau) %in% c(0L, 1L)) {
      stop(name, " must have length 0 or 1 when n = 1.", call. = FALSE)
    }
    if (length(p_tau) == 1L) {
      .profit_check_probability(p_tau, name)
    }
    return(1)
  }
  p_tau <- .profit_check_probability(p_tau, name)
  if (!length(p_tau) %in% c(n - 1L, n)) {
    stop(name, " must have length n - 1 or n.", call. = FALSE)
  }
  c(1, cumprod(p_tau[seq_len(n - 1L)]))
}

# -------------------------------------------------------------------------
# Profit vector and profit signature
# -------------------------------------------------------------------------

#' Profit vector for a discrete profit-analysis model
#'
#' Computes expected profit by policy year for a discrete contract with up to
#' two decrements. The first element is the negative pre-contract expense.
#'
#' For policy year \eqn{k}, the expected profit is
#' \deqn{
#' [V_{k-1} + G_k(1-r_k)-e_k](1+i_k)
#' -
#' [(b_k^{(1)}+s_k^{(1)})q_k^{(1)}
#' +(b_k^{(2)}+s_k^{(2)})q_k^{(2)}
#' +V_kp_k^{(\tau)}].
#' }
#'
#' Scalar yearly inputs are recycled to the number of policy years determined
#' by \code{length(V) - 1}.
#'
#' @param V Numeric vector of gross premium reserves with length
#'   \code{n + 1}, including the issue-time and terminal reserves.
#' @param G Gross premium by policy year.
#' @param i Annual effective interest rate by policy year. Values must be
#'   greater than \code{-1}.
#' @param r Percent-of-premium expense rate by policy year. Values must lie in
#'   \code{[0, 1]}.
#' @param e Fixed expense by policy year.
#' @param q1 Probability of the first decrement by policy year.
#' @param q2 Probability of the second decrement by policy year.
#' @param b1 Benefit payable on the first decrement.
#' @param b2 Benefit payable on the second decrement.
#' @param s1 Settlement expense associated with the first decrement.
#' @param s2 Settlement expense associated with the second decrement.
#' @param p_tau Optional in-force probability by policy year. If omitted, it is
#'   calculated as \code{1 - q1 - q2}.
#' @param pre_contract_expense Nonnegative scalar pre-contract expense.
#'
#' @return A named numeric vector of length \code{n + 1}.
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
  V <- .profit_check_numeric(V, "V")
  G <- .profit_check_nonnegative(G, "G")
  i <- .profit_check_rate(i, "i")
  r <- .profit_check_probability(r, "r")
  e <- .profit_check_nonnegative(e, "e")
  q1 <- .profit_check_probability(q1, "q1")
  q2 <- .profit_check_probability(q2, "q2")
  b1 <- .profit_check_nonnegative(b1, "b1")
  b2 <- .profit_check_nonnegative(b2, "b2")
  s1 <- .profit_check_nonnegative(s1, "s1")
  s2 <- .profit_check_nonnegative(s2, "s2")
  pre_contract_expense <- .profit_check_scalar(
    pre_contract_expense, "pre_contract_expense", nonnegative = TRUE
  )

  if (length(V) < 2L) {
    stop("V must have length at least 2.", call. = FALSE)
  }

  n <- length(V) - 1L
  G <- .profit_recycle_to(G, n, "G")
  i <- .profit_recycle_to(i, n, "i")
  r <- .profit_recycle_to(r, n, "r")
  e <- .profit_recycle_to(e, n, "e")
  q1 <- .profit_recycle_to(q1, n, "q1")
  q2 <- .profit_recycle_to(q2, n, "q2")
  b1 <- .profit_recycle_to(b1, n, "b1")
  b2 <- .profit_recycle_to(b2, n, "b2")
  s1 <- .profit_recycle_to(s1, n, "s1")
  s2 <- .profit_recycle_to(s2, n, "s2")

  if (is.null(p_tau)) {
    p_tau <- 1 - q1 - q2
    if (any(p_tau < 0)) {
      stop("q1 + q2 must not exceed 1.", call. = FALSE)
    }
  } else {
    p_tau <- .profit_check_probability(p_tau, "p_tau")
    p_tau <- .profit_recycle_to(p_tau, n, "p_tau")
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

#' Profit signature
#'
#' Converts policy-year expected profits into a profit signature by weighting
#' each future expected profit by the probability that the contract is in
#' force at the start of that policy year.
#'
#' @param Pr Profit vector of length \code{n + 1}.
#' @param p_tau One-year in-force probabilities. For \code{n > 1}, this may
#'   have length \code{n - 1} or \code{n}; the final value is ignored when
#'   length \code{n}. For a one-year contract, use \code{numeric(0)} or a
#'   single probability.
#'
#' @return A named numeric vector with the same length as \code{Pr}.
#'
#' @examples
#' Pr <- c(-15.00, 8.42, 8.40, 8.61)
#' Pi_signature(Pr, p_tau = c(0.99858, 0.99847, 0.99834))
#'
#' @export
Pi_signature <- function(Pr, p_tau) {
  Pr <- .profit_check_numeric(Pr, "Pr")
  if (length(Pr) < 2L) {
    stop("Pr must have length at least 2.", call. = FALSE)
  }
  n <- length(Pr) - 1L
  start_probabilities <- .profit_start_probabilities(p_tau = p_tau, n = n)
  out <- c(Pr[1L], Pr[-1L] * start_probabilities)
  names(out) <- paste0("Pi", 0:n)
  out
}

# -------------------------------------------------------------------------
# Profitability measures
# -------------------------------------------------------------------------

#' Net present value of a profit signature
#'
#' Discounts a profit signature at one or more annual effective risk discount
#' rates.
#'
#' @param Pi Numeric profit-signature vector.
#' @param r Annual effective risk discount rate. May be scalar or vector;
#'   values must be greater than \code{-1}.
#'
#' @return A numeric vector with one value for each rate in \code{r}.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' NPV_profit(Pi, r = 0.10)
#' NPV_profit(Pi, r = c(0.08, 0.10, 0.12))
#'
#' @export
NPV_profit <- function(Pi, r) {
  Pi <- .profit_check_numeric(Pi, "Pi")
  r <- .profit_check_rate(r, "r")
  times <- seq_along(Pi) - 1L

  vapply(
    r,
    function(rate) sum(Pi / (1 + rate)^times),
    numeric(1)
  )
}

#' Partial net present values
#'
#' Computes cumulative discounted profits through each duration. A scalar
#' discount rate returns a named vector; vectorized rates return a matrix.
#'
#' @inheritParams NPV_profit
#'
#' @return A named numeric vector for scalar \code{r}, or a numeric matrix for
#'   vectorized \code{r}.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' NPV_partial(Pi, r = 0.10)
#'
#' @export
NPV_partial <- function(Pi, r) {
  Pi <- .profit_check_numeric(Pi, "Pi")
  r <- .profit_check_rate(r, "r")
  times <- seq_along(Pi) - 1L

  rows <- lapply(r, function(rate) cumsum(Pi / (1 + rate)^times))
  column_names <- paste0("NPV(", times, ")")

  if (length(rows) == 1L) {
    return(stats::setNames(rows[[1L]], column_names))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  colnames(out) <- column_names
  out
}

#' Discounted payback period
#'
#' Returns the first duration at which cumulative discounted profit is
#' nonnegative.
#'
#' @inheritParams NPV_profit
#'
#' @return An integer vector with one value for each discount rate. An element
#'   is \code{NA_integer_} when payback is not reached.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' discounted_payback_period(Pi, r = 0.10)
#'
#' @export
discounted_payback_period <- function(Pi, r) {
  Pi <- .profit_check_numeric(Pi, "Pi")
  r <- .profit_check_rate(r, "r")
  partial <- NPV_partial(Pi = Pi, r = r)

  if (is.null(dim(partial))) {
    index <- which(partial >= 0)
    if (length(index) == 0L) {
      return(NA_integer_)
    }
    return(as.integer(index[1L] - 1L))
  }

  apply(partial, 1L, function(x) {
    index <- which(x >= 0)
    if (length(index) == 0L) {
      return(NA_integer_)
    }
    as.integer(index[1L] - 1L)
  })
}

#' Internal rate of return
#'
#' Computes a root of the profit-signature net present value.
#'
#' @param Pi Numeric profit-signature vector.
#' @param interval Numeric vector of length two giving the root-search
#'   interval. Its lower endpoint must be greater than \code{-1}.
#' @param tol Positive scalar tolerance passed to \code{stats::uniroot()}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' Pi <- c(-15.00, 8.42, 8.39, 8.58)
#' IRR_profit(Pi)
#'
#' @export
IRR_profit <- function(Pi, interval = c(0, 1),
                       tol = .Machine$double.eps^0.5) {
  Pi <- .profit_check_numeric(Pi, "Pi")
  interval <- .profit_check_numeric(interval, "interval")
  tol <- .profit_check_scalar(tol, "tol", positive = TRUE)

  if (length(interval) != 2L ||
      interval[1L] <= -1 ||
      interval[1L] >= interval[2L]) {
    stop(
      "interval must have length 2 and satisfy -1 < interval[1] < interval[2].",
      call. = FALSE
    )
  }

  objective <- function(rate) NPV_profit(Pi = Pi, r = rate)
  lower_value <- objective(interval[1L])
  upper_value <- objective(interval[2L])

  if (lower_value == 0) return(interval[1L])
  if (upper_value == 0) return(interval[2L])

  if (lower_value * upper_value > 0) {
    stop("NPV does not change sign over the supplied interval.", call. = FALSE)
  }

  stats::uniroot(objective, interval = interval, tol = tol)$root
}

#' Actuarial present value of gross premiums
#'
#' Computes the actuarial present value of premiums weighted by contract
#' persistency at the start of each policy year.
#'
#' @param G Nonnegative gross premium vector.
#' @param r Annual effective risk discount rate. May be scalar or vector;
#'   values must be greater than \code{-1}.
#' @param p_tau One-year in-force probabilities. For \code{n > 1}, this may
#'   have length \code{n - 1} or \code{n}; the final value is ignored when
#'   length \code{n}. For one premium, use \code{numeric(0)} or one value.
#'
#' @return A numeric vector with one value for each rate in \code{r}.
#'
#' @examples
#' APV_gross_premiums(
#'   G = rep(95, 3),
#'   r = 0.10,
#'   p_tau = c(0.99858, 0.99847, 0.99834)
#' )
#'
#' @export
APV_gross_premiums <- function(G, r, p_tau) {
  G <- .profit_check_nonnegative(G, "G")
  r <- .profit_check_rate(r, "r")
  n <- length(G)
  start_probabilities <- .profit_start_probabilities(p_tau = p_tau, n = n)
  times <- 0:(n - 1L)

  vapply(
    r,
    function(rate) sum(G * start_probabilities / (1 + rate)^times),
    numeric(1)
  )
}

#' Profit margin
#'
#' Computes net present value divided by the actuarial present value of gross
#' premiums. Scalar arguments are recycled to the common length.
#'
#' @param NPV Net present value of profits.
#' @param APV_GP Positive actuarial present value of gross premiums.
#'
#' @return A numeric vector.
#'
#' @examples
#' profit_margin(6.03, 259.52)
#'
#' @export
profit_margin <- function(NPV, APV_GP) {
  NPV <- .profit_check_numeric(NPV, "NPV")
  APV_GP <- .profit_check_nonnegative(APV_GP, "APV_GP")

  values <- .profit_recycle_common(
    NPV, APV_GP, .names = c("NPV", "APV_GP")
  )
  NPV <- values[[1]]
  APV_GP <- values[[2]]

  if (any(APV_GP <= 0)) {
    stop("APV_GP must contain positive values.", call. = FALSE)
  }

  NPV / APV_GP
}

# -------------------------------------------------------------------------
# Zeroized reserves
# -------------------------------------------------------------------------

#' Zeroized reserves for a discrete death-benefit contract
#'
#' Computes reserves backward by setting expected profit in each policy year
#' equal to zero. Negative reserves may optionally be floored at zero.
#'
#' @param qx Mortality probability by policy year.
#' @param i Annual effective interest rate by policy year. Values must be
#'   greater than \code{-1}.
#' @param G Gross premium by policy year.
#' @param benefit Death benefit by policy year.
#' @param r Percent-of-premium expense rate by policy year. Values must lie in
#'   \code{[0, 1]}.
#' @param e Fixed expense by policy year.
#' @param V_terminal Nonnegative scalar terminal reserve.
#' @param floor_zero Logical scalar. If \code{TRUE}, negative reserves are
#'   replaced by zero.
#'
#' @return A named numeric vector of length \code{length(qx) + 1}.
#'
#' @examples
#' V_zeroized(
#'   qx = c(0.015, 0.017, 0.019, 0.021, 0.024),
#'   i = 0.06,
#'   G = 19279,
#'   benefit = 1000000,
#'   e = 240
#' )
#'
#' @export
V_zeroized <- function(qx, i, G, benefit, r = 0, e = 0,
                       V_terminal = 0, floor_zero = TRUE) {
  qx <- .profit_check_probability(qx, "qx")
  i <- .profit_check_rate(i, "i")
  G <- .profit_check_nonnegative(G, "G")
  benefit <- .profit_check_nonnegative(benefit, "benefit")
  r <- .profit_check_probability(r, "r")
  e <- .profit_check_nonnegative(e, "e")
  V_terminal <- .profit_check_scalar(
    V_terminal, "V_terminal", nonnegative = TRUE
  )
  floor_zero <- .profit_check_logical_scalar(floor_zero, "floor_zero")

  n <- length(qx)
  i <- .profit_recycle_to(i, n, "i")
  G <- .profit_recycle_to(G, n, "G")
  benefit <- .profit_recycle_to(benefit, n, "benefit")
  r <- .profit_recycle_to(r, n, "r")
  e <- .profit_recycle_to(e, n, "e")

  V <- numeric(n + 1L)
  V[n + 1L] <- V_terminal

  for (k in seq.int(n, 1L)) {
    net_premium <- G[k] * (1 - r[k]) - e[k]
    reserve <- (
      benefit[k] * qx[k] +
        V[k + 1L] * (1 - qx[k])
    ) / (1 + i[k]) - net_premium

    if (floor_zero) {
      reserve <- max(reserve, 0)
    }

    V[k] <- reserve
  }

  names(V) <- paste0("V", 0:n)
  V
}
