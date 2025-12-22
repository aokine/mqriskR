#' Convert between i, d, and delta (compound interest)
#'
#' Provides consistent conversions between:
#' - effective interest rate i
#' - effective discount rate d
#' - force of interest delta
#'
#' Exactly one of i, d, or delta must be provided.
#'
#' @param i Effective interest rate.
#' @param d Effective discount rate.
#' @param delta Force of interest.
#' @param m Optional compounding frequency for nominal rate i^(m).
#'
#' @return A list with elements i, d, delta and (if m is supplied) im (nominal).
#' @export
#'
#' @examples
#' interest_convert(i = 0.05)
#' interest_convert(d = 0.04761905)
#' interest_convert(delta = log(1.05))
interest_convert <- function(i = NULL, d = NULL, delta = NULL, m = NULL) {
  supplied <- sum(!is.null(i), !is.null(d), !is.null(delta))
  if (supplied != 1) {
    stop("Provide exactly one of i, d, or delta.", call. = FALSE)
  }

  if (!is.null(i)) {
    if (i <= -1) stop("i must be > -1.", call. = FALSE)
    d <- i / (1 + i)
    delta <- log(1 + i)
  } else if (!is.null(d)) {
    if (d >= 1) stop("d must be < 1.", call. = FALSE)
    if (d < 0) stop("d must be >= 0 for standard discounting.", call. = FALSE)
    i <- d / (1 - d)
    delta <- log(1 + i)
  } else {
    # delta supplied
    i <- exp(delta) - 1
    d <- i / (1 + i)
  }

  out <- list(i = i, d = d, delta = delta)

  if (!is.null(m)) {
    if (length(m) != 1 || m <= 0 || m != as.integer(m)) {
      stop("m must be a positive integer.", call. = FALSE)
    }
    # Nominal annual rate convertible m-thly:
    # i^(m) = m * ((1+i)^(1/m) - 1)
    out$im <- m * ((1 + i)^(1 / m) - 1)
    out$m <- m
  }

  out
}

#' Discount factor (1+i)^(-t)
#'
#' @param i Effective interest rate.
#' @param t Time (can be vector).
#' @return Discount factor v^t.
#' @export
#' @examples
#' discount(0.05, 0:5)
discount <- function(i, t) {
  if (any(i <= -1)) stop("i must be > -1.", call. = FALSE)
  (1 + i)^(-t)
}

#' Present value of cash flows at time 0
#'
#' @param cf Cash flow amounts (positive = inflow, negative = outflow).
#' @param t Times of cash flows (same length as cf).
#' @param i Effective interest rate.
#' @return Present value at time 0.
#' @export
#' @examples
#' pv_cashflows(c(-100, 60, 60), c(0, 1, 2), i = 0.10)
pv_cashflows <- function(cf, t, i) {
  if (length(cf) != length(t)) stop("cf and t must have the same length.", call. = FALSE)
  sum(cf * discount(i, t))
}

#' Solve yield rate (IRR) via an equation of value
#'
#' Finds i such that pv_cashflows(cf, t, i) = 0.
#'
#' @param cf Cash flows.
#' @param t Times.
#' @param interval Two-length numeric vector bracketing the root.
#' @param tol Tolerance passed to uniroot.
#' @return Yield rate i.
#' @export
#' @examples
#' solve_yield(c(-100, 60, 60), c(0, 1, 2), interval = c(-0.5, 1))
solve_yield <- function(cf, t, interval = c(-0.99, 1), tol = 1e-10) {
  f <- function(i) pv_cashflows(cf, t, i)

  a <- interval[1]
  b <- interval[2]
  fa <- f(a)
  fb <- f(b)

  if (!is.finite(fa) || !is.finite(fb)) {
    stop("NPV not finite at interval endpoints. Adjust interval.", call. = FALSE)
  }
  if (fa == 0) return(a)
  if (fb == 0) return(b)
  if (fa * fb > 0) {
    stop("No sign change in NPV over interval; cannot bracket root. Try a wider interval.", call. = FALSE)
  }

  stats::uniroot(f, interval = interval, tol = tol)$root
}

#' Level annuity-certain present value
#'
#' Computes PV of an n-period annuity certain with level payments of 1 per period.
#'
#' @param n Number of payments/periods.
#' @param i Effective interest rate per period.
#' @param due If TRUE, annuity-due; otherwise annuity-immediate.
#' @param m Payment frequency per period (m=1 means annual).
#' @param cont If TRUE, continuous payment model.
#' @return Present value.
#' @export
#' @examples
#' annuity_certain(n = 10, i = 0.05)
#' annuity_certain(n = 10, i = 0.05, due = TRUE)
#' annuity_certain(n = 10, i = 0.05, cont = TRUE)
annuity_certain <- function(n, i, due = FALSE, m = 1, cont = FALSE) {
  if (length(n) != 1 || n < 0) stop("n must be a nonnegative scalar.", call. = FALSE)
  if (length(m) != 1 || m <= 0 || m != as.integer(m)) stop("m must be a positive integer.", call. = FALSE)
  if (i <= -1) stop("i must be > -1.", call. = FALSE)

  if (cont) {
    delta <- log(1 + i)
    if (delta == 0) return(n)  # limit as delta -> 0
    return((1 - exp(-delta * n)) / delta)
  }

  # m-thly: treat as nm payments of 1/m each at spacing 1/m
  nm <- n * m
  if (nm == 0) return(0)

  j <- (1 + i)^(1 / m) - 1  # effective per subperiod
  v <- 1 / (1 + j)

  if (!due) {
    # immediate: payments at times 1/m,...,nm/m
    pv <- (1 / m) * sum(v^(1:nm))
  } else {
    # due: payments at times 0,...,(nm-1)/m
    pv <- (1 / m) * sum(v^(0:(nm - 1)))
  }

  pv
}
