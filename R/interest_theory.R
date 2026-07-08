# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.recycle2_interest <- function(a, b, a_name = "a", b_name = "b") {
  n <- max(length(a), length(b))

  if (!length(a) %in% c(1L, n)) {
    stop(a_name, " must have length 1 or common length.", call. = FALSE)
  }
  if (!length(b) %in% c(1L, n)) {
    stop(b_name, " must have length 1 or common length.", call. = FALSE)
  }

  if (length(a) == 1L) a <- rep(a, n)
  if (length(b) == 1L) b <- rep(b, n)

  list(a = a, b = b)
}

#' Convert between compound-interest quantities
#'
#' Provides consistent conversions between effective interest rate, effective
#' discount rate, force of interest, and optional nominal interest rate
#' convertible m-thly.
#'
#' Exactly one of `i`, `d`, or `delta` must be provided.
#'
#' @param i Effective interest rate. May be scalar or vector.
#' @param d Effective discount rate. May be scalar or vector.
#' @param delta Force of interest. May be scalar or vector.
#' @param m Optional positive integer compounding frequency for the nominal
#'   rate convertible m-thly.
#'
#' @return A list with elements `i`, `d`, `delta`, and, if `m` is supplied,
#'   `im` and `m`.
#' @export
#'
#' @examples
#' interest_convert(i = 0.05)
#' interest_convert(i = c(0.03, 0.05, 0.07))
#' interest_convert(d = 0.04761905)
#' interest_convert(delta = log(1.05))
interest_convert <- function(i = NULL, d = NULL, delta = NULL, m = NULL) {
  supplied <- sum(!is.null(i), !is.null(d), !is.null(delta))
  if (supplied != 1) {
    stop("Provide exactly one of i, d, or delta.", call. = FALSE)
  }

  if (!is.null(i)) {
    i <- as.numeric(i)
    if (length(i) == 0L || any(!is.finite(i))) {
      stop("i must contain finite numeric values.", call. = FALSE)
    }
    if (any(i <= -1)) {
      stop("i must be greater than -1.", call. = FALSE)
    }

    d <- i / (1 + i)
    delta <- log(1 + i)

  } else if (!is.null(d)) {
    d <- as.numeric(d)
    if (length(d) == 0L || any(!is.finite(d))) {
      stop("d must contain finite numeric values.", call. = FALSE)
    }
    if (any(d >= 1)) {
      stop("d must be less than 1.", call. = FALSE)
    }
    if (any(d < 0)) {
      stop("d must be nonnegative for standard discounting.", call. = FALSE)
    }

    i <- d / (1 - d)
    delta <- log(1 + i)

  } else {
    delta <- as.numeric(delta)
    if (length(delta) == 0L || any(!is.finite(delta))) {
      stop("delta must contain finite numeric values.", call. = FALSE)
    }

    i <- exp(delta) - 1
    d <- i / (1 + i)
  }

  out <- list(i = i, d = d, delta = delta)

  if (!is.null(m)) {
    m <- as.numeric(m)
    if (length(m) != 1L || !is.finite(m) || m <= 0 ||
        abs(m - round(m)) > 1e-10) {
      stop("m must be a positive integer.", call. = FALSE)
    }

    m <- as.integer(round(m))

    out$im <- m * ((1 + i)^(1 / m) - 1)
    out$m <- m
  }

  out
}

#' Discount factor for compound interest
#'
#' Computes the discount factor \eqn{v^t = (1+i)^{-t}}.
#'
#' @param i Effective interest rate. May be scalar or vector.
#' @param t Time. May be scalar or vector.
#'
#' @return Numeric vector of discount factors.
#' @export
#'
#' @examples
#' discount(0.05, 0:5)
#' discount(c(0.03, 0.05), 1)
discount <- function(i, t) {
  i <- as.numeric(i)
  t <- as.numeric(t)

  if (length(i) == 0L || length(t) == 0L) {
    stop("i and t must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(i)) || any(!is.finite(t))) {
    stop("i and t must contain finite numeric values.", call. = FALSE)
  }
  if (any(i <= -1)) {
    stop("i must be greater than -1.", call. = FALSE)
  }

  it <- .recycle2_interest(i, t, "i", "t")
  i <- it$a
  t <- it$b

  (1 + i)^(-t)
}

#' Present value of cash flows at time 0
#'
#' @param cf Cash flow amounts.
#' @param t Times of cash flows.
#' @param i Effective interest rate.
#'
#' @return Present value at time 0.
#' @export
#'
#' @examples
#' pv_cashflows(c(-100, 60, 60), c(0, 1, 2), i = 0.10)
pv_cashflows <- function(cf, t, i) {
  cf <- as.numeric(cf)
  t <- as.numeric(t)
  i <- as.numeric(i)

  if (length(cf) != length(t)) {
    stop("cf and t must have the same length.", call. = FALSE)
  }
  if (length(cf) == 0L) {
    stop("cf and t must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(cf)) || any(!is.finite(t))) {
    stop("cf and t must contain finite numeric values.", call. = FALSE)
  }
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single effective interest rate greater than -1.", call. = FALSE)
  }

  sum(cf * discount(i, t))
}

#' Solve the yield rate by the equation of value
#'
#' Finds the interest rate `i` such that the present value of the cash flows is 0.
#'
#' @param cf Cash flows.
#' @param t Times.
#' @param interval Two-length numeric vector bracketing the root.
#' @param tol Tolerance passed to `uniroot()`.
#'
#' @return Yield rate `i`.
#' @export
#'
#' @examples
#' solve_yield(c(-100, 60, 60), c(0, 1, 2), interval = c(-0.5, 1))
solve_yield <- function(cf, t, interval = c(-0.99, 1), tol = 1e-10) {
  cf <- as.numeric(cf)
  t <- as.numeric(t)
  interval <- as.numeric(interval)

  if (length(cf) != length(t)) {
    stop("cf and t must have the same length.", call. = FALSE)
  }
  if (length(interval) != 2L || any(!is.finite(interval))) {
    stop("interval must be a finite numeric vector of length 2.", call. = FALSE)
  }
  if (interval[1] <= -1 || interval[1] >= interval[2]) {
    stop("interval must satisfy -1 < interval[1] < interval[2].", call. = FALSE)
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive numeric scalar.", call. = FALSE)
  }

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
    stop("No sign change in NPV over interval; cannot bracket root. Try a wider interval.",
         call. = FALSE)
  }

  stats::uniroot(f, interval = interval, tol = tol)$root
}

#' Present value of a level annuity-certain
#'
#' Computes the present value of an n-period annuity certain with level
#' payments of 1 per period.
#'
#' @param n Number of payments or periods. May be scalar or vector.
#' @param i Effective interest rate per period. May be scalar or vector.
#' @param due If `TRUE`, annuity-due; otherwise annuity-immediate.
#' @param m Payment frequency per period. `m = 1` means annual.
#' @param cont If `TRUE`, use the continuous payment model.
#'
#' @return Numeric vector of present values.
#' @export
#'
#' @examples
#' annuity_certain(n = 10, i = 0.05)
#' annuity_certain(n = c(5, 10), i = 0.05)
#' annuity_certain(n = 10, i = c(0.03, 0.05))
#' annuity_certain(n = 10, i = 0.05, due = TRUE)
#' annuity_certain(n = 10, i = 0.05, cont = TRUE)
annuity_certain <- function(n, i, due = FALSE, m = 1, cont = FALSE) {
  n <- as.numeric(n)
  i <- as.numeric(i)
  m <- as.numeric(m)

  if (length(n) == 0L || any(!is.finite(n)) || any(n < 0)) {
    stop("n must contain nonnegative finite values.", call. = FALSE)
  }
  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain effective interest rates greater than -1.", call. = FALSE)
  }
  if (length(m) != 1L || !is.finite(m) || m <= 0 ||
      abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(due) || length(due) != 1L || is.na(due)) {
    stop("due must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(cont) || length(cont) != 1L || is.na(cont)) {
    stop("cont must be TRUE or FALSE.", call. = FALSE)
  }

  ni <- .recycle2_interest(n, i, "n", "i")
  n <- ni$a
  i <- ni$b
  m <- as.integer(round(m))

  out <- numeric(length(n))

  for (j in seq_along(n)) {
    if (cont) {
      delta <- log(1 + i[j])
      out[j] <- if (abs(delta) < 1e-14) {
        n[j]
      } else {
        (1 - exp(-delta * n[j])) / delta
      }
      next
    }

    nm <- n[j] * m
    if (abs(nm - round(nm)) > 1e-10) {
      stop("Each value of n * m must be an integer for m-thly annuity payments.",
           call. = FALSE)
    }

    nm <- as.integer(round(nm))
    if (nm == 0L) {
      out[j] <- 0
      next
    }

    sub_i <- (1 + i[j])^(1 / m) - 1
    v <- 1 / (1 + sub_i)

    out[j] <- if (!due) {
      (1 / m) * sum(v^(1:nm))
    } else {
      (1 / m) * sum(v^(0:(nm - 1L)))
    }
  }

  out
}
