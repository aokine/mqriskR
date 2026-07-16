# =========================================================
# Multiple-decrement applications
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.md_check_numeric <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    stop(name, " must contain finite numeric values.", call. = FALSE)
  }
  as.numeric(x)
}

#' @noRd
.md_check_nonnegative <- function(x, name) {
  x <- .md_check_numeric(x, name)
  if (any(x < 0)) stop(name, " must contain nonnegative values.", call. = FALSE)
  x
}

#' @noRd
.md_check_probability <- function(x, name) {
  x <- .md_check_numeric(x, name)
  if (any(x < 0 | x > 1)) stop(name, " must contain values in [0, 1].", call. = FALSE)
  x
}

#' @noRd
.md_check_interest <- function(i, name = "i") {
  i <- .md_check_numeric(i, name)
  if (any(i <= -1)) stop(name, " must contain values greater than -1.", call. = FALSE)
  i
}

#' @noRd
.md_check_scalar <- function(x, name, lower = -Inf, strict = FALSE) {
  x <- .md_check_numeric(x, name)
  if (length(x) != 1L) stop(name, " must be a finite numeric scalar.", call. = FALSE)
  bad <- if (strict) x <= lower else x < lower
  if (bad) stop(name, if (strict) " must be greater than " else " must be at least ", lower, ".", call. = FALSE)
  x
}

#' @noRd
.md_check_integer_scalar <- function(x, name, lower = 0L) {
  x <- .md_check_numeric(x, name)
  if (length(x) != 1L || x < lower || abs(x - round(x)) > 1e-10) {
    stop(name, " must be a single integer greater than or equal to ", lower, ".", call. = FALSE)
  }
  as.integer(round(x))
}

#' @noRd
.md_check_function <- function(x, name) {
  if (!is.function(x)) stop(name, " must be a function.", call. = FALSE)
  x
}

#' @noRd
.md_recycle <- function(..., .names = NULL) {
  values <- list(...)
  lengths <- vapply(values, length, integer(1))
  if (any(lengths == 0L)) stop("Inputs must have positive length.", call. = FALSE)
  common <- max(lengths)
  if (is.null(.names)) .names <- paste0("Argument ", seq_along(values))
  for (j in seq_along(values)) {
    if (!lengths[j] %in% c(1L, common)) {
      stop(.names[j], " must have length 1 or the common length ", common, ".", call. = FALSE)
    }
    values[[j]] <- rep(values[[j]], length.out = common)
  }
  values
}

#' @noRd
.md_check_grid <- function(t) {
  t <- .md_check_nonnegative(t, "t")
  if (length(t) < 2L) stop("t must contain at least two time points.", call. = FALSE)
  if (is.unsorted(t, strictly = TRUE)) stop("t must be strictly increasing.", call. = FALSE)
  t
}

#' @noRd
.md_grid <- function(n, h) {
  n <- .md_check_scalar(n, "n", lower = 0)
  h <- .md_check_scalar(h, "h", lower = 0, strict = TRUE)
  if (n == 0) return(0)
  out <- seq(0, n, by = h)
  if (tail(out, 1L) < n) out <- c(out, n)
  out[length(out)] <- n
  unique(out)
}

#' @noRd
.md_trapz <- function(t, y) {
  sum(diff(t) * (head(y, -1L) + tail(y, -1L)) / 2)
}

# -------------------------------------------------------------------------
# Insurance present values
# -------------------------------------------------------------------------

#' Discrete multiple-decrement insurance present value
#'
#' Computes the actuarial present value of a benefit payable at the end of
#' the year of decrement from a specified cause.
#'
#' @param qj Cause-specific decrement probabilities by policy year.
#' @param ptau In-force probabilities at the beginning of each policy year.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param benefit Benefit payable on decrement. May be scalar or vector.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @examples
#' qj <- c(0.02, 0.02, 0.02, 0.02, 0.02)
#' ptau <- c(1, 0.95, 0.89, 0.82, 0.74)
#' Axj_md(qj, ptau, i = 0.06, benefit = 1000)
#'
#' @export
Axj_md <- function(qj, ptau, i, benefit = 1) {
  qj <- .md_check_probability(qj, "qj")
  ptau <- .md_check_probability(ptau, "ptau")
  if (length(qj) != length(ptau)) stop("qj and ptau must have the same length.", call. = FALSE)
  i <- .md_check_interest(i)
  benefit <- .md_check_nonnegative(benefit, "benefit")
  vals <- .md_recycle(i, benefit, .names = c("i", "benefit"))
  times <- seq_along(qj)
  vapply(seq_along(vals[[1L]]), function(j) {
    sum(vals[[2L]][j] * (1 + vals[[1L]][j])^(-times) * ptau * qj)
  }, numeric(1))
}

#' Continuous multiple-decrement insurance present value
#'
#' Approximates the actuarial present value of a benefit payable at the moment
#' of decrement using the trapezoidal rule.
#'
#' @param t Strictly increasing nonnegative time points.
#' @param ptau In-force probabilities at the supplied time points.
#' @param muj Cause-specific decrement intensities at the supplied time points.
#' @param delta Force of interest. May be scalar or vector.
#' @param benefit Benefit payable on decrement. May be scalar or vector.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @examples
#' t <- seq(0, 20, by = 0.1)
#' ptau <- exp(-0.012 * t)
#' muj <- rep(0.002, length(t))
#' Abarxj_md(t, ptau, muj, delta = 0.05, benefit = 2000)
#'
#' @export
Abarxj_md <- function(t, ptau, muj, delta, benefit = 1) {
  t <- .md_check_grid(t)
  ptau <- .md_check_probability(ptau, "ptau")
  muj <- .md_check_nonnegative(muj, "muj")
  if (length(ptau) != length(t) || length(muj) != length(t)) {
    stop("t, ptau, and muj must have the same length.", call. = FALSE)
  }
  delta <- .md_check_numeric(delta, "delta")
  benefit <- .md_check_nonnegative(benefit, "benefit")
  vals <- .md_recycle(delta, benefit, .names = c("delta", "benefit"))
  vapply(seq_along(vals[[1L]]), function(j) {
    .md_trapz(t, vals[[2L]][j] * exp(-vals[[1L]][j] * t) * ptau * muj)
  }, numeric(1))
}

# -------------------------------------------------------------------------
# Asset shares
# -------------------------------------------------------------------------

#' General projected asset-share path
#'
#' Computes projected asset shares for a multiple-decrement contract. Rows of
#' \code{b_mat} and \code{q_mat} represent policy years and columns represent
#' decrement causes.
#'
#' @param AS0 Initial asset share.
#' @param G Premium amount by policy year.
#' @param r Percent-of-premium expense rate by policy year.
#' @param e Fixed expense by policy year.
#' @param b_mat Matrix of decrement benefits.
#' @param q_mat Matrix of decrement probabilities.
#' @param p_tau In-force probability by policy year.
#' @param i Effective annual interest rate by policy year.
#' @param b_surv Survival benefit by policy year.
#'
#' @return A data frame with policy year \code{k} and asset share \code{AS}.
#'
#' @export
AS_path_md <- function(AS0, G, r, e, b_mat, q_mat, p_tau, i, b_surv = NULL) {
  AS0 <- .md_check_scalar(AS0, "AS0")
  if (!is.matrix(b_mat) || !is.matrix(q_mat)) stop("b_mat and q_mat must be matrices.", call. = FALSE)
  if (!identical(dim(b_mat), dim(q_mat))) stop("b_mat and q_mat must have the same dimensions.", call. = FALSE)
  if (nrow(b_mat) == 0L || ncol(b_mat) == 0L) stop("b_mat and q_mat must be nonempty.", call. = FALSE)
  storage.mode(b_mat) <- "double"
  storage.mode(q_mat) <- "double"
  if (any(!is.finite(b_mat)) || any(b_mat < 0)) stop("b_mat must contain nonnegative finite values.", call. = FALSE)
  if (any(!is.finite(q_mat)) || any(q_mat < 0) || any(q_mat > 1)) stop("q_mat must contain values in [0, 1].", call. = FALSE)
  if (any(rowSums(q_mat) > 1 + 1e-12)) stop("Each row sum of q_mat must not exceed 1.", call. = FALSE)

  n_years <- nrow(b_mat)
  if (is.null(b_surv)) b_surv <- 0
  vals <- .md_recycle(
    .md_check_nonnegative(G, "G"),
    .md_check_probability(r, "r"),
    .md_check_nonnegative(e, "e"),
    .md_check_probability(p_tau, "p_tau"),
    .md_check_interest(i),
    .md_check_nonnegative(b_surv, "b_surv"),
    .names = c("G", "r", "e", "p_tau", "i", "b_surv")
  )
  common <- length(vals[[1L]])
  if (common == 1L) vals <- lapply(vals, rep, length.out = n_years)
  if (length(vals[[1L]]) != n_years) stop("Yearly inputs must have length 1 or nrow(b_mat).", call. = FALSE)
  if (any(vals[[4L]] <= 0)) stop("p_tau must be positive in every policy year.", call. = FALSE)

  AS <- numeric(n_years + 1L)
  AS[1L] <- AS0
  for (k in seq_len(n_years)) {
    claims <- sum(b_mat[k, ] * q_mat[k, ])
    AS[k + 1L] <- ((AS[k] + vals[[1L]][k] * (1 - vals[[2L]][k]) - vals[[3L]][k]) *
                     (1 + vals[[5L]][k]) - claims) / vals[[4L]][k] - vals[[6L]][k]
  }
  data.frame(k = 0:n_years, AS = AS)
}

#' Projected asset-share path for two decrement causes
#'
#' Convenience interface to \code{AS_path_md()} for two causes.
#'
#' @param AS0 Initial asset share.
#' @param G Premium amount by policy year.
#' @param r Percent-of-premium expense rate by policy year.
#' @param e Fixed expense by policy year.
#' @param b1 Cause 1 benefit by policy year.
#' @param b2 Cause 2 benefit by policy year.
#' @param q1 Cause 1 probability by policy year.
#' @param q2 Cause 2 probability by policy year.
#' @param p_tau In-force probability by policy year.
#' @param i Effective annual interest rate by policy year.
#' @param b3 Survival benefit by policy year.
#'
#' @return A data frame with policy year \code{k} and asset share \code{AS}.
#'
#' @export
AS_path <- function(AS0, G, r, e, b1, b2, q1, q2, p_tau, i, b3 = NULL) {
  vals <- .md_recycle(
    .md_check_nonnegative(b1, "b1"),
    .md_check_nonnegative(b2, "b2"),
    .md_check_probability(q1, "q1"),
    .md_check_probability(q2, "q2"),
    .names = c("b1", "b2", "q1", "q2")
  )
  if (any(vals[[3L]] + vals[[4L]] > 1 + 1e-12)) stop("q1 + q2 must not exceed 1.", call. = FALSE)
  AS_path_md(
    AS0 = AS0, G = G, r = r, e = e,
    b_mat = cbind(cause1 = vals[[1L]], cause2 = vals[[2L]]),
    q_mat = cbind(cause1 = vals[[3L]], cause2 = vals[[4L]]),
    p_tau = p_tau, i = i, b_surv = b3
  )
}

# -------------------------------------------------------------------------
# Disability-state probabilities and premiums
# -------------------------------------------------------------------------

#' Euler approximation of disability-state probabilities
#'
#' Approximates probabilities of being healthy, disabled, or deceased in a
#' three-state model that allows recovery from disability. The final time is
#' always included, with a shorter final step when needed.
#'
#' @param h Positive step size.
#' @param n Nonnegative final time.
#' @param mu01 Healthy-to-disabled intensity function.
#' @param mu02 Healthy-to-deceased intensity function.
#' @param mu10 Disabled-to-healthy intensity function.
#' @param mu12 Disabled-to-deceased intensity function.
#' @param p00_0 Initial healthy-state probability.
#' @param p01_0 Initial disabled-state probability.
#'
#' @return A data frame with columns \code{t}, \code{tp00}, \code{tp01},
#'   and \code{tp02}.
#'
#' @examples
#' mu01 <- function(t) 0.10 * t + 0.20
#' mu02 <- function(t) 0.20
#' mu10 <- function(t) 0.50
#' mu12 <- function(t) 0.125 * t + 0.20
#' tp00_tp01_euler(0.10, 2, mu01, mu02, mu10, mu12)
#'
#' @export
tp00_tp01_euler <- function(h, n, mu01, mu02, mu10, mu12,
                            p00_0 = 1, p01_0 = 0) {
  grid <- .md_grid(n, h)
  mu01 <- .md_check_function(mu01, "mu01")
  mu02 <- .md_check_function(mu02, "mu02")
  mu10 <- .md_check_function(mu10, "mu10")
  mu12 <- .md_check_function(mu12, "mu12")
  p00_0 <- .md_check_probability(p00_0, "p00_0")
  p01_0 <- .md_check_probability(p01_0, "p01_0")
  if (length(p00_0) != 1L || length(p01_0) != 1L) stop("p00_0 and p01_0 must be scalar.", call. = FALSE)
  if (p00_0 + p01_0 > 1 + 1e-12) stop("p00_0 + p01_0 must not exceed 1.", call. = FALSE)

  p00 <- numeric(length(grid)); p01 <- numeric(length(grid))
  p00[1L] <- p00_0; p01[1L] <- p01_0
  if (length(grid) > 1L) {
    for (idx in seq_len(length(grid) - 1L)) {
      time <- grid[idx]; step <- grid[idx + 1L] - grid[idx]
      rates <- c(mu01(time), mu02(time), mu10(time), mu12(time))
      if (length(rates) != 4L || any(!is.finite(rates)) || any(rates < 0)) {
        stop("Intensity functions must return one nonnegative finite value.", call. = FALSE)
      }
      dp00 <- p01[idx] * rates[3L] - p00[idx] * (rates[1L] + rates[2L])
      dp01 <- p00[idx] * rates[1L] - p01[idx] * (rates[3L] + rates[4L])
      p00[idx + 1L] <- p00[idx] + step * dp00
      p01[idx + 1L] <- p01[idx] + step * dp01
    }
  }
  data.frame(t = grid, tp00 = p00, tp01 = p01, tp02 = 1 - p00 - p01)
}

#' Continuous premium approximation in a disability model
#'
#' Approximates a continuous premium rate by trapezoidal integration.
#'
#' @param t Strictly increasing nonnegative time points.
#' @param tp00 Healthy-state probabilities.
#' @param tp01 Disabled-state probabilities.
#' @param delta Force of interest.
#' @param mu02 Healthy-to-deceased intensity function.
#' @param mu12 Disabled-to-deceased intensity function.
#' @param B02 Benefit on death while healthy.
#' @param B12 Benefit on death while disabled.
#' @param R Continuous disability income rate.
#'
#' @return A numeric scalar.
#'
#' @examples
#' mu01 <- function(t) 0.10 * t + 0.20
#' mu02 <- function(t) 0.20
#' mu10 <- function(t) 0.50
#' mu12 <- function(t) 0.125 * t + 0.20
#'
#' probs <- tp00_tp01_euler(
#'   h = 0.10,
#'   n = 2,
#'   mu01 = mu01,
#'   mu02 = mu02,
#'   mu10 = mu10,
#'   mu12 = mu12
#' )
#'
#' Pbar_trapz_ms(
#'   t = probs$t,
#'   tp00 = probs$tp00,
#'   tp01 = probs$tp01,
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
                          mu02, mu12, B02 = 1, B12 = 1, R = 0) {
  t <- .md_check_grid(t)
  tp00 <- .md_check_probability(tp00, "tp00")
  tp01 <- .md_check_probability(tp01, "tp01")
  if (length(tp00) != length(t) || length(tp01) != length(t)) stop("t, tp00, and tp01 must have the same length.", call. = FALSE)
  if (any(tp00 + tp01 > 1 + 1e-10)) stop("tp00 + tp01 must not exceed 1.", call. = FALSE)
  delta <- .md_check_scalar(delta, "delta")
  B02 <- .md_check_scalar(B02, "B02", lower = 0)
  B12 <- .md_check_scalar(B12, "B12", lower = 0)
  R <- .md_check_scalar(R, "R", lower = 0)
  mu02 <- .md_check_function(mu02, "mu02")
  mu12 <- .md_check_function(mu12, "mu12")
  m02 <- rep(as.numeric(mu02(t)), length.out = length(t))
  m12 <- rep(as.numeric(mu12(t)), length.out = length(t))
  if (!length(mu02(t)) %in% c(1L, length(t)) || !length(mu12(t)) %in% c(1L, length(t))) stop("mu02 and mu12 must return length 1 or length(t).", call. = FALSE)
  if (any(!is.finite(m02)) || any(m02 < 0) || any(!is.finite(m12)) || any(m12 < 0)) stop("mu02 and mu12 must return nonnegative finite values.", call. = FALSE)
  discount <- exp(-delta * t)
  num <- .md_trapz(t, discount * (tp00 * m02 * B02 + tp01 * m12 * B12 + tp01 * R))
  den <- .md_trapz(t, discount * tp00)
  if (!is.finite(den) || den <= 0) stop("The premium present-value denominator must be positive.", call. = FALSE)
  num / den
}

# -------------------------------------------------------------------------
# Thiele reserve equations
# -------------------------------------------------------------------------

#' Reserve derivatives for a disability model with recovery
#'
#' Computes the coupled Thiele reserve derivatives. Numeric arguments may be
#' scalar or compatible vectors.
#'
#' @param t Time.
#' @param V0 Healthy-state reserve.
#' @param V1 Disabled-state reserve.
#' @param delta Force of interest.
#' @param Pbar Continuous premium rate.
#' @param B Death benefit.
#' @param R Continuous disability income rate.
#' @param mu01 Healthy-to-disabled intensity function.
#' @param mu02 Healthy-to-deceased intensity function.
#' @param mu10 Disabled-to-healthy intensity function.
#' @param mu12 Disabled-to-deceased intensity function.
#'
#' @return A named vector for scalar input or a two-column matrix for
#'   vectorized input.
#'
#' @export
thiele_dVdt_01 <- function(t, V0, V1, delta, Pbar, B, R,
                           mu01, mu02, mu10, mu12) {
  vals <- .md_recycle(
    .md_check_nonnegative(t, "t"), .md_check_numeric(V0, "V0"),
    .md_check_numeric(V1, "V1"), .md_check_numeric(delta, "delta"),
    .md_check_numeric(Pbar, "Pbar"), .md_check_nonnegative(B, "B"),
    .md_check_nonnegative(R, "R"),
    .names = c("t", "V0", "V1", "delta", "Pbar", "B", "R")
  )
  t <- vals[[1L]]; V0 <- vals[[2L]]; V1 <- vals[[3L]]
  delta <- vals[[4L]]; Pbar <- vals[[5L]]; B <- vals[[6L]]; R <- vals[[7L]]
  funs <- list(mu01 = mu01, mu02 = mu02, mu10 = mu10, mu12 = mu12)
  rates <- lapply(names(funs), function(nm) {
    fun <- .md_check_function(funs[[nm]], nm)
    out <- as.numeric(fun(t))
    if (!length(out) %in% c(1L, length(t))) stop(nm, " must return length 1 or the common length.", call. = FALSE)
    out <- rep(out, length.out = length(t))
    if (any(!is.finite(out)) || any(out < 0)) stop(nm, " must return nonnegative finite values.", call. = FALSE)
    out
  })
  dV0 <- Pbar + delta * V0 - rates[[2L]] * (B - V0) - rates[[1L]] * (V1 - V0)
  dV1 <- delta * V1 - R - rates[[4L]] * (B - V1) - rates[[3L]] * (V0 - V1)
  if (length(dV0) == 1L) return(c(dV0 = dV0, dV1 = dV1))
  cbind(dV0 = dV0, dV1 = dV1)
}

#' Backward reserve path for a disability model with recovery
#'
#' Computes a backward Euler reserve path from terminal healthy-state and
#' disabled-state reserves.
#'
#' @inheritParams thiele_dVdt_01
#' @param h Positive step size.
#' @param n Nonnegative final time.
#' @param V0_n Terminal healthy-state reserve.
#' @param V1_n Terminal disabled-state reserve.
#'
#' @return A data frame with columns \code{t}, \code{tV0}, and \code{tV1}.
#'
#' @export
thiele_path_01 <- function(h, n, delta, Pbar, B, R,
                           mu01, mu02, mu10, mu12,
                           V0_n = 0, V1_n = 0) {
  grid <- .md_grid(n, h); back <- rev(grid)
  delta <- .md_check_scalar(delta, "delta")
  Pbar <- .md_check_scalar(Pbar, "Pbar")
  B <- .md_check_scalar(B, "B", lower = 0)
  R <- .md_check_scalar(R, "R", lower = 0)
  V0_n <- .md_check_scalar(V0_n, "V0_n")
  V1_n <- .md_check_scalar(V1_n, "V1_n")
  V0 <- numeric(length(back)); V1 <- numeric(length(back))
  V0[1L] <- V0_n; V1[1L] <- V1_n
  if (length(back) > 1L) {
    for (idx in seq_len(length(back) - 1L)) {
      step <- back[idx] - back[idx + 1L]
      dV <- thiele_dVdt_01(back[idx], V0[idx], V1[idx], delta, Pbar, B, R,
                           mu01, mu02, mu10, mu12)
      V0[idx + 1L] <- V0[idx] - step * unname(dV["dV0"])
      V1[idx + 1L] <- V1[idx] - step * unname(dV["dV1"])
    }
  }
  data.frame(t = grid, tV0 = rev(V0), tV1 = rev(V1))
}

# -------------------------------------------------------------------------
# Discrete-time Markov chains
# -------------------------------------------------------------------------

#' Multi-step transition probability
#'
#' Computes an entry of the matrix power \eqn{P^n}.
#'
#' @param P Square transition-probability matrix.
#' @param n Nonnegative integer number of steps.
#' @param i Starting-state index.
#' @param j Ending-state index.
#'
#' @return A numeric scalar.
#'
#' @examples
#' P <- matrix(c(0.9, 0.1, 0, 1), nrow = 2, byrow = TRUE)
#' markov_nstep_prob(P, n = 3, i = 1, j = 2)
#'
#' @export
markov_nstep_prob <- function(P, n, i, j) {
  if (!is.matrix(P) || !is.numeric(P)) stop("P must be a numeric matrix.", call. = FALSE)
  if (nrow(P) == 0L || nrow(P) != ncol(P)) stop("P must be a nonempty square matrix.", call. = FALSE)
  if (any(!is.finite(P)) || any(P < 0) || any(P > 1)) stop("P must contain values in [0, 1].", call. = FALSE)
  if (any(abs(rowSums(P) - 1) > 1e-10)) stop("Each row of P must sum to 1.", call. = FALSE)
  n <- .md_check_integer_scalar(n, "n", 0L)
  i <- .md_check_integer_scalar(i, "i", 1L)
  j <- .md_check_integer_scalar(j, "j", 1L)
  if (i > nrow(P) || j > ncol(P)) stop("i and j must be valid state indices.", call. = FALSE)
  Pn <- diag(nrow(P))
  if (n > 0L) for (step in seq_len(n)) Pn <- Pn %*% P
  unname(Pn[i, j])
}

# -------------------------------------------------------------------------
# Multiple-decrement gain or loss
# -------------------------------------------------------------------------

#' Gain or loss in a two-cause multiple-decrement model
#'
#' Computes one-year gain or loss under simultaneous within-year decrement
#' probabilities or an ordered case in which Cause 2 occurs at year-end.
#' Numeric arguments may be scalar or compatible vectors.
#'
#' @param Vt Gross reserve at the beginning of the year.
#' @param G Gross premium.
#' @param r Percent-of-premium expense rate.
#' @param e Fixed beginning-of-year expense.
#' @param i Earned effective annual interest rate.
#' @param b1 Cause 1 benefit.
#' @param b2 Cause 2 benefit.
#' @param s1 Cause 1 settlement expense.
#' @param s2 Cause 2 settlement expense.
#' @param q1 Cause 1 decrement probability.
#' @param q2 Cause 2 decrement probability.
#' @param Vt1 Gross reserve at the end of the year.
#' @param year_end_cause2 Whether Cause 2 occurs only at year-end.
#' @param q1prime Single-decrement Cause 1 probability for the ordered case.
#' @param q2prime Single-decrement Cause 2 probability for the ordered case.
#'
#' @return A numeric vector of gains or losses.
#'
#' @examples
#' gain_loss_md(
#'   Vt = 115, G = 16, r = 0, e = 3, i = 0.06,
#'   b1 = 1000, b2 = 110, q1 = 0.01, q2 = 0.10, Vt1 = 128.83
#' )
#'
#' @export
gain_loss_md <- function(Vt, G, r, e, i,
                         b1, b2, s1 = 0, s2 = 0,
                         q1, q2, Vt1,
                         year_end_cause2 = FALSE,
                         q1prime = NULL, q2prime = NULL) {
  if (!is.logical(year_end_cause2) || length(year_end_cause2) == 0L || anyNA(year_end_cause2)) {
    stop("year_end_cause2 must be a non-missing logical value or vector.", call. = FALSE)
  }
  q1prime <- if (is.null(q1prime)) NA_real_ else .md_check_probability(q1prime, "q1prime")
  q2prime <- if (is.null(q2prime)) NA_real_ else .md_check_probability(q2prime, "q2prime")
  vals <- .md_recycle(
    .md_check_numeric(Vt, "Vt"), .md_check_nonnegative(G, "G"),
    .md_check_probability(r, "r"), .md_check_nonnegative(e, "e"),
    .md_check_interest(i), .md_check_nonnegative(b1, "b1"),
    .md_check_nonnegative(b2, "b2"), .md_check_nonnegative(s1, "s1"),
    .md_check_nonnegative(s2, "s2"), .md_check_probability(q1, "q1"),
    .md_check_probability(q2, "q2"), .md_check_numeric(Vt1, "Vt1"),
    year_end_cause2, q1prime, q2prime,
    .names = c("Vt", "G", "r", "e", "i", "b1", "b2", "s1", "s2",
               "q1", "q2", "Vt1", "year_end_cause2", "q1prime", "q2prime")
  )
  Vt <- vals[[1L]]; G <- vals[[2L]]; r <- vals[[3L]]; e <- vals[[4L]]
  i <- vals[[5L]]; b1 <- vals[[6L]]; b2 <- vals[[7L]]; s1 <- vals[[8L]]
  s2 <- vals[[9L]]; q1 <- vals[[10L]]; q2 <- vals[[11L]]; Vt1 <- vals[[12L]]
  ordered <- vals[[13L]]; q1prime <- vals[[14L]]; q2prime <- vals[[15L]]
  if (any(!ordered & q1 + q2 > 1 + 1e-12)) stop("q1 + q2 must not exceed 1.", call. = FALSE)
  if (any(ordered & (is.na(q1prime) | is.na(q2prime)))) {
    stop("q1prime and q2prime must be provided when year_end_cause2 = TRUE.", call. = FALSE)
  }
  lhs <- (Vt + G * (1 - r) - e) * (1 + i)
  rhs <- numeric(length(lhs))
  simultaneous <- !ordered
  rhs[simultaneous] <- (b1[simultaneous] + s1[simultaneous]) * q1[simultaneous] +
    (b2[simultaneous] + s2[simultaneous]) * q2[simultaneous] +
    (1 - q1[simultaneous] - q2[simultaneous]) * Vt1[simultaneous]
  rhs[ordered] <- (b1[ordered] + s1[ordered]) * q1prime[ordered] +
    (b2[ordered] + s2[ordered]) * (1 - q1prime[ordered]) * q2prime[ordered] +
    (1 - q1prime[ordered]) * (1 - q2prime[ordered]) * Vt1[ordered]
  lhs - rhs
}
