# Chapter 13 extension: multiple-decrement models
# Suggested file: R/decrement_ch13_ext.R

# ---------------------------
# Basic validators
# ---------------------------

.check_prob_vec <- function(x, nm = "x") {
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop(nm, " must be finite.", call. = FALSE)
  if (any(x < 0 | x > 1)) stop(nm, " must lie in [0, 1].", call. = FALSE)
  x
}

.check_nonneg_scalar <- function(x, nm = "x") {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stop(nm, " must be a single nonnegative finite number.", call. = FALSE)
  }
  x
}

# ---------------------------
# Discrete multiple-decrement functions
# ---------------------------

#' Total probability of decrement q_x^(tau)
#'
#' @param qxj Numeric vector of cause-specific decrement probabilities.
#' @return Numeric scalar.
#' @examples
#' qxtau(c(0.011, 0.100))
#' @export
qxtau <- function(qxj) {
  qxj <- .check_prob_vec(qxj, "qxj")
  out <- sum(qxj)
  if (out > 1) stop("sum(qxj) cannot exceed 1.", call. = FALSE)
  out
}

#' Survival probability p_x^(tau)
#'
#' @param qxj Numeric vector of cause-specific decrement probabilities.
#' @return Numeric scalar.
#' @examples
#' pxtau(c(0.011, 0.100))
#' @export
pxtau <- function(qxj) {
  1 - qxtau(qxj)
}

#' Cause-specific decrements d_x^(j)
#'
#' @param lxtau Number alive at age x in the multiple-decrement table.
#' @param qxj Numeric vector of cause-specific decrement probabilities.
#' @return Numeric vector.
#' @examples
#' dxj(1000, c(0.011, 0.100))
#' @export
dxj <- function(lxtau, qxj) {
  lxtau <- .check_nonneg_scalar(lxtau, "lxtau")
  qxj <- .check_prob_vec(qxj, "qxj")
  if (sum(qxj) > 1) stop("sum(qxj) cannot exceed 1.", call. = FALSE)
  lxtau * qxj
}

#' Total decrements d_x^(tau)
#'
#' @param lxtau Number alive at age x in the multiple-decrement table.
#' @param qxj Numeric vector of cause-specific decrement probabilities.
#' @return Numeric scalar.
#' @examples
#' dxtau(1000, c(0.011, 0.100))
#' @export
dxtau <- function(lxtau, qxj) {
  sum(dxj(lxtau = lxtau, qxj = qxj))
}

#' Build a multiple-decrement table
#'
#' @param x Integer vector of ages or durations.
#' @param qxj Matrix/data.frame of cause-specific decrement probabilities.
#'   Rows correspond to ages in x, columns correspond to causes.
#' @param radix Starting l_x^(tau).
#' @return Data frame containing q^(j), q^(tau), p^(tau), l^(tau), d^(j), d^(tau).
#' @examples
#' x <- 45:50
#' qmat <- cbind(
#'   q1 = c(.011, .012, .013, .014, .015, .016),
#'   q2 = c(.100, .100, .100, .100, .100, .100)
#' )
#' md_table(x, qmat, radix = 1000)
#' @export
md_table <- function(x, qxj, radix = 100000) {
  x <- as.integer(x)
  radix <- .check_nonneg_scalar(radix, "radix")

  qxj <- as.matrix(qxj)
  storage.mode(qxj) <- "numeric"
  if (nrow(qxj) != length(x)) {
    stop("nrow(qxj) must equal length(x).", call. = FALSE)
  }
  if (any(!is.finite(qxj)) || any(qxj < 0) || any(qxj > 1)) {
    stop("qxj entries must lie in [0, 1].", call. = FALSE)
  }

  qtau <- rowSums(qxj)
  if (any(qtau > 1)) stop("Each row sum of qxj must not exceed 1.", call. = FALSE)
  ptau <- 1 - qtau

  n <- length(x)
  lx <- numeric(n)
  lx[1] <- radix
  if (n >= 2) {
    for (i in 2:n) {
      lx[i] <- lx[i - 1] * ptau[i - 1]
    }
  }

  dx_mat <- lx * qxj
  dxt <- lx * qtau

  out <- data.frame(x = x, stringsAsFactors = FALSE)

  for (j in seq_len(ncol(qxj))) {
    out[[paste0("q", j)]] <- qxj[, j]
  }
  out$qtau <- qtau
  out$ptau <- ptau
  out$ltau <- lx
  for (j in seq_len(ncol(qxj))) {
    out[[paste0("d", j)]] <- dx_mat[, j]
  }
  out$dtau <- dxt
  out
}

#' n-year survival probability {}_n p_x^(tau) from a multiple-decrement table
#'
#' @param tbl Output from md_table().
#' @param x Starting age.
#' @param n Number of years.
#' @return Numeric scalar.
#' @examples
#' x <- 45:50
#' qmat <- cbind(q1 = c(.011, .012, .013, .014, .015, .016), q2 = rep(.1, 6))
#' tbl <- md_table(x, qmat, radix = 1000)
#' npxtau_md(tbl, x = 46, n = 3)
#' @export
npxtau_md <- function(tbl, x, n) {
  n <- as.integer(n)
  if (n < 0) stop("n must be nonnegative.", call. = FALSE)
  idx0 <- match(x, tbl$x)
  idx1 <- match(x + n, tbl$x)
  if (is.na(idx0) || is.na(idx1)) stop("x and x+n must be in tbl$x.", call. = FALSE)
  tbl$ltau[idx1] / tbl$ltau[idx0]
}

#' n-year probability of decrement due to cause j from a multiple-decrement table
#'
#' @param tbl Output from md_table().
#' @param x Starting age.
#' @param n Number of years.
#' @param j Cause index.
#' @return Numeric scalar.
#' @examples
#' x <- 45:50
#' qmat <- cbind(q1 = c(.011, .012, .013, .014, .015, .016), q2 = rep(.1, 6))
#' tbl <- md_table(x, qmat, radix = 1000)
#' nqxj_md(tbl, x = 46, n = 2, j = 1)
#' @export
nqxj_md <- function(tbl, x, n, j) {
  n <- as.integer(n)
  j <- as.integer(j)
  if (n < 0) stop("n must be nonnegative.", call. = FALSE)
  idx0 <- match(x, tbl$x)
  if (is.na(idx0)) stop("x must be in tbl$x.", call. = FALSE)

  rows <- which(tbl$x >= x & tbl$x < x + n)
  dcol <- paste0("d", j)
  if (!dcol %in% names(tbl)) stop("Cause j not found in table.", call. = FALSE)

  sum(tbl[[dcol]][rows]) / tbl$ltau[idx0]
}

#' n-year probability of decrement for all causes from a multiple-decrement table
#'
#' @param tbl Output from md_table().
#' @param x Starting age.
#' @param n Number of years.
#' @return Numeric scalar.
#' @examples
#' x <- 45:50
#' qmat <- cbind(q1 = c(.011, .012, .013, .014, .015, .016), q2 = rep(.1, 6))
#' tbl <- md_table(x, qmat, radix = 1000)
#' nqxtau_md(tbl, x = 46, n = 2)
#' @export
nqxtau_md <- function(tbl, x, n) {
  1 - npxtau_md(tbl = tbl, x = x, n = n)
}

# ---------------------------
# Continuous multiple-decrement functions under constant force
# ---------------------------

#' Single-decrement survival probability {}_t p_x^{'(j)} under constant force
#'
#' @param mu Force of decrement for cause j.
#' @param t Time.
#' @return Numeric scalar/vector.
#' @examples
#' tpxprimej_cf(0.10, 5)
#' @export
tpxprimej_cf <- function(mu, t) {
  mu <- .check_nonneg_scalar(mu, "mu")
  t <- as.numeric(t)
  if (any(!is.finite(t)) || any(t < 0)) stop("t must be nonnegative.", call. = FALSE)
  exp(-mu * t)
}

#' Single-decrement failure probability {}_t q_x^{'(j)} under constant force
#'
#' @param mu Force of decrement for cause j.
#' @param t Time.
#' @return Numeric scalar/vector.
#' @examples
#' tqxprimej_cf(0.10, 5)
#' @export
tqxprimej_cf <- function(mu, t) {
  1 - tpxprimej_cf(mu = mu, t = t)
}

#' Total survival probability {}_t p_x^(tau) under constant forces
#'
#' @param mu Numeric vector of cause-specific forces.
#' @param t Time.
#' @return Numeric scalar/vector.
#' @examples
#' tpx_tau_cf(c(0.10, 0.20), 5)
#' @export
tpx_tau_cf <- function(mu, t) {
  mu <- as.numeric(mu)
  if (any(!is.finite(mu)) || any(mu < 0)) stop("mu must be nonnegative.", call. = FALSE)
  tpxprimej_cf(mu = sum(mu), t = t)
}

#' Cause-specific probability {}_t q_x^(j) under constant forces
#'
#' @param mu Numeric vector of cause-specific forces.
#' @param j Cause index.
#' @param t Time.
#' @return Numeric scalar/vector.
#' @examples
#' tqxj_cf(c(0.10, 0.20), j = 1, t = 5)
#' @export
tqxj_cf <- function(mu, j, t) {
  mu <- as.numeric(mu)
  if (any(!is.finite(mu)) || any(mu < 0)) stop("mu must be nonnegative.", call. = FALSE)
  j <- as.integer(j)
  if (j < 1 || j > length(mu)) stop("Invalid cause index j.", call. = FALSE)
  mu_tau <- sum(mu)
  if (mu_tau == 0) return(rep(0, length(t)))
  (mu[j] / mu_tau) * (1 - exp(-mu_tau * t))
}

# ---------------------------
# MUDD relationships
# ---------------------------

#' Independent probabilities q_x^{'(j)} from dependent probabilities q_x^(j) under MUDD
#'
#' @param qxj Numeric vector of dependent probabilities q_x^(j).
#' @return Numeric vector of independent probabilities q_x^{'(j)}.
#' @examples
#' qxprime_mudd(c(.20, .10))
#' @export
qxprime_mudd <- function(qxj) {
  qxj <- .check_prob_vec(qxj, "qxj")
  qtau <- sum(qxj)
  if (qtau <= 0 || qtau >= 1) stop("sum(qxj) must lie in (0, 1).", call. = FALSE)
  1 - (1 - qtau)^(qxj / qtau)
}

#' Independent probabilities {}_t q_x^{'(j)} under MUDD
#'
#' @param qxj Numeric vector of dependent probabilities q_x^(j).
#' @param t Time in [0,1].
#' @return Numeric vector.
#' @examples
#' tqxprime_mudd(c(.20, .10), t = 0.5)
#' @export
tqxprime_mudd <- function(qxj, t) {
  qxj <- .check_prob_vec(qxj, "qxj")
  t <- as.numeric(t)
  if (any(!is.finite(t)) || any(t < 0) || any(t > 1)) {
    stop("t must lie in [0, 1].", call. = FALSE)
  }
  qtau <- sum(qxj)
  if (qtau <= 0 || qtau >= 1) stop("sum(qxj) must lie in (0, 1).", call. = FALSE)
  1 - (1 - t * qtau)^(qxj / qtau)
}

# ---------------------------
# Constant-force relationships from independent q'
# ---------------------------

#' Dependent probabilities q_x^(j) from independent probabilities q_x^{'(j)} under constant force
#'
#' @param qxprime Numeric vector of independent probabilities.
#' @return Numeric vector of dependent probabilities.
#' @examples
#' qx_dep_cf(c(0.20, 0.10))
#' @export
qx_dep_cf <- function(qxprime) {
  qxprime <- .check_prob_vec(qxprime, "qxprime")
  mu <- -log(1 - qxprime)
  mu_tau <- sum(mu)
  qtau <- 1 - exp(-mu_tau)
  if (mu_tau == 0) return(rep(0, length(qxprime)))
  (mu / mu_tau) * qtau
}

# ---------------------------
# SUDD relationships for two decrements
# ---------------------------

#' Dependent probabilities q_x^(j) from independent probabilities q_x^{'(j)} under SUDD
#'
#' Two-decrement case only.
#'
#' @param q1prime Independent probability for decrement 1.
#' @param q2prime Independent probability for decrement 2.
#' @return Numeric vector c(q1, q2).
#' @examples
#' qx_dep_sudd(0.20, 0.10)
#' @export
qx_dep_sudd <- function(q1prime, q2prime) {
  q1prime <- .check_nonneg_scalar(q1prime, "q1prime")
  q2prime <- .check_nonneg_scalar(q2prime, "q2prime")
  if (q1prime > 1 || q2prime > 1) stop("Probabilities must lie in [0,1].", call. = FALSE)

  q1 <- q1prime * (1 - 0.5 * q2prime)
  q2 <- q2prime * (1 - 0.5 * q1prime)
  c(q1 = q1, q2 = q2)
}

#' Independent probabilities q_x^{'(j)} from dependent probabilities q_x^(j) under SUDD
#'
#' Two-decrement case only.
#'
#' @param q1 Dependent probability for decrement 1.
#' @param q2 Dependent probability for decrement 2.
#' @return Numeric vector c(q1prime, q2prime).
#' @examples
#' qxprime_sudd(0.20, 0.10)
#' @export
qxprime_sudd <- function(q1, q2) {
  q1 <- .check_nonneg_scalar(q1, "q1")
  q2 <- .check_nonneg_scalar(q2, "q2")
  if (q1 > 1 || q2 > 1 || q1 + q2 > 1) {
    stop("q1 and q2 must lie in [0,1] and satisfy q1 + q2 <= 1.", call. = FALSE)
  }

  # Solve quadratic from Example 13.10 / Case F
  # z^2 - (q1 - q2 + 2) z + 2 q1 = 0
  a <- 1
  b <- -(q1 - q2 + 2)
  c0 <- 2 * q1
  disc <- b^2 - 4 * a * c0
  if (disc < 0) stop("No real SUDD solution.", call. = FALSE)

  roots <- c((-b - sqrt(disc)) / (2 * a), (-b + sqrt(disc)) / (2 * a))
  roots <- roots[roots >= 0 & roots <= 1]

  if (length(roots) == 0L) stop("No admissible SUDD solution in [0,1].", call. = FALSE)

  q1prime <- min(roots)
  q2prime <- q2 / (1 - 0.5 * q1prime)

  if (!is.finite(q2prime) || q2prime < 0 || q2prime > 1) {
    stop("Derived q2prime is outside [0,1].", call. = FALSE)
  }

  c(q1prime = q1prime, q2prime = q2prime)
}
