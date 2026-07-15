# -------------------------------------------------------------------------
# Internal validation and recycling helpers
# -------------------------------------------------------------------------

#' @noRd
.decrement_check_probability <- function(x, name = "x") {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop(name, " must contain finite values.", call. = FALSE)
  }

  if (any(x < 0 | x > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }

  x
}

#' @noRd
.decrement_check_nonnegative <- function(x, name = "x") {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  if (any(!is.finite(x)) || any(x < 0)) {
    stop(
      name,
      " must contain nonnegative finite values.",
      call. = FALSE
    )
  }

  x
}

#' @noRd
.decrement_check_nonnegative_scalar <- function(x, name = "x") {
  x <- .decrement_check_nonnegative(x, name)

  if (length(x) != 1L) {
    stop(
      name,
      " must be a single nonnegative finite number.",
      call. = FALSE
    )
  }

  x
}

#' @noRd
.decrement_check_integerish <- function(x, name = "x", positive = FALSE) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  invalid_sign <- if (positive) {
    x <= 0
  } else {
    x < 0
  }

  if (any(!is.finite(x)) ||
      any(invalid_sign) ||
      any(abs(x - round(x)) > 1e-10)) {
    descriptor <- if (positive) "positive" else "nonnegative"

    stop(
      name,
      " must contain ",
      descriptor,
      " integer-like values.",
      call. = FALSE
    )
  }

  as.integer(round(x))
}

#' @noRd
.decrement_recycle <- function(..., .names = NULL) {
  values <- list(...)

  if (length(values) == 0L) {
    return(values)
  }

  lengths <- vapply(values, length, integer(1))

  if (any(lengths == 0L)) {
    stop("Inputs must have positive length.", call. = FALSE)
  }

  common_length <- max(lengths)

  if (is.null(.names)) {
    .names <- paste0("Argument ", seq_along(values))
  }

  if (length(.names) != length(values)) {
    stop(
      ".names must have the same length as the supplied arguments.",
      call. = FALSE
    )
  }

  for (k in seq_along(values)) {
    if (!lengths[k] %in% c(1L, common_length)) {
      stop(
        .names[k],
        " must have length 1 or the common length ",
        common_length,
        ".",
        call. = FALSE
      )
    }

    values[[k]] <- rep(values[[k]], length.out = common_length)
  }

  values
}

#' @noRd
.decrement_scalar_or_matrix <- function(rows, column_names) {
  out <- do.call(rbind, rows)
  colnames(out) <- column_names
  rownames(out) <- NULL

  if (nrow(out) == 1L) {
    return(stats::setNames(as.numeric(out[1L, ]), column_names))
  }

  out
}

#' @noRd
.decrement_validate_table <- function(tbl) {
  if (!is.data.frame(tbl)) {
    stop(
      "tbl must be a data frame produced by md_table().",
      call. = FALSE
    )
  }

  required <- c("x", "qtau", "ptau", "ltau", "dtau")
  missing_columns <- setdiff(required, names(tbl))

  if (length(missing_columns) > 0L) {
    stop(
      "tbl is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (nrow(tbl) == 0L) {
    stop("tbl must contain at least one row.", call. = FALSE)
  }

  x <- as.numeric(tbl$x)

  if (any(!is.finite(x)) ||
      any(abs(x - round(x)) > 1e-10)) {
    stop(
      "tbl$x must contain finite integer-like values.",
      call. = FALSE
    )
  }

  if (anyDuplicated(x)) {
    stop("tbl$x must not contain duplicate values.", call. = FALSE)
  }

  if (length(x) > 1L && any(diff(x) != 1)) {
    stop(
      "tbl$x must contain consecutive increasing ages or durations.",
      call. = FALSE
    )
  }

  q_columns <- grep("^q[0-9]+$", names(tbl), value = TRUE)
  d_columns <- grep("^d[0-9]+$", names(tbl), value = TRUE)

  if (length(q_columns) == 0L) {
    stop(
      "tbl must contain at least one cause-specific q column.",
      call. = FALSE
    )
  }

  if (length(d_columns) != length(q_columns)) {
    stop(
      "tbl must contain matching cause-specific q and d columns.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @noRd
.decrement_table_rows <- function(tbl, x, n) {
  start <- match(x, tbl$x)

  if (is.na(start)) {
    stop("x must be present in tbl$x.", call. = FALSE)
  }

  if (n == 0L) {
    return(integer(0))
  }

  end_age <- x + n - 1L
  end <- match(end_age, tbl$x)

  if (is.na(end)) {
    stop(
      "The interval from x through x + n - 1 must be contained in tbl$x.",
      call. = FALSE
    )
  }

  start:end
}

# -------------------------------------------------------------------------
# Discrete multiple-decrement quantities
# -------------------------------------------------------------------------

#' Total one-year decrement probability
#'
#' Computes the total one-year multiple-decrement probability
#' \deqn{q_x^{(\tau)}=\sum_j q_x^{(j)}.}
#'
#' @param qxj Numeric vector of cause-specific decrement probabilities.
#'
#' @return A numeric scalar.
#'
#' @examples
#' qxtau(c(0.011, 0.100))
#'
#' @export
qxtau <- function(qxj) {
  qxj <- .decrement_check_probability(qxj, "qxj")

  total <- sum(qxj)

  if (total > 1 + sqrt(.Machine$double.eps)) {
    stop("sum(qxj) cannot exceed 1.", call. = FALSE)
  }

  min(total, 1)
}

#' Total one-year survival probability
#'
#' Computes
#' \deqn{p_x^{(\tau)}=1-q_x^{(\tau)}.}
#'
#' @inheritParams qxtau
#'
#' @return A numeric scalar.
#'
#' @examples
#' pxtau(c(0.011, 0.100))
#'
#' @export
pxtau <- function(qxj) {
  1 - qxtau(qxj)
}

#' Cause-specific numbers of decrements
#'
#' Computes
#' \deqn{d_x^{(j)}=l_x^{(\tau)}q_x^{(j)}.}
#'
#' @param lxtau Number alive at age or duration \eqn{x} in the
#'   multiple-decrement table.
#' @param qxj Numeric vector of cause-specific decrement probabilities.
#'
#' @return A numeric vector containing the number of decrements from each
#'   cause.
#'
#' @examples
#' dxj(1000, c(0.011, 0.100))
#'
#' @export
dxj <- function(lxtau, qxj) {
  lxtau <- .decrement_check_nonnegative_scalar(lxtau, "lxtau")
  qxj <- .decrement_check_probability(qxj, "qxj")

  if (sum(qxj) > 1 + sqrt(.Machine$double.eps)) {
    stop("sum(qxj) cannot exceed 1.", call. = FALSE)
  }

  lxtau * qxj
}

#' Total number of decrements
#'
#' Computes
#' \deqn{d_x^{(\tau)}=\sum_j d_x^{(j)}.}
#'
#' @inheritParams dxj
#'
#' @return A numeric scalar.
#'
#' @examples
#' dxtau(1000, c(0.011, 0.100))
#'
#' @export
dxtau <- function(lxtau, qxj) {
  sum(dxj(lxtau = lxtau, qxj = qxj))
}

#' Construct a multiple-decrement table
#'
#' Constructs a discrete multiple-decrement table from cause-specific
#' one-year decrement probabilities.
#'
#' @param x Integer vector of consecutive ages or durations.
#' @param qxj Numeric matrix or data frame of cause-specific decrement
#'   probabilities. Rows correspond to values in \code{x}, and columns
#'   correspond to causes.
#' @param radix Initial value of \eqn{l_x^{(\tau)}}.
#'
#' @return An object of classes \code{"md_table"} and \code{"data.frame"}
#'   containing age or duration, cause-specific decrement probabilities,
#'   total decrement and survival probabilities, numbers alive, and
#'   cause-specific and total decrements.
#'
#' @examples
#' ages <- 45:50
#' qmat <- cbind(
#'   withdrawal = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
#'   retirement = rep(0.100, 6)
#' )
#'
#' md_table(ages, qmat, radix = 1000)
#'
#' @export
md_table <- function(x, qxj, radix = 100000) {
  x <- .decrement_check_integerish(x, "x")
  radix <- .decrement_check_nonnegative_scalar(radix, "radix")

  if (anyDuplicated(x)) {
    stop("x must not contain duplicate values.", call. = FALSE)
  }

  if (length(x) > 1L && any(diff(x) != 1L)) {
    stop(
      "x must contain consecutive increasing ages or durations.",
      call. = FALSE
    )
  }

  qxj <- as.matrix(qxj)

  if (length(qxj) == 0L || ncol(qxj) == 0L) {
    stop("qxj must contain at least one cause column.", call. = FALSE)
  }

  suppressWarnings(storage.mode(qxj) <- "double")

  if (nrow(qxj) != length(x)) {
    stop("nrow(qxj) must equal length(x).", call. = FALSE)
  }

  if (any(!is.finite(qxj)) ||
      any(qxj < 0) ||
      any(qxj > 1)) {
    stop("qxj entries must lie in [0, 1].", call. = FALSE)
  }

  qtau <- rowSums(qxj)

  if (any(qtau > 1 + sqrt(.Machine$double.eps))) {
    stop(
      "Each row sum of qxj must not exceed 1.",
      call. = FALSE
    )
  }

  qtau <- pmin(qtau, 1)
  ptau <- 1 - qtau

  row_count <- length(x)
  ltau <- numeric(row_count)
  ltau[1L] <- radix

  if (row_count > 1L) {
    for (k in 2:row_count) {
      ltau[k] <- ltau[k - 1L] * ptau[k - 1L]
    }
  }

  decrement_matrix <- qxj * ltau
  dtau <- ltau * qtau

  cause_count <- ncol(qxj)
  cause_names <- colnames(qxj)

  if (is.null(cause_names) ||
      any(!nzchar(cause_names)) ||
      anyDuplicated(cause_names)) {
    cause_names <- paste0("cause", seq_len(cause_count))
  }

  out <- data.frame(
    x = x,
    stringsAsFactors = FALSE
  )

  for (k in seq_len(cause_count)) {
    out[[paste0("q", k)]] <- qxj[, k]
  }

  out$qtau <- qtau
  out$ptau <- ptau
  out$ltau <- ltau

  for (k in seq_len(cause_count)) {
    out[[paste0("d", k)]] <- decrement_matrix[, k]
  }

  out$dtau <- dtau

  attr(out, "cause_names") <- cause_names
  class(out) <- c("md_table", "data.frame")

  out
}

#' Multiple-decrement survival probability from a table
#'
#' Computes
#' \deqn{{}_np_x^{(\tau)}
#' =\prod_{k=0}^{n-1}p_{x+k}^{(\tau)}.}
#'
#' @param tbl A multiple-decrement table produced by \code{md_table()}.
#' @param x Starting age or duration. May be scalar or vector.
#' @param n Nonnegative integer term. May be scalar or vector.
#'
#' @return A numeric vector of survival probabilities.
#'
#' @examples
#' ages <- 45:50
#' qmat <- cbind(
#'   q1 = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
#'   q2 = rep(0.100, 6)
#' )
#' tbl <- md_table(ages, qmat, radix = 1000)
#'
#' npxtau_md(tbl, x = 46, n = 3)
#'
#' @export
npxtau_md <- function(tbl, x, n) {
  .decrement_validate_table(tbl)

  x <- .decrement_check_integerish(x, "x")
  n <- .decrement_check_integerish(n, "n")

  values <- .decrement_recycle(
    x,
    n,
    .names = c("x", "n")
  )

  x <- values[[1]]
  n <- values[[2]]

  vapply(seq_along(x), function(k) {
    if (n[k] == 0L) {
      return(1)
    }

    rows <- .decrement_table_rows(
      tbl = tbl,
      x = x[k],
      n = n[k]
    )

    prod(tbl$ptau[rows])
  }, numeric(1))
}

#' Cause-specific multiple-decrement probability from a table
#'
#' Computes the probability of decrement from cause \code{j} within
#' \code{n} years:
#' \deqn{
#' {}_nq_x^{(j)}
#' =
#' \sum_{k=0}^{n-1}
#' {}_kp_x^{(\tau)}q_{x+k}^{(j)}.
#' }
#'
#' @inheritParams npxtau_md
#' @param j Positive integer cause index. May be scalar or vector.
#'
#' @return A numeric vector of cause-specific decrement probabilities.
#'
#' @examples
#' ages <- 45:50
#' qmat <- cbind(
#'   q1 = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
#'   q2 = rep(0.100, 6)
#' )
#' tbl <- md_table(ages, qmat, radix = 1000)
#'
#' nqxj_md(tbl, x = 46, n = 2, j = 1)
#'
#' @export
nqxj_md <- function(tbl, x, n, j) {
  .decrement_validate_table(tbl)

  x <- .decrement_check_integerish(x, "x")
  n <- .decrement_check_integerish(n, "n")
  j <- .decrement_check_integerish(j, "j", positive = TRUE)

  values <- .decrement_recycle(
    x,
    n,
    j,
    .names = c("x", "n", "j")
  )

  x <- values[[1]]
  n <- values[[2]]
  j <- values[[3]]

  cause_count <- length(grep("^q[0-9]+$", names(tbl)))

  if (any(j > cause_count)) {
    stop(
      "j must identify an available cause in tbl.",
      call. = FALSE
    )
  }

  vapply(seq_along(x), function(k) {
    if (n[k] == 0L) {
      return(0)
    }

    rows <- .decrement_table_rows(
      tbl = tbl,
      x = x[k],
      n = n[k]
    )

    cause_column <- paste0("q", j[k])
    survival_to_start <- c(
      1,
      cumprod(tbl$ptau[rows])[-length(rows)]
    )

    sum(survival_to_start * tbl[[cause_column]][rows])
  }, numeric(1))
}

#' Total multiple-decrement probability from a table
#'
#' Computes
#' \deqn{{}_nq_x^{(\tau)}=1-{}_np_x^{(\tau)}.}
#'
#' @inheritParams npxtau_md
#'
#' @return A numeric vector of total decrement probabilities.
#'
#' @examples
#' ages <- 45:50
#' qmat <- cbind(
#'   q1 = c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016),
#'   q2 = rep(0.100, 6)
#' )
#' tbl <- md_table(ages, qmat, radix = 1000)
#'
#' nqxtau_md(tbl, x = 46, n = 2)
#'
#' @export
nqxtau_md <- function(tbl, x, n) {
  1 - npxtau_md(tbl = tbl, x = x, n = n)
}

# -------------------------------------------------------------------------
# Constant-force multiple-decrement quantities
# -------------------------------------------------------------------------

#' Single-decrement survival under a constant force
#'
#' Computes
#' \deqn{{}_tp_x^{\prime(j)}=\exp(-\mu_jt).}
#'
#' @param mu Nonnegative constant force of decrement. May be scalar or
#'   vector.
#' @param t Nonnegative time. May be scalar or vector.
#'
#' @return A numeric vector.
#'
#' @examples
#' tpxprimej_cf(0.10, 5)
#' tpxprimej_cf(0.10, c(1, 5, 10))
#'
#' @export
tpxprimej_cf <- function(mu, t) {
  mu <- .decrement_check_nonnegative(mu, "mu")
  t <- .decrement_check_nonnegative(t, "t")

  values <- .decrement_recycle(
    mu,
    t,
    .names = c("mu", "t")
  )

  mu <- values[[1]]
  t <- values[[2]]

  exp(-mu * t)
}

#' Single-decrement failure under a constant force
#'
#' Computes
#' \deqn{{}_tq_x^{\prime(j)}
#' =1-\exp(-\mu_jt).}
#'
#' @inheritParams tpxprimej_cf
#'
#' @return A numeric vector.
#'
#' @examples
#' tqxprimej_cf(0.10, 5)
#'
#' @export
tqxprimej_cf <- function(mu, t) {
  1 - tpxprimej_cf(mu = mu, t = t)
}

#' Total survival under constant cause-specific forces
#'
#' Computes
#' \deqn{
#' {}_tp_x^{(\tau)}
#' =
#' \exp\left(-t\sum_j\mu_j\right).
#' }
#'
#' @param mu Numeric vector of nonnegative cause-specific forces.
#' @param t Nonnegative time. May be scalar or vector.
#'
#' @return A numeric vector.
#'
#' @examples
#' tpx_tau_cf(c(0.10, 0.20), 5)
#' tpx_tau_cf(c(0.10, 0.20), c(1, 5, 10))
#'
#' @export
tpx_tau_cf <- function(mu, t) {
  mu <- .decrement_check_nonnegative(mu, "mu")
  t <- .decrement_check_nonnegative(t, "t")

  exp(-sum(mu) * t)
}

#' Cause-specific decrement probability under constant forces
#'
#' Computes
#' \deqn{
#' {}_tq_x^{(j)}
#' =
#' \frac{\mu_j}{\sum_k\mu_k}
#' \left[
#' 1-\exp\left(-t\sum_k\mu_k\right)
#' \right].
#' }
#'
#' @param mu Numeric vector of nonnegative cause-specific forces.
#' @param j Positive integer cause index.
#' @param t Nonnegative time. May be scalar or vector.
#'
#' @return A numeric vector.
#'
#' @examples
#' tqxj_cf(c(0.10, 0.20), j = 1, t = 5)
#'
#' @export
tqxj_cf <- function(mu, j, t) {
  mu <- .decrement_check_nonnegative(mu, "mu")
  j <- .decrement_check_integerish(j, "j", positive = TRUE)
  t <- .decrement_check_nonnegative(t, "t")

  if (length(j) != 1L) {
    stop("j must be a single cause index.", call. = FALSE)
  }

  if (j > length(mu)) {
    stop("j must identify an available cause in mu.", call. = FALSE)
  }

  total_force <- sum(mu)

  if (total_force == 0) {
    return(rep(0, length(t)))
  }

  (mu[j] / total_force) *
    (1 - exp(-total_force * t))
}

# -------------------------------------------------------------------------
# Multiple-decrement uniform distribution relationships
# -------------------------------------------------------------------------

#' Associated single-decrement probabilities under MUDD
#'
#' Converts dependent multiple-decrement probabilities to associated
#' single-decrement probabilities under the multiple-decrement uniform
#' distribution assumption.
#'
#' @param qxj Numeric vector of dependent cause-specific decrement
#'   probabilities.
#'
#' @return A numeric vector of associated single-decrement probabilities.
#'
#' @examples
#' qxprime_mudd(c(0.20, 0.10))
#'
#' @export
qxprime_mudd <- function(qxj) {
  qxj <- .decrement_check_probability(qxj, "qxj")
  qtau <- sum(qxj)

  if (qtau > 1 + sqrt(.Machine$double.eps)) {
    stop("sum(qxj) cannot exceed 1.", call. = FALSE)
  }

  if (qtau == 0) {
    return(rep(0, length(qxj)))
  }

  if (qtau >= 1) {
    stop(
      "sum(qxj) must be less than 1 for the MUDD conversion.",
      call. = FALSE
    )
  }

  1 - (1 - qtau)^(qxj / qtau)
}

#' Fractional-year associated single-decrement probabilities under MUDD
#'
#' Computes the associated single-decrement probabilities through time
#' \code{t} under the multiple-decrement uniform distribution assumption.
#'
#' @param qxj Numeric vector of dependent cause-specific decrement
#'   probabilities.
#' @param t A single time in \eqn{[0,1]}.
#'
#' @return A numeric vector of associated single-decrement probabilities.
#'
#' @examples
#' tqxprime_mudd(c(0.20, 0.10), t = 0.5)
#'
#' @export
tqxprime_mudd <- function(qxj, t) {
  qxj <- .decrement_check_probability(qxj, "qxj")
  t <- .decrement_check_probability(t, "t")

  if (length(t) != 1L) {
    stop("t must be a single value in [0, 1].", call. = FALSE)
  }

  qtau <- sum(qxj)

  if (qtau > 1 + sqrt(.Machine$double.eps)) {
    stop("sum(qxj) cannot exceed 1.", call. = FALSE)
  }

  if (qtau == 0 || t == 0) {
    return(rep(0, length(qxj)))
  }

  if (qtau >= 1 && t >= 1) {
    stop(
      "t * sum(qxj) must be less than 1 for the MUDD conversion.",
      call. = FALSE
    )
  }

  1 - (1 - t * qtau)^(qxj / qtau)
}

# -------------------------------------------------------------------------
# Constant-force conversion from associated single-decrement probabilities
# -------------------------------------------------------------------------

#' Multiple-decrement probabilities under constant forces
#'
#' Converts associated single-decrement probabilities to dependent
#' cause-specific decrement probabilities under constant forces.
#'
#' @param qxprime Numeric vector of associated single-decrement
#'   probabilities. Each value must be less than one.
#'
#' @return A numeric vector of dependent cause-specific decrement
#'   probabilities.
#'
#' @examples
#' qx_dep_cf(c(0.20, 0.10))
#'
#' @export
qx_dep_cf <- function(qxprime) {
  qxprime <- .decrement_check_probability(qxprime, "qxprime")

  if (any(qxprime >= 1)) {
    stop(
      "qxprime values must be less than 1 under a finite constant-force model.",
      call. = FALSE
    )
  }

  mu <- -log1p(-qxprime)
  total_force <- sum(mu)

  if (total_force == 0) {
    return(rep(0, length(qxprime)))
  }

  total_probability <- -expm1(-total_force)

  (mu / total_force) * total_probability
}

# -------------------------------------------------------------------------
# Single-decrement uniform distribution relationships
# -------------------------------------------------------------------------

#' Multiple-decrement probabilities under SUDD
#'
#' Converts two associated single-decrement probabilities to dependent
#' multiple-decrement probabilities under the single-decrement uniform
#' distribution assumption.
#'
#' @param q1prime Associated single-decrement probability for cause 1.
#'   May be scalar or vector.
#' @param q2prime Associated single-decrement probability for cause 2.
#'   May be scalar or vector.
#'
#' @return For scalar input, a named numeric vector containing \code{q1}
#'   and \code{q2}. For vectorized input, a numeric matrix with columns
#'   \code{q1} and \code{q2}.
#'
#' @examples
#' qx_dep_sudd(0.20, 0.10)
#'
#' @export
qx_dep_sudd <- function(q1prime, q2prime) {
  q1prime <- .decrement_check_probability(q1prime, "q1prime")
  q2prime <- .decrement_check_probability(q2prime, "q2prime")

  values <- .decrement_recycle(
    q1prime,
    q2prime,
    .names = c("q1prime", "q2prime")
  )

  q1prime <- values[[1]]
  q2prime <- values[[2]]

  rows <- lapply(seq_along(q1prime), function(k) {
    c(
      q1 = q1prime[k] * (1 - 0.5 * q2prime[k]),
      q2 = q2prime[k] * (1 - 0.5 * q1prime[k])
    )
  })

  .decrement_scalar_or_matrix(
    rows,
    column_names = c("q1", "q2")
  )
}

#' Associated single-decrement probabilities under SUDD
#'
#' Converts two dependent multiple-decrement probabilities to associated
#' single-decrement probabilities under the single-decrement uniform
#' distribution assumption.
#'
#' @param q1 Dependent decrement probability for cause 1. May be scalar or
#'   vector.
#' @param q2 Dependent decrement probability for cause 2. May be scalar or
#'   vector.
#'
#' @return For scalar input, a named numeric vector containing
#'   \code{q1prime} and \code{q2prime}. For vectorized input, a numeric
#'   matrix with columns \code{q1prime} and \code{q2prime}.
#'
#' @examples
#' qxprime_sudd(0.20, 0.10)
#'
#' @export
qxprime_sudd <- function(q1, q2) {
  q1 <- .decrement_check_probability(q1, "q1")
  q2 <- .decrement_check_probability(q2, "q2")

  values <- .decrement_recycle(
    q1,
    q2,
    .names = c("q1", "q2")
  )

  q1 <- values[[1]]
  q2 <- values[[2]]

  if (any(q1 + q2 > 1 + sqrt(.Machine$double.eps))) {
    stop(
      "q1 and q2 must satisfy q1 + q2 <= 1.",
      call. = FALSE
    )
  }

  rows <- lapply(seq_along(q1), function(k) {
    if (q1[k] == 0 && q2[k] == 0) {
      return(c(q1prime = 0, q2prime = 0))
    }

    coefficient_b <- -(q1[k] - q2[k] + 2)
    coefficient_c <- 2 * q1[k]
    discriminant <- coefficient_b^2 - 4 * coefficient_c

    if (discriminant < -sqrt(.Machine$double.eps)) {
      stop(
        "No real SUDD solution exists for one or more input pairs.",
        call. = FALSE
      )
    }

    discriminant <- max(discriminant, 0)

    roots <- c(
      (-coefficient_b - sqrt(discriminant)) / 2,
      (-coefficient_b + sqrt(discriminant)) / 2
    )

    roots <- roots[
      roots >= -sqrt(.Machine$double.eps) &
        roots <= 1 + sqrt(.Machine$double.eps)
    ]

    if (length(roots) == 0L) {
      stop(
        "No admissible SUDD solution exists in [0, 1].",
        call. = FALSE
      )
    }

    roots <- pmin(pmax(roots, 0), 1)

    candidate_rows <- lapply(roots, function(q1prime) {
      denominator <- 1 - 0.5 * q1prime

      if (denominator <= 0) {
        return(NULL)
      }

      q2prime <- q2[k] / denominator

      if (!is.finite(q2prime) ||
          q2prime < -sqrt(.Machine$double.eps) ||
          q2prime > 1 + sqrt(.Machine$double.eps)) {
        return(NULL)
      }

      q2prime <- min(max(q2prime, 0), 1)

      reconstructed <- qx_dep_sudd(q1prime, q2prime)

      error <- max(abs(
        reconstructed -
          c(q1 = q1[k], q2 = q2[k])
      ))

      list(
        q1prime = q1prime,
        q2prime = q2prime,
        error = error
      )
    })

    candidate_rows <- Filter(Negate(is.null), candidate_rows)

    if (length(candidate_rows) == 0L) {
      stop(
        "No admissible SUDD solution exists in [0, 1].",
        call. = FALSE
      )
    }

    errors <- vapply(
      candidate_rows,
      function(candidate) candidate$error,
      numeric(1)
    )

    best <- candidate_rows[[which.min(errors)]]

    c(
      q1prime = best$q1prime,
      q2prime = best$q2prime
    )
  })

  .decrement_scalar_or_matrix(
    rows,
    column_names = c("q1prime", "q2prime")
  )
}
