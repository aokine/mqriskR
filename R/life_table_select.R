# life_table_select.R
#
# Select and ultimate life table utilities for Chapter 6 of
# Models for Quantifying Risk.
#
# The functions here support select table quantities such as:
#   l_[x]+t
#   {}_n p_[x]+t
#   {}_n q_[x]+t
#   {}_{n|m} q_[x]+t
#
# We assume the select table is supplied in "long" form with columns:
#   x_sel         age at selection
#   duration      duration since selection
#   attained_age  x_sel + duration
#   lx            select-table survivor value
#
# For attained ages beyond the select period, the table should contain
# rows with duration == select_period and lx equal to the corresponding
# ultimate l_x value at attained age x_sel + select_period.
#
# Example interpretation (4-year select period):
#   x_sel = 20, duration = 0,1,2,3,4
#   attained_age = 20,21,22,23,24
#   lx = l_[20], l_[20]+1, l_[20]+2, l_[20]+3, l_24

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.validate_select_life_table <- function(tbl) {
  required_cols <- c("x_sel", "duration", "attained_age", "lx")
  if (!all(required_cols %in% names(tbl))) {
    stop(
      "select_life_table object must contain columns: x_sel, duration, attained_age, lx.",
      call. = FALSE
    )
  }

  if (any(!is.finite(tbl$x_sel)) ||
      any(!is.finite(tbl$duration)) ||
      any(!is.finite(tbl$attained_age)) ||
      any(!is.finite(tbl$lx))) {
    stop("All select table columns must be finite.", call. = FALSE)
  }

  if (any(tbl$duration < 0)) {
    stop("duration must be nonnegative.", call. = FALSE)
  }

  if (any(tbl$lx < 0)) {
    stop("lx values must be nonnegative.", call. = FALSE)
  }

  invisible(tbl)
}

.get_select_row <- function(tbl, x_sel, duration) {
  idx <- which(tbl$x_sel == x_sel & tbl$duration == duration)
  if (length(idx) == 0L) return(NA_integer_)
  idx[1L]
}

.get_select_lx <- function(tbl, x_sel, duration) {
  idx <- .get_select_row(tbl, x_sel, duration)
  if (is.na(idx)) return(NA_real_)
  tbl$lx[idx]
}

.recycle_select_args <- function(x_sel, t) {
  x_sel <- as.numeric(x_sel)
  t <- as.numeric(t)

  if (length(x_sel) == 0L || length(t) == 0L) {
    stop("x_sel and t must have positive length.", call. = FALSE)
  }

  if (length(x_sel) == 1L && length(t) > 1L) x_sel <- rep(x_sel, length(t))
  if (length(t) == 1L && length(x_sel) > 1L) t <- rep(t, length(x_sel))

  if (length(x_sel) != length(t)) {
    stop("x_sel and t must have the same length, or one must have length 1.", call. = FALSE)
  }

  list(x_sel = x_sel, t = t)
}

# -------------------------------------------------------------------------
# Constructor
# -------------------------------------------------------------------------

#' Construct a select life table
#'
#' Builds a select-life-table object from vectors of selection age,
#' duration since selection, attained age, and survivor values.
#'
#' @param x_sel Numeric vector of ages at selection.
#' @param duration Numeric vector of durations since selection.
#' @param attained_age Numeric vector of attained ages.
#' @param lx Numeric vector of select-table survivor values.
#'
#' @return A data.frame with class \code{"select_life_table"}.
#' @export
select_life_table <- function(x_sel, duration, attained_age, lx) {
  x_sel <- as.numeric(x_sel)
  duration <- as.numeric(duration)
  attained_age <- as.numeric(attained_age)
  lx <- as.numeric(lx)

  n <- length(x_sel)
  if (!all(length(duration) == n, length(attained_age) == n, length(lx) == n)) {
    stop("x_sel, duration, attained_age, and lx must have the same length.", call. = FALSE)
  }

  if (any(duration < 0)) {
    stop("duration must be nonnegative.", call. = FALSE)
  }

  if (any(attained_age != x_sel + duration)) {
    stop("attained_age must equal x_sel + duration row by row.", call. = FALSE)
  }

  if (any(lx < 0)) {
    stop("lx values must be nonnegative.", call. = FALSE)
  }

  out <- data.frame(
    x_sel = x_sel,
    duration = duration,
    attained_age = attained_age,
    lx = lx,
    row.names = NULL
  )

  out <- out[order(out$x_sel, out$duration), , drop = FALSE]
  class(out) <- c("select_life_table", class(out))
  .validate_select_life_table(out)
  out
}

# -------------------------------------------------------------------------
# Core accessor
# -------------------------------------------------------------------------

#' Extract select-table survivor value
#'
#' Returns \eqn{l_{[x]+t}} from a select life table.
#'
#' @param tbl A select_life_table object.
#' @param x_sel Numeric vector of ages at selection.
#' @param t Numeric vector of durations since selection.
#'
#' @return Numeric vector of survivor values.
#' @export
lx_select <- function(tbl, x_sel, t) {
  .validate_select_life_table(tbl)
  xt <- .recycle_select_args(x_sel, t)
  x_sel <- xt$x_sel
  t <- xt$t

  out <- numeric(length(x_sel))
  for (i in seq_along(x_sel)) {
    out[i] <- .get_select_lx(tbl, x_sel[i], t[i])
  }
  out
}

# -------------------------------------------------------------------------
# Select probabilities
# -------------------------------------------------------------------------

#' Select-life survival probability
#'
#' Computes
#' \eqn{{}_n p_{[x]+t} = l_{[x]+t+n} / l_{[x]+t}}
#'
#' in the discrete select-table setting.
#'
#' @param tbl A select_life_table object.
#' @param x_sel Numeric vector of ages at selection.
#' @param t Numeric vector of current durations since selection.
#' @param n Numeric vector of nonnegative integer future durations.
#'
#' @return Numeric vector of survival probabilities.
#' @export
npx_select <- function(tbl, x_sel, t, n) {
  .validate_select_life_table(tbl)

  xt <- .recycle_select_args(x_sel, t)
  x_sel <- xt$x_sel
  t <- xt$t
  n <- rep_len(as.numeric(n), length(x_sel))

  if (any(n < 0 | n != floor(n))) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }

  lx_now <- lx_select(tbl, x_sel, t)
  lx_future <- lx_select(tbl, x_sel, t + n)

  out <- lx_future / lx_now
  out[is.na(lx_now) | is.na(lx_future) | lx_now <= 0] <- NA_real_
  out
}

#' Select-life death probability
#'
#' Computes
#' \eqn{{}_n q_{[x]+t} = 1 - {}_n p_{[x]+t}}.
#'
#' @inheritParams npx_select
#' @return Numeric vector of death probabilities.
#' @export
nqx_select <- function(tbl, x_sel, t, n) {
  1 - npx_select(tbl, x_sel, t, n)
}

#' Deferred select-life death probability
#'
#' Computes
#' \eqn{{}_{n|m} q_{[x]+t} = {}_n p_{[x]+t} \cdot {}_m q_{[x]+t+n}}.
#'
#' @param tbl A select_life_table object.
#' @param x_sel Numeric vector of ages at selection.
#' @param t Numeric vector of current durations since selection.
#' @param n Numeric vector of nonnegative integer deferred periods.
#' @param m Numeric vector of nonnegative integer death windows.
#'
#' @return Numeric vector of deferred death probabilities.
#' @export
nmxq_select <- function(tbl, x_sel, t, n, m) {
  .validate_select_life_table(tbl)

  xt <- .recycle_select_args(x_sel, t)
  x_sel <- xt$x_sel
  t <- xt$t
  n <- rep_len(as.numeric(n), length(x_sel))
  m <- rep_len(as.numeric(m), length(x_sel))

  if (any(n < 0 | n != floor(n))) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }
  if (any(m < 0 | m != floor(m))) {
    stop("m must be a nonnegative integer.", call. = FALSE)
  }

  npx_select(tbl, x_sel, t, n) * nqx_select(tbl, x_sel, t + n, m)
}

# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

#' @export
print.select_life_table <- function(x, ...) {
  cat("<select_life_table>\n")
  print.data.frame(unclass(x), row.names = FALSE, ...)
  invisible(x)
}
