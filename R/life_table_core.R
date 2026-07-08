# life_table_core.R
#
# Core life table utilities.
# These functions work in the discrete tabular setting.

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_scalar_or_vector_lengths <- function(...) {
  lens <- vapply(list(...), length, integer(1))
  nz <- lens[lens > 0]

  if (length(nz) == 0L) return(TRUE)

  big <- nz[nz != 1L]
  length(unique(big)) <= 1L
}

.recycle_life_table_args <- function(x, n = NULL) {
  x <- as.numeric(x)

  if (length(x) == 0L || any(!is.finite(x))) {
    stop("x must contain finite numeric values.", call. = FALSE)
  }

  if (is.null(n)) {
    return(list(x = x))
  }

  n <- as.numeric(n)

  if (length(n) == 0L || any(!is.finite(n))) {
    stop("n must contain finite numeric values.", call. = FALSE)
  }

  if (!.check_scalar_or_vector_lengths(x, n)) {
    stop("x and n must have compatible lengths.", call. = FALSE)
  }

  len <- max(length(x), length(n))

  list(
    x = rep_len(x, len),
    n = rep_len(n, len)
  )
}

.validate_life_table <- function(tbl) {
  required_cols <- c("x", "lx", "dx", "qx", "px")

  if (!all(required_cols %in% names(tbl))) {
    stop("life_table object must contain columns: x, lx, dx, qx, px.",
         call. = FALSE)
  }

  if (any(!is.finite(tbl$x))) {
    stop("x values must be finite.", call. = FALSE)
  }

  if (any(diff(tbl$x) <= 0)) {
    stop("x values must be strictly increasing.", call. = FALSE)
  }

  if (any(!is.finite(tbl$lx), na.rm = TRUE)) {
    stop("lx values must be finite.", call. = FALSE)
  }

  if (any(tbl$lx < 0, na.rm = TRUE)) {
    stop("lx values must be nonnegative.", call. = FALSE)
  }

  invisible(tbl)
}

.is_closed_life_table <- function(tbl) {
  n <- nrow(tbl)
  isTRUE(tbl$lx[n] == 0)
}

.get_tbl_index <- function(tbl, x) {
  match(x, tbl$x)
}

.get_lx_at <- function(tbl, x) {
  x <- as.numeric(x)

  idx <- .get_tbl_index(tbl, x)
  out <- rep(NA_real_, length(x))

  ok <- !is.na(idx)
  out[ok] <- tbl$lx[idx[ok]]

  beyond <- is.na(idx) & x > max(tbl$x)

  if (any(beyond) && .is_closed_life_table(tbl)) {
    out[beyond] <- 0
  }

  out
}

# -------------------------------------------------------------------------
# Constructor
# -------------------------------------------------------------------------

#' Construct a life table
#'
#' Build a discrete life table from one of \code{lx}, \code{qx}, \code{px},
#' or \code{S0}.
#'
#' @param x Numeric vector of ages.
#' @param lx Numeric vector of \eqn{l_x} values.
#' @param qx Numeric vector of \eqn{q_x} values.
#' @param px Numeric vector of \eqn{p_x} values.
#' @param S0 Numeric vector of \eqn{S_0(x)} values.
#' @param radix Radix used when converting \code{S0} to \code{lx}, or when
#'   building from \code{qx} or \code{px}.
#'
#' @return A data frame with class \code{"life_table"}.
#' @export
life_table <- function(x,
                       lx = NULL,
                       qx = NULL,
                       px = NULL,
                       S0 = NULL,
                       radix = 100000) {
  supplied <- c(
    !is.null(lx),
    !is.null(qx),
    !is.null(px),
    !is.null(S0)
  )

  if (sum(supplied) != 1L) {
    stop("Supply exactly one of lx, qx, px, or S0.", call. = FALSE)
  }

  x <- as.numeric(x)

  if (length(x) == 0L || any(!is.finite(x))) {
    stop("x must contain finite numeric values.", call. = FALSE)
  }

  if (is.null(radix) || length(radix) != 1L ||
      !is.finite(radix) || radix <= 0) {
    stop("radix must be a positive finite scalar.", call. = FALSE)
  }

  if (!all(diff(x) > 0)) {
    stop("x must be strictly increasing.", call. = FALSE)
  }

  if (!is.null(lx)) {
    lx <- as.numeric(lx)

    if (length(lx) != length(x)) {
      stop("When lx is supplied, length(lx) must equal length(x).",
           call. = FALSE)
    }

    if (any(!is.finite(lx))) {
      stop("lx must contain finite numeric values.", call. = FALSE)
    }

    if (any(lx < 0)) {
      stop("lx must be nonnegative.", call. = FALSE)
    }

    if (is.unsorted(-lx)) {
      stop("lx must be nonincreasing with age.", call. = FALSE)
    }

    tbl_x <- x
    tbl_lx <- lx
  }

  if (!is.null(S0)) {
    S0 <- as.numeric(S0)

    if (length(S0) != length(x)) {
      stop("When S0 is supplied, length(S0) must equal length(x).",
           call. = FALSE)
    }

    if (any(!is.finite(S0))) {
      stop("S0 must contain finite numeric values.", call. = FALSE)
    }

    if (any(S0 < 0 | S0 > 1)) {
      stop("S0 values must lie in [0, 1].", call. = FALSE)
    }

    if (is.unsorted(-S0)) {
      stop("S0 must be nonincreasing with age.", call. = FALSE)
    }

    tbl_x <- x
    tbl_lx <- radix * S0
  }

  if (!is.null(qx)) {
    qx <- as.numeric(qx)

    if (length(qx) != length(x)) {
      stop("When qx is supplied, length(qx) must equal length(x).",
           call. = FALSE)
    }

    if (any(!is.finite(qx))) {
      stop("qx must contain finite numeric values.", call. = FALSE)
    }

    if (any(qx < 0 | qx > 1)) {
      stop("qx values must lie in [0, 1].", call. = FALSE)
    }

    lx_vals <- numeric(length(qx) + 1L)
    lx_vals[1] <- radix

    for (i in seq_along(qx)) {
      lx_vals[i + 1L] <- lx_vals[i] * (1 - qx[i])
    }

    tbl_x <- c(x, max(x) + 1)
    tbl_lx <- lx_vals
  }

  if (!is.null(px)) {
    px <- as.numeric(px)

    if (length(px) != length(x)) {
      stop("When px is supplied, length(px) must equal length(x).",
           call. = FALSE)
    }

    if (any(!is.finite(px))) {
      stop("px must contain finite numeric values.", call. = FALSE)
    }

    if (any(px < 0 | px > 1)) {
      stop("px values must lie in [0, 1].", call. = FALSE)
    }

    lx_vals <- numeric(length(px) + 1L)
    lx_vals[1] <- radix

    for (i in seq_along(px)) {
      lx_vals[i + 1L] <- lx_vals[i] * px[i]
    }

    tbl_x <- c(x, max(x) + 1)
    tbl_lx <- lx_vals
  }

  n <- length(tbl_lx)

  dx_vals <- c(tbl_lx[-n] - tbl_lx[-1L], NA_real_)

  qx_vals <- rep(NA_real_, n)
  px_vals <- rep(NA_real_, n)

  denom <- tbl_lx[-n]
  positive <- denom > 0

  qx_vals[-n][positive] <- dx_vals[-n][positive] / denom[positive]
  px_vals[-n][positive] <- tbl_lx[-1L][positive] / denom[positive]

  out <- data.frame(
    x = tbl_x,
    lx = tbl_lx,
    dx = dx_vals,
    qx = qx_vals,
    px = px_vals,
    row.names = NULL
  )

  class(out) <- c("life_table", class(out))
  .validate_life_table(out)
  out
}

# -------------------------------------------------------------------------
# Accessors and derived quantities
# -------------------------------------------------------------------------

#' Extract life-table survivor values
#'
#' @param tbl A \code{life_table} object.
#' @param x Ages.
#'
#' @return Numeric vector of \eqn{l_x} values.
#' @export
lx <- function(tbl, x) {
  .validate_life_table(tbl)
  args <- .recycle_life_table_args(x)
  .get_lx_at(tbl, args$x)
}

#' Compute deaths between ages x and x+1
#'
#' @param tbl A \code{life_table} object.
#' @param x Ages.
#'
#' @return Numeric vector of \eqn{d_x} values.
#' @export
dx <- function(tbl, x) {
  .validate_life_table(tbl)
  args <- .recycle_life_table_args(x)

  lx_x <- .get_lx_at(tbl, args$x)
  lx_x1 <- .get_lx_at(tbl, args$x + 1)

  lx_x - lx_x1
}

#' Compute deaths over an n-year interval from a life table
#'
#' @param tbl A \code{life_table} object.
#' @param x Ages.
#' @param n Nonnegative integer durations.
#'
#' @return Numeric vector of \eqn{{}_n d_x} values.
#' @export
ndx <- function(tbl, x, n) {
  .validate_life_table(tbl)

  args <- .recycle_life_table_args(x, n)
  x <- args$x
  n <- args$n

  if (any(n < 0 | n != floor(n), na.rm = TRUE)) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }

  lx_x <- .get_lx_at(tbl, x)
  lx_xn <- .get_lx_at(tbl, x + n)

  lx_x - lx_xn
}

#' Compute one-year death probability from a life table
#'
#' @param tbl A \code{life_table} object.
#' @param x Ages.
#'
#' @return Numeric vector of \eqn{q_x} values.
#' @export
qx_tab <- function(tbl, x) {
  .validate_life_table(tbl)
  args <- .recycle_life_table_args(x)

  lx_x <- lx(tbl, args$x)
  dx_x <- dx(tbl, args$x)

  out <- dx_x / lx_x
  out[!is.finite(out)] <- NA_real_

  out
}

#' Compute n-year survival probability from a life table
#'
#' @param tbl A \code{life_table} object.
#' @param x Ages.
#' @param n Nonnegative integer durations.
#'
#' @return Numeric vector of \eqn{{}_n p_x} values.
#' @export
npx <- function(tbl, x, n) {
  .validate_life_table(tbl)

  args <- .recycle_life_table_args(x, n)
  x <- args$x
  n <- args$n

  if (any(n < 0 | n != floor(n), na.rm = TRUE)) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }

  lx_x <- .get_lx_at(tbl, x)
  lx_xn <- .get_lx_at(tbl, x + n)

  out <- lx_xn / lx_x

  zero_duration <- n == 0 & !is.na(lx_x) & lx_x > 0
  out[zero_duration] <- 1

  out[!is.finite(out)] <- NA_real_

  out
}

#' Compute n-year death probability from a life table
#'
#' @param tbl A \code{life_table} object.
#' @param x Ages.
#' @param n Nonnegative integer durations.
#'
#' @return Numeric vector of \eqn{{}_n q_x} values.
#' @export
nqx <- function(tbl, x, n) {
  .validate_life_table(tbl)

  p <- npx(tbl, x, n)
  out <- 1 - p
  out[is.na(p)] <- NA_real_

  out
}

# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

#' @export
print.life_table <- function(x, ...) {
  cat("<life_table>\n")
  print.data.frame(as.data.frame(x), row.names = FALSE, ...)
  invisible(x)
}
