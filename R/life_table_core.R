# life_table_core.R
#
# Core life table utilities for Chapter 6 of Models for Quantifying Risk.
# These functions follow the chapter notation as closely as possible:
#   l_x, d_x, {}_n d_x, q_x, {}_n q_x, {}_n p_x
#
# The functions below work in the discrete tabular setting only.

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_scalar_or_vector_lengths <- function(...) {
  lens <- vapply(list(...), length, integer(1))
  nz <- lens[lens > 0]
  length(unique(nz)) <= 1L
}

.validate_life_table <- function(tbl) {
  required_cols <- c("x", "lx", "dx", "qx", "px")
  if (!all(required_cols %in% names(tbl))) {
    stop("life_table object must contain columns: x, lx, dx, qx, px.", call. = FALSE)
  }

  if (any(diff(tbl$x) <= 0)) {
    stop("x values must be strictly increasing.", call. = FALSE)
  }

  if (any(tbl$lx < 0, na.rm = TRUE)) {
    stop("lx values must be nonnegative.", call. = FALSE)
  }

  invisible(tbl)
}

.get_tbl_index <- function(tbl, x) {
  match(x, tbl$x)
}

.get_lx_at <- function(tbl, x) {
  idx <- .get_tbl_index(tbl, x)
  out <- rep(NA_real_, length(x))
  ok <- !is.na(idx)
  out[ok] <- tbl$lx[idx[ok]]
  out
}

# -------------------------------------------------------------------------
# Constructor
# -------------------------------------------------------------------------

#' Construct a life table
#'
#' Build a discrete life table from one of lx, qx, px, or S0.
#'
#' @param x Numeric vector of ages.
#' @param lx Numeric vector of l_x values.
#' @param qx Numeric vector of q_x values.
#' @param px Numeric vector of p_x values.
#' @param S0 Numeric vector of S_0(x) values.
#' @param radix Radix used when converting S0 to lx, or when building from qx/px.
#'
#' @return A data.frame with class "life_table".
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

  if (is.null(radix) || length(radix) != 1L || !is.finite(radix) || radix <= 0) {
    stop("radix must be a positive finite scalar.", call. = FALSE)
  }

  if (!all(diff(x) > 0)) {
    stop("x must be strictly increasing.", call. = FALSE)
  }

  # -----------------------------------------------------------------------
  # Case 1: lx supplied directly
  # -----------------------------------------------------------------------
  if (!is.null(lx)) {
    lx <- as.numeric(lx)

    if (length(lx) != length(x)) {
      stop("When lx is supplied, length(lx) must equal length(x).", call. = FALSE)
    }

    if (any(lx < 0, na.rm = TRUE)) {
      stop("lx must be nonnegative.", call. = FALSE)
    }

    if (is.unsorted(-lx)) {
      stop("lx must be nonincreasing with age.", call. = FALSE)
    }

    tbl_x <- x
    tbl_lx <- lx
  }

  # -----------------------------------------------------------------------
  # Case 2: S0 supplied
  # -----------------------------------------------------------------------
  if (!is.null(S0)) {
    S0 <- as.numeric(S0)

    if (length(S0) != length(x)) {
      stop("When S0 is supplied, length(S0) must equal length(x).", call. = FALSE)
    }

    if (any(S0 < 0 | S0 > 1, na.rm = TRUE)) {
      stop("S0 values must lie in [0, 1].", call. = FALSE)
    }

    if (is.unsorted(-S0)) {
      stop("S0 must be nonincreasing with age.", call. = FALSE)
    }

    tbl_x <- x
    tbl_lx <- radix * S0
  }

  # -----------------------------------------------------------------------
  # Case 3: qx supplied
  # x corresponds to ages where q_x is defined
  # Returns lx on ages c(x, max(x)+1)
  # -----------------------------------------------------------------------
  if (!is.null(qx)) {
    qx <- as.numeric(qx)

    if (length(qx) != length(x)) {
      stop("When qx is supplied, length(qx) must equal length(x).", call. = FALSE)
    }

    if (any(qx < 0 | qx > 1, na.rm = TRUE)) {
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

  # -----------------------------------------------------------------------
  # Case 4: px supplied
  # x corresponds to ages where p_x is defined
  # Returns lx on ages c(x, max(x)+1)
  # -----------------------------------------------------------------------
  if (!is.null(px)) {
    px <- as.numeric(px)

    if (length(px) != length(x)) {
      stop("When px is supplied, length(px) must equal length(x).", call. = FALSE)
    }

    if (any(px < 0 | px > 1, na.rm = TRUE)) {
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

  # Derive dx, qx, px
  n <- length(tbl_lx)

  dx_vals <- c(tbl_lx[-n] - tbl_lx[-1L], NA_real_)
  qx_vals <- c(dx_vals[-n] / tbl_lx[-n], NA_real_)
  px_vals <- c(tbl_lx[-1L] / tbl_lx[-n], NA_real_)

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

#' Extract l_x
#' @param tbl A life_table object.
#' @param x Ages.
#' @return Numeric vector of l_x values.
#' @export
lx <- function(tbl, x) {
  .validate_life_table(tbl)
  .get_lx_at(tbl, x)
}

#' Compute d_x = l_x - l_{x+1}
#' @param tbl A life_table object.
#' @param x Ages.
#' @return Numeric vector of d_x values.
#' @export
dx <- function(tbl, x) {
  .validate_life_table(tbl)
  lx_x <- .get_lx_at(tbl, x)
  lx_x1 <- .get_lx_at(tbl, x + 1)
  lx_x - lx_x1
}

#' Compute {}_n d_x = l_x - l_{x+n}
#' @param tbl A life_table object.
#' @param x Ages.
#' @param n Nonnegative integer durations.
#' @return Numeric vector of {}_n d_x values.
#' @export
ndx <- function(tbl, x, n) {
  .validate_life_table(tbl)

  if (!.check_scalar_or_vector_lengths(x, n)) {
    stop("x and n must have compatible lengths.", call. = FALSE)
  }

  len <- max(length(x), length(n))
  x <- rep_len(x, len)
  n <- rep_len(n, len)

  if (any(n < 0 | n != floor(n), na.rm = TRUE)) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }

  lx_x <- .get_lx_at(tbl, x)
  lx_xn <- .get_lx_at(tbl, x + n)
  lx_x - lx_xn
}

#' Compute q_x = d_x / l_x
#' @param tbl A life_table object.
#' @param x Ages.
#' @return Numeric vector of q_x values.
#' @export
qx_tab <- function(tbl, x) {
  .validate_life_table(tbl)
  dx(tbl, x) / lx(tbl, x)
}

#' Compute {}_n p_x = l_{x+n} / l_x
#' @param tbl A life_table object.
#' @param x Ages.
#' @param n Nonnegative integer durations.
#' @return Numeric vector of {}_n p_x values.
#' @export
npx <- function(tbl, x, n) {
  .validate_life_table(tbl)

  if (!.check_scalar_or_vector_lengths(x, n)) {
    stop("x and n must have compatible lengths.", call. = FALSE)
  }

  len <- max(length(x), length(n))
  x <- rep_len(x, len)
  n <- rep_len(n, len)

  if (any(n < 0 | n != floor(n), na.rm = TRUE)) {
    stop("n must be a nonnegative integer.", call. = FALSE)
  }

  lx_x <- .get_lx_at(tbl, x)
  lx_xn <- .get_lx_at(tbl, x + n)
  lx_xn / lx_x
}

#' Compute {}_n q_x = {}_n d_x / l_x = 1 - {}_n p_x
#' @param tbl A life_table object.
#' @param x Ages.
#' @param n Nonnegative integer durations.
#' @return Numeric vector of {}_n q_x values.
#' @export
nqx <- function(tbl, x, n) {
  .validate_life_table(tbl)
  1 - npx(tbl, x, n)
}

# -------------------------------------------------------------------------
# Print method
# -------------------------------------------------------------------------

#' @export
print.life_table <- function(x, ...) {
  cat("<life_table>\n")
  print.data.frame(unclass(x), row.names = FALSE, ...)
  invisible(x)
}
