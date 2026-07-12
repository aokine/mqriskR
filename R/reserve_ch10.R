# reserve_premium_functions.R
#
# Net level premium reserve, duration-t loss, and gain-analysis functions.
#
# The exported functions in this file support parametric survival models and,
# where the underlying insurance, annuity, and premium functions permit it,
# life table objects supplied through `tbl`.
#
# Existing exported function names and the traditional argument order are
# preserved for backward compatibility.

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_nonneg_integerish <- function(x, name) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0) || any(abs(x - round(x)) > 1e-10)) {
    stop(name, " must contain nonnegative integer-like values.", call. = FALSE)
  }

  as.integer(round(x))
}

#' @noRd
.check_nonneg_numeric <- function(x, name) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop(name, " must contain nonnegative finite values.", call. = FALSE)
  }

  x
}

#' @noRd
.recycle2_ch10 <- function(a, b, a_name = "a", b_name = "b") {
  if (length(a) == 1L && length(b) > 1L) a <- rep(a, length(b))
  if (length(b) == 1L && length(a) > 1L) b <- rep(b, length(a))

  if (length(a) != length(b)) {
    stop(
      a_name, " and ", b_name,
      " must have the same length, or one must have length 1.",
      call. = FALSE
    )
  }

  list(a = a, b = b)
}

#' @noRd
.recycle3_ch10 <- function(a, b, c,
                           a_name = "a", b_name = "b", c_name = "c") {
  n <- max(length(a), length(b), length(c))

  if (!length(a) %in% c(1L, n)) {
    stop(a_name, " must have length 1 or common length.", call. = FALSE)
  }
  if (!length(b) %in% c(1L, n)) {
    stop(b_name, " must have length 1 or common length.", call. = FALSE)
  }
  if (!length(c) %in% c(1L, n)) {
    stop(c_name, " must have length 1 or common length.", call. = FALSE)
  }

  if (length(a) == 1L) a <- rep(a, n)
  if (length(b) == 1L) b <- rep(b, n)
  if (length(c) == 1L) c <- rep(c, n)

  list(a = a, b = b, c = c)
}

#' @noRd
.check_t_le_n <- function(t, n, allow_equal = TRUE) {
  if (allow_equal) {
    if (any(t > n)) stop("t must satisfy t <= n.", call. = FALSE)
  } else {
    if (any(t >= n)) stop("t must satisfy t < n.", call. = FALSE)
  }
}

#' @noRd
.check_m_scalar_ch10 <- function(m) {
  m <- as.numeric(m)

  if (length(m) != 1L || !is.finite(m) || m <= 0 || abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }

  as.integer(round(m))
}

#' @noRd
.check_reserve_source <- function(tbl, model) {
  if (!is.null(tbl) && !is.null(model)) {
    stop("Supply either tbl or model, not both.", call. = FALSE)
  }
  if (is.null(tbl) && is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.check_reserve_interest <- function(i) {
  i <- as.numeric(i)
  if (length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("i must be a single finite number greater than -1.", call. = FALSE)
  }
  i
}

# -------------------------------------------------------------------------
# Annual premium reserves
# -------------------------------------------------------------------------

#' Whole life net level premium reserve
#'
#' Computes the prospective reserve
#' for a whole life insurance with annual premiums:
#' reserve at duration t equals future APV of benefits minus future APV
#' of net premiums.
#'
#' @param x Issue age.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param model Optional parametric survival model name.
#' @param ... Additional model parameters.
#' @param tbl Optional life table object.
#'
#' @return A numeric vector of values.
#' @examples
#' tVx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVx <- function(x, t, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")

  xt <- .recycle2_ch10(x, t, "x", "t")
  x <- xt$a
  t <- xt$b

  sapply(seq_along(x), function(j) {
    Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...) -
      Px(x = x[j], i = i, tbl = tbl, model = model, ...) *
      adotx(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)
  })
}

#' Term insurance net level premium reserve
#'
#' Computes the prospective reserve for an n-year term insurance.
#'
#' @param x Issue age.
#' @param n Term in years.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param model Optional parametric survival model name.
#' @param ... Additional model parameters.
#' @param tbl Optional life table object.
#'
#' @return A numeric vector of values.
#' @examples
#' tVxn1(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVxn1 <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  .check_t_le_n(t, n, allow_equal = TRUE)

  sapply(seq_along(x), function(j) {
    if (t[j] == n[j]) return(0)

    Axn1(x = x[j] + t[j], n = n[j] - t[j], i = i, tbl = tbl, model = model, ...) -
      Pxn1(x = x[j], n = n[j], i = i, tbl = tbl, model = model, ...) *
      adotxn(x = x[j] + t[j], n = n[j] - t[j], i = i, tbl = tbl, model = model, ...)
  })
}

#' Pure endowment net level premium reserve
#'
#' Computes the prospective reserve for an n-year pure endowment.
#'
#' @inheritParams tVxn1
#' @return A numeric vector of values.
#' @examples
#' tVnEx(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVnEx <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  .check_t_le_n(t, n, allow_equal = TRUE)

  sapply(seq_along(x), function(j) {
    if (t[j] == n[j]) return(1)

    nEx(x = x[j] + t[j], n = n[j] - t[j], i = i, tbl = tbl, model = model, ...) -
      PnEx(x = x[j], n = n[j], i = i, tbl = tbl, model = model, ...) *
      adotxn(x = x[j] + t[j], n = n[j] - t[j], i = i, tbl = tbl, model = model, ...)
  })
}

#' Endowment insurance net level premium reserve
#'
#' Computes the prospective reserve for an n-year endowment insurance.
#'
#' @inheritParams tVxn1
#' @return A numeric vector of values.
#' @examples
#' tVxn(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVxn <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")

  xnt <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnt$a
  n <- xnt$b
  t <- xnt$c

  .check_t_le_n(t, n, allow_equal = TRUE)

  sapply(seq_along(x), function(j) {
    if (t[j] == n[j]) return(1)

    Axn(x = x[j] + t[j], n = n[j] - t[j], i = i, tbl = tbl, model = model, ...) -
      Pxn(x = x[j], n = n[j], i = i, tbl = tbl, model = model, ...) *
      adotxn(x = x[j] + t[j], n = n[j] - t[j], i = i, tbl = tbl, model = model, ...)
  })
}

#' h-pay whole life net level premium reserve
#'
#' Computes the prospective reserve for an h-pay whole life policy.
#'
#' @param h Premium-paying period in years.
#' @inheritParams tVx
#'
#' @return A numeric vector of values.
#' @examples
#' htVx(40, h = 10, t = 5, i = 0.05, model = "uniform", omega = 100)
#' @export
htVx <- function(x, h, t, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  h <- .check_nonneg_integerish(h, "h")
  t <- .check_nonneg_integerish(t, "t")

  xht <- .recycle3_ch10(x, h, t, "x", "h", "t")
  x <- xht$a
  h <- xht$b
  t <- xht$c

  sapply(seq_along(x), function(j) {
    if (t[j] < h[j]) {
      Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...) -
        tPx(x = x[j], t = h[j], i = i, tbl = tbl, model = model, ...) *
        adotxn(x = x[j] + t[j], n = h[j] - t[j], i = i, tbl = tbl, model = model, ...)
    } else {
      Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)
    }
  })
}

# -------------------------------------------------------------------------
# Duration-t loss moments
# -------------------------------------------------------------------------

#' Mean present value of loss at duration t for whole life insurance
#'
#' Computes the conditional mean
#' \eqn{E[{}_tL_x \mid K_x \ge t]} for a fully discrete whole life insurance.
#'
#' Under the equivalence-principle premium, this equals the prospective
#' reserve \eqn{{}_tV_x}.
#'
#' @param x Issue age.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param P Annual premium.
#' @param model Optional parametric survival model name.
#' @param ... Additional model parameters.
#' @param tbl Optional life table object.
#'
#' @return A numeric vector of values.
#' @examples
#' prem <- Px(40, i = 0.05, model = "uniform", omega = 100)
#' ELtx(40, t = 10, i = 0.05, P = prem, model = "uniform", omega = 100)
#' @export
ELtx <- function(x, t, i, P, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")
  P <- .check_nonneg_numeric(P, "P")

  xtp <- .recycle3_ch10(x, t, P, "x", "t", "P")
  x <- xtp$a
  t <- xtp$b
  P <- xtp$c

  sapply(seq_along(x), function(j) {
    Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...) -
      P[j] * adotx(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)
  })
}

#' Variance of present value of loss at duration t for whole life insurance
#'
#' Computes the conditional variance
#' \eqn{\mathrm{Var}({}_tL_x \mid K_x \ge t)}
#' for a fully discrete whole life insurance.
#'
#' @inheritParams ELtx
#' @return A numeric vector of values.
#' @examples
#' prem <- Px(40, i = 0.05, model = "uniform", omega = 100)
#' varLtx(40, t = 10, i = 0.05, P = prem, model = "uniform", omega = 100)
#' @export
varLtx <- function(x, t, i, P, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")
  P <- .check_nonneg_numeric(P, "P")

  xtp <- .recycle3_ch10(x, t, P, "x", "t", "P")
  x <- xtp$a
  t <- xtp$b
  P <- xtp$c

  d <- i / (1 + i)

  sapply(seq_along(x), function(j) {
    (1 + P[j] / d)^2 *
      (A2x(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...) -
         Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)^2)
  })
}

# -------------------------------------------------------------------------
# Continuous and m-thly premium reserves
# -------------------------------------------------------------------------

#' Whole life reserve with continuous premiums
#'
#' Computes the reserve for a discrete whole life insurance
#' funded by continuous premiums.
#'
#' @inheritParams tVx
#' @return A numeric vector of values.
#' @examples
#' tVbarx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
tVbarx <- function(x, t, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")

  xt <- .recycle2_ch10(x, t, "x", "t")
  x <- xt$a
  t <- xt$b

  sapply(seq_along(x), function(j) {
    Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...) -
      Pbarx(x = x[j], i = i, tbl = tbl, model = model, ...) *
      abarx(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)
  })
}

#' Fully continuous whole life reserve
#'
#' Computes the reserve for a whole life insurance with
#' continuous premiums and immediate payment of claims.
#'
#' In the fully continuous setting, reserve time \eqn{t} may be any
#' nonnegative real value.
#'
#' @param x Issue age.
#' @param t Duration, allowed to be any nonnegative numeric value.
#' @param i Effective annual interest rate.
#' @param model Optional parametric survival model name.
#' @param ... Additional model parameters.
#' @param tbl Optional life table object.
#'
#' @return A numeric vector of values.
#' @examples
#' tVbarAbarx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#' tVbarAbarx(40, t = c(19, 19.25, 19.5, 19.75, 20), i = 0.06,
#'            model = "uniform", omega = 100)
#' @export
tVbarAbarx <- function(x, t, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_numeric(t, "t")

  xt <- .recycle2_ch10(x, t, "x", "t")
  x <- xt$a
  t <- xt$b

  sapply(seq_along(x), function(j) {
    Abarx(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...) -
      PbarAbarx(x = x[j], i = i, tbl = tbl, model = model, ...) *
      abarx(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)
  })
}

#' Whole life reserve with m-thly premiums
#'
#' Computes the reserve for a whole life insurance funded
#' by true m-thly premiums.
#'
#' @param m Number of premium payments per year.
#' @inheritParams tVx
#'
#' @return A numeric vector of values.
#' @examples
#' tVx_m(40, t = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#' @export
tVx_m <- function(x, t, m, i, model = NULL, ..., tbl = NULL) {
  .check_reserve_source(tbl, model)
  i <- .check_reserve_interest(i)
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")
  m <- .check_m_scalar_ch10(m)

  xt <- .recycle2_ch10(x, t, "x", "t")
  x <- xt$a
  t <- xt$b

  sapply(seq_along(x), function(j) {
    Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...) -
      Px_m(x = x[j], m = m, i = i, tbl = tbl, model = model, ...) *
      adotx_m(x = x[j] + t[j], m = m, i = i, tbl = tbl, model = model, ...)
  })
}

# -------------------------------------------------------------------------
# Gain / loss analysis
# -------------------------------------------------------------------------

#' Total gain for a discrete insurance contract
#'
#' Computes the total gain:
#' amount on hand at year-end minus amount required.
#'
#' @param Vt Reserve at duration t.
#' @param Vt1 Reserve at duration t+1.
#' @param P Net premium for the year.
#' @param i_actual Actual annual effective interest rate.
#' @param q_actual Actual mortality rate for the year.
#' @param B Benefit amount. Defaults to 1.
#'
#' @return A numeric vector of values.
#' @examples
#' GT_disc(Vt = 0.1, Vt1 = 0.11, P = 0.02, i_actual = 0.05, q_actual = 0.01)
#' @export
GT_disc <- function(Vt, Vt1, P, i_actual, q_actual, B = 1) {
  Vt <- .check_nonneg_numeric(Vt, "Vt")
  Vt1 <- .check_nonneg_numeric(Vt1, "Vt1")
  P <- .check_nonneg_numeric(P, "P")
  i_actual <- as.numeric(i_actual)
  q_actual <- as.numeric(q_actual)
  B <- .check_nonneg_numeric(B, "B")

  n <- max(length(Vt), length(Vt1), length(P), length(i_actual), length(q_actual), length(B))

  rep_len_checked <- function(z, nm) {
    if (!length(z) %in% c(1L, n)) stop(nm, " must have length 1 or common length.", call. = FALSE)
    rep(z, length.out = n)
  }

  Vt <- rep_len_checked(Vt, "Vt")
  Vt1 <- rep_len_checked(Vt1, "Vt1")
  P <- rep_len_checked(P, "P")
  i_actual <- rep_len_checked(i_actual, "i_actual")
  q_actual <- rep_len_checked(q_actual, "q_actual")
  B <- rep_len_checked(B, "B")

  if (any(!is.finite(i_actual)) || any(i_actual <= -1)) {
    stop("i_actual must be finite and greater than -1.", call. = FALSE)
  }
  if (any(!is.finite(q_actual)) || any(q_actual < 0) || any(q_actual > 1)) {
    stop("q_actual must lie in [0, 1].", call. = FALSE)
  }

  (Vt + P) * (1 + i_actual) - (q_actual * B + (1 - q_actual) * Vt1)
}

#' Mortality gain for a discrete insurance contract
#'
#' @param i_assumed Assumed annual effective interest rate.
#' @inheritParams GT_disc
#'
#' @return A numeric vector of values.
#' @examples
#' GM_disc(Vt = 0.1, Vt1 = 0.11, P = 0.02, i_assumed = 0.04, q_actual = 0.01)
#' @export
GM_disc <- function(Vt, Vt1, P, i_assumed, q_actual, B = 1) {
  GT_disc(
    Vt = Vt,
    Vt1 = Vt1,
    P = P,
    i_actual = i_assumed,
    q_actual = q_actual,
    B = B
  )
}

#' Interest gain for a discrete insurance contract
#'
#' @param i_actual Actual annual effective interest rate.
#' @param q_assumed Assumed mortality rate for the year.
#' @inheritParams GT_disc
#'
#' @return A numeric vector of values.
#' @examples
#' GI_disc(Vt = 0.1, Vt1 = 0.11, P = 0.02, i_actual = 0.05, q_assumed = 0.01)
#' @export
GI_disc <- function(Vt, Vt1, P, i_actual, q_assumed, B = 1) {
  GT_disc(
    Vt = Vt,
    Vt1 = Vt1,
    P = P,
    i_actual = i_actual,
    q_actual = q_assumed,
    B = B
  )
}

