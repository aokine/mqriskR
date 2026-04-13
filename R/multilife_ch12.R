#' Multi-life functions for Chapter 12
#'
#' Core Chapter 12 functions for two-life joint-life and last-survivor models,
#' including survival probabilities, annuities, insurances, and reversionary annuities.
#'
#' @name multilife_ch12
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.recycle3_ch12 <- function(a, b, c,
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
.check_i_ch12 <- function(i) {
  i <- as.numeric(i)
  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }
  i
}

#' @noRd
.sum_until_tol_ch12 <- function(fun, k_max = 5000L, tol = 1e-12) {
  ks <- 0:(k_max - 1L)
  vals <- fun(ks)
  vals[!is.finite(vals)] <- 0
  s <- sum(vals)

  if (length(vals) > 0L && abs(tail(vals, 1L)) > tol) {
    warning("Series may not have converged; consider increasing k_max.", call. = FALSE)
  }

  s
}

# -------------------------------------------------------------------------
# Joint-life and last-survivor survival probabilities
# -------------------------------------------------------------------------

#' Joint-life survival probability
#'
#' Computes \eqn{{}_tp_{xy} = {}_tp_x\,{}_tp_y} under independence.
#'
#' @param x Age of first life.
#' @param y Age of second life.
#' @param t Time.
#' @param tbl Life table.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tpxy(40, 50, t = 10, model = "uniform", omega = 100)
#' @export
tpxy <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  t <- .check_nonneg_numeric(t, "t")

  xyt <- .recycle3_ch12(x, y, t, "x", "y", "t")
  x <- xyt$a
  y <- xyt$b
  t <- xyt$c

  tpx(x = x, t = t, tbl = tbl, model = model, ...) *
    tpx(x = y, t = t, tbl = tbl, model = model, ...)
}

#' Joint-life failure probability
#'
#' Computes \eqn{{}_tq_{xy} = 1 - {}_tp_{xy}}.
#'
#' @inheritParams tpxy
#' @return Numeric vector.
#' @examples
#' tqxy(40, 50, t = 10, model = "uniform", omega = 100)
#' @export
tqxy <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  1 - tpxy(x = x, y = y, t = t, tbl = tbl, model = model, ...)
}

#' Last-survivor survival probability
#'
#' Computes \eqn{{}_tp_{\overline{xy}} = {}_tp_x + {}_tp_y - {}_tp_{xy}}.
#'
#' @inheritParams tpxy
#' @return Numeric vector.
#' @examples
#' tpxybar(40, 50, t = 10, model = "uniform", omega = 100)
#' @export
tpxybar <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  t <- .check_nonneg_numeric(t, "t")

  xyt <- .recycle3_ch12(x, y, t, "x", "y", "t")
  x <- xyt$a
  y <- xyt$b
  t <- xyt$c

  tpx(x = x, t = t, tbl = tbl, model = model, ...) +
    tpx(x = y, t = t, tbl = tbl, model = model, ...) -
    tpxy(x = x, y = y, t = t, tbl = tbl, model = model, ...)
}

#' Last-survivor failure probability
#'
#' Computes \eqn{{}_tq_{\overline{xy}} = 1 - {}_tp_{\overline{xy}}}.
#'
#' @inheritParams tpxy
#' @return Numeric vector.
#' @examples
#' tqxybar(40, 50, t = 10, model = "uniform", omega = 100)
#' @export
tqxybar <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  1 - tpxybar(x = x, y = y, t = t, tbl = tbl, model = model, ...)
}

# -------------------------------------------------------------------------
# Pure endowments
# -------------------------------------------------------------------------

#' Joint-life pure endowment
#'
#' Computes \eqn{{}_nE_{xy} = v^n\,{}_np_{xy}}.
#'
#' @param n Term.
#' @inheritParams tpxy
#' @param i Effective annual interest rate.
#'
#' @return Numeric vector.
#' @examples
#' nExy(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
nExy <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  n <- .check_nonneg_integerish(n, "n")
  i <- .check_i_ch12(i)

  xyn <- .recycle3_ch12(x, y, n, "x", "y", "n")
  x <- xyn$a
  y <- xyn$b
  n <- xyn$c

  (1 / (1 + i))^n * tpxy(x = x, y = y, t = n, tbl = tbl, model = model, ...)
}

#' Last-survivor pure endowment
#'
#' Computes \eqn{{}_nE_{\overline{xy}} = v^n\,{}_np_{\overline{xy}}}.
#'
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' nExybar(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
nExybar <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  n <- .check_nonneg_integerish(n, "n")
  i <- .check_i_ch12(i)

  xyn <- .recycle3_ch12(x, y, n, "x", "y", "n")
  x <- xyn$a
  y <- xyn$b
  n <- xyn$c

  (1 / (1 + i))^n * tpxybar(x = x, y = y, t = n, tbl = tbl, model = model, ...)
}

# -------------------------------------------------------------------------
# Joint-life annuities and insurances
# -------------------------------------------------------------------------

#' Joint-life temporary annuity-due
#'
#' Computes \eqn{\ddot{a}_{xy:\overline{n}|}}.
#'
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' adotxyn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
adotxyn <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  n <- .check_nonneg_integerish(n, "n")
  i <- .check_i_ch12(i)

  xyn <- .recycle3_ch12(x, y, n, "x", "y", "n")
  x <- xyn$a
  y <- xyn$b
  n <- xyn$c

  sapply(seq_along(x), function(j) {
    k <- 0:(n[j] - 1L)
    sum((1 / (1 + i[j]))^k *
          tpxy(x = x[j], y = y[j], t = k, tbl = tbl, model = model, ...))
  })
}

#' Joint-life temporary annuity-immediate
#'
#' Computes \eqn{a_{xy:\overline{n}|}}.
#'
#' Shared documentation topic used to avoid filename collisions with
#' case-distinct function names on case-insensitive file systems.
#'
#' @name axyn_ch12
#' @aliases axyn
#' @inheritParams adotxyn
#' @return Numeric vector.
#' @examples
#' axyn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
axyn <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  adotxyn(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...) -
    nExy(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...)
}

#' Joint-life whole life annuity-due
#'
#' Computes \eqn{\ddot{a}_{xy}}.
#'
#' @param k_max Maximum number of terms.
#' @param tol Convergence tolerance.
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' adotxy(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
adotxy <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  i <- .check_i_ch12(i)

  xyi <- .recycle3_ch12(x, y, i, "x", "y", "i")
  x <- xyi$a
  y <- xyi$b
  i <- xyi$c

  sapply(seq_along(x), function(j) {
    .sum_until_tol_ch12(
      fun = function(k) {
        (1 / (1 + i[j]))^k *
          tpxy(x = x[j], y = y[j], t = k, tbl = tbl, model = model, ...)
      },
      k_max = k_max,
      tol = tol
    )
  })
}

#' Joint-life whole life annuity-immediate
#'
#' Computes \eqn{a_{xy}}.
#'
#' Shared documentation topic used to avoid filename collisions with
#' case-distinct function names on case-insensitive file systems.
#'
#' @name axy_ch12
#' @aliases axy
#' @inheritParams adotxy
#' @return Numeric vector.
#' @examples
#' axy(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
axy <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  adotxy(
    x = x, y = y, i = i, tbl = tbl, model = model, ...,
    k_max = k_max, tol = tol
  ) - 1
}

#' Joint-life term insurance
#'
#' Computes \eqn{A^{1}_{xy:\overline{n}|}}.
#'
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' Axyn1(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
Axyn1 <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  d <- i / (1 + i)
  1 - d * adotxyn(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...) -
    nExy(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...)
}

#' Joint-life endowment insurance
#'
#' Computes \eqn{A_{xy:\overline{n}|}}.
#'
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' Axyn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
Axyn <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  Axyn1(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...) +
    nExy(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...)
}

#' Joint-life whole life insurance
#'
#' Computes \eqn{A_{xy}}.
#'
#' @inheritParams adotxy
#' @return Numeric vector.
#' @examples
#' Axy(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Axy <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  d <- i / (1 + i)
  1 - d * adotxy(
    x = x, y = y, i = i, tbl = tbl, model = model, ...,
    k_max = k_max, tol = tol
  )
}

# -------------------------------------------------------------------------
# Last-survivor annuities and insurances
# -------------------------------------------------------------------------

#' Last-survivor temporary annuity-due
#'
#' Computes \eqn{\ddot{a}_{\overline{xy}:\overline{n}|}}.
#'
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' adotxybarn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
adotxybarn <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  adotxn(x = x, n = n, i = i, tbl = tbl, model = model, ...) +
    adotxn(x = y, n = n, i = i, tbl = tbl, model = model, ...) -
    adotxyn(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...)
}

#' Last-survivor temporary annuity-immediate
#'
#' Computes \eqn{a_{\overline{xy}:\overline{n}|}}.
#'
#' Shared documentation topic used to avoid filename collisions with
#' case-distinct function names on case-insensitive file systems.
#'
#' @name axybarn_ch12
#' @aliases axybarn
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' axybarn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
axybarn <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  axn(x = x, n = n, i = i, tbl = tbl, model = model, ...) +
    axn(x = y, n = n, i = i, tbl = tbl, model = model, ...) -
    axyn(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...)
}

#' Last-survivor whole life annuity-due
#'
#' Computes \eqn{\ddot{a}_{\overline{xy}}}.
#'
#' @inheritParams adotxy
#' @return Numeric vector.
#' @examples
#' adotxybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
adotxybar <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  adotx(x = x, i = i, tbl = tbl, model = model, ...) +
    adotx(x = y, i = i, tbl = tbl, model = model, ...) -
    adotxy(x = x, y = y, i = i, tbl = tbl, model = model, ...,
           k_max = k_max, tol = tol)
}

#' Last-survivor whole life annuity-immediate
#'
#' Computes \eqn{a_{\overline{xy}}}.
#'
#' Shared documentation topic used to avoid filename collisions with
#' case-distinct function names on case-insensitive file systems.
#'
#' @name axybar_ch12
#' @aliases axybar
#' @inheritParams adotxybar
#' @return Numeric vector.
#' @examples
#' axybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
axybar <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  ax(x = x, i = i, tbl = tbl, model = model, ...) +
    ax(x = y, i = i, tbl = tbl, model = model, ...) -
    axy(x = x, y = y, i = i, tbl = tbl, model = model, ...,
        k_max = k_max, tol = tol)
}

#' Last-survivor term insurance
#'
#' Computes \eqn{A^{1}_{\overline{xy}:\overline{n}|}}.
#'
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' Axybarn1(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
Axybarn1 <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  Axn1(x = x, n = n, i = i, tbl = tbl, model = model, ...) +
    Axn1(x = y, n = n, i = i, tbl = tbl, model = model, ...) -
    Axyn1(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...)
}

#' Last-survivor endowment insurance
#'
#' Computes \eqn{A_{\overline{xy}:\overline{n}|}}.
#'
#' @inheritParams nExy
#' @return Numeric vector.
#' @examples
#' Axybarn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
Axybarn <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  Axybarn1(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...) +
    nExybar(x = x, y = y, n = n, i = i, tbl = tbl, model = model, ...)
}

#' Last-survivor whole life insurance
#'
#' Computes \eqn{A_{\overline{xy}} = A_x + A_y - A_{xy}}.
#'
#' @inheritParams adotxy
#' @return Numeric vector.
#' @examples
#' Axybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Axybar <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  Ax(x = x, i = i, tbl = tbl, model = model, ...) +
    Ax(x = y, i = i, tbl = tbl, model = model, ...) -
    Axy(x = x, y = y, i = i, tbl = tbl, model = model, ...,
        k_max = k_max, tol = tol)
}

# -------------------------------------------------------------------------
# Reversionary annuities
# -------------------------------------------------------------------------

#' Reversionary annuity to (y) after death of (x)
#'
#' Computes \eqn{a_{x|y} = a_y - a_{xy}}.
#'
#' @inheritParams adotxy
#' @return Numeric vector.
#' @examples
#' ax_y(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
ax_y <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  ax(x = y, i = i, tbl = tbl, model = model, ...) -
    axy(x = x, y = y, i = i, tbl = tbl, model = model, ...,
        k_max = k_max, tol = tol)
}

#' Reversionary annuity to (x) after death of (y)
#'
#' Computes \eqn{a_{y|x} = a_x - a_{xy}}.
#'
#' @inheritParams adotxy
#' @return Numeric vector.
#' @examples
#' ay_x(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
ay_x <- function(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12) {
  ax(x = x, i = i, tbl = tbl, model = model, ...) -
    axy(x = x, y = y, i = i, tbl = tbl, model = model, ...,
        k_max = k_max, tol = tol)
}
