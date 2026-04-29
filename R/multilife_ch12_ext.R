#' Extended multi-life functions for Chapter 12
#'
#' Additional Chapter 12 functions for contingent probabilities,
#' continuous contingent insurances, and continuous multi-life annuities.
#'
#' @name multilife_ch12_ext
#' @keywords internal
NULL

# -------------------------------------------------------------------------
# Shared documentation block
# -------------------------------------------------------------------------

#' Shared parameters for Chapter 12 continuous multi-life functions
#'
#' @param x Age of first life.
#' @param y Age of second life.
#' @param i Effective annual interest rate.
#' @param n Term in years.
#' @param tbl Life table.
#' @param model Survival model.
#' @param ... Additional model parameters.
#'
#' @name ch12_multilife_cont_params
#' @keywords internal
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.check_delta_ch12 <- function(i) {
  i <- as.numeric(i)
  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }
  log(1 + i)
}

# -------------------------------------------------------------------------
# Contingent probability functions
# -------------------------------------------------------------------------

#' Probability that (x) fails before (y) within n years
#'
#' Computes \eqn{{}_n q_{xy}^{1} = \int_0^n {}_tp_{xy}\mu_{x+t}\,dt}
#' under independence.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' tqxy1(40, 50, n = 10, model = "uniform", omega = 100)
#' @export
tqxy1 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  n <- .check_nonneg_numeric(n, "n")

  xyn <- .recycle3_ch12(x, y, n, "x", "y", "n")
  x <- xyn$a
  y <- xyn$b
  n <- xyn$c

  sapply(seq_along(x), function(j) {
    integrate(
      function(t) {
        tpxy(x = x[j], y = y[j], t = t, tbl = tbl, model = model, ...) *
          hazard0(x[j] + t, tbl = tbl, model = model, ...)
      },
      lower = 0, upper = n[j]
    )$value
  })
}

#' Probability that (y) fails before (x) within n years
#'
#' Computes \eqn{{}_n q_{xy}^{\hspace{1mm}1} = \int_0^n {}_tp_{xy}\mu_{y+t}\,dt}
#' under independence.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' tqyx1(40, 50, n = 10, model = "uniform", omega = 100)
#' @export
tqyx1 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  n <- .check_nonneg_numeric(n, "n")

  xyn <- .recycle3_ch12(x, y, n, "x", "y", "n")
  x <- xyn$a
  y <- xyn$b
  n <- xyn$c

  sapply(seq_along(x), function(j) {
    integrate(
      function(t) {
        tpxy(x = x[j], y = y[j], t = t, tbl = tbl, model = model, ...) *
          hazard0(y[j] + t, tbl = tbl, model = model, ...)
      },
      lower = 0, upper = n[j]
    )$value
  })
}

#' Probability that (x) fails after (y) within n years
#'
#' Computes \eqn{{}_n q_{xy}^{2} = {}_n q_x - {}_n q_{xy}^{1}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' tqxy2(40, 50, n = 10, model = "uniform", omega = 100)
#' @export
tqxy2 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  tqx(x = x, t = n, tbl = tbl, model = model, ...) -
    tqxy1(x = x, y = y, n = n, tbl = tbl, model = model, ...)
}

#' Probability that (y) fails after (x) within n years
#'
#' Computes \eqn{{}_n q_{xy}^{\hspace{1mm}2} = {}_n q_y - {}_n q_{xy}^{\hspace{1mm}1}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' tqyx2(40, 50, n = 10, model = "uniform", omega = 100)
#' @export
tqyx2 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  tqx(x = y, t = n, tbl = tbl, model = model, ...) -
    tqyx1(x = x, y = y, n = n, tbl = tbl, model = model, ...)
}

# -------------------------------------------------------------------------
# Continuous annuities
# -------------------------------------------------------------------------

#' Continuous joint-life whole life annuity
#'
#' Computes \eqn{\overline{a}_{xy} = \int_0^\infty v^t\,{}_tp_{xy}\,dt}.
#'
#' Shared documentation topic used to avoid filename collisions with
#' case-distinct function names on case-insensitive file systems.
#'
#' @name abarxy_ch12
#' @aliases abarxy
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' abarxy(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
abarxy <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  i <- .check_i_ch12(i)
  delta <- .check_delta_ch12(i)

  xyi <- .recycle3_ch12(x, y, delta, "x", "y", "i")
  x <- xyi$a
  y <- xyi$b
  delta <- xyi$c

  sapply(seq_along(x), function(j) {
    integrate(
      function(t) exp(-delta[j] * t) *
        tpxy(x = x[j], y = y[j], t = t, tbl = tbl, model = model, ...),
      lower = 0, upper = Inf
    )$value
  })
}

#' Continuous last-survivor whole life annuity
#'
#' Computes \eqn{\overline{a}_{\overline{xy}} = \overline{a}_x + \overline{a}_y - \overline{a}_{xy}}.
#'
#' Shared documentation topic used to avoid filename collisions with
#' case-distinct function names on case-insensitive file systems.
#'
#' @name abarxybar_ch12
#' @aliases abarxybar
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' abarxybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
abarxybar <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  abarx(x = x, i = i, tbl = tbl, model = model, ...) +
    abarx(x = y, i = i, tbl = tbl, model = model, ...) -
    abarxy(x = x, y = y, i = i, tbl = tbl, model = model, ...)
}

#' Continuous reversionary annuity to (y) after death of (x)
#'
#' Computes \eqn{\overline{a}_{x|y} = \overline{a}_y - \overline{a}_{xy}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' abarx_y(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
abarx_y <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  abarx(x = y, i = i, tbl = tbl, model = model, ...) -
    abarxy(x = x, y = y, i = i, tbl = tbl, model = model, ...)
}

#' Continuous reversionary annuity to (x) after death of (y)
#'
#' Computes \eqn{\overline{a}_{y|x} = \overline{a}_x - \overline{a}_{xy}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' abary_x(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
abary_x <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  abarx(x = x, i = i, tbl = tbl, model = model, ...) -
    abarxy(x = x, y = y, i = i, tbl = tbl, model = model, ...)
}

# -------------------------------------------------------------------------
# Continuous insurances
# -------------------------------------------------------------------------

#' Continuous joint-life whole life insurance
#'
#' Computes \eqn{\overline{A}_{xy} = 1 - \delta \overline{a}_{xy}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' Abarxy(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Abarxy <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  delta <- .check_delta_ch12(i)
  1 - delta * abarxy(x = x, y = y, i = i, tbl = tbl, model = model, ...)
}

#' Continuous last-survivor whole life insurance
#'
#' Computes \eqn{\overline{A}_{\overline{xy}} = \overline{A}_x + \overline{A}_y - \overline{A}_{xy}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' Abarxybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Abarxybar <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  Abarx(x = x, i = i, tbl = tbl, model = model, ...) +
    Abarx(x = y, i = i, tbl = tbl, model = model, ...) -
    Abarxy(x = x, y = y, i = i, tbl = tbl, model = model, ...)
}

#' Continuous contingent insurance: benefit on death of (x) if before (y)
#'
#' Computes \eqn{\overline{A}_{xy}^{1} = \int_0^\infty v^t\,{}_tp_{xy}\mu_{x+t}\,dt}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' Abarxy1(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Abarxy1 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  i <- .check_i_ch12(i)
  delta <- .check_delta_ch12(i)

  dots <- list(...)
  xyi <- .recycle3_ch12(x, y, delta, "x", "y", "i")
  x <- xyi$a
  y <- xyi$b
  delta <- xyi$c

  sapply(seq_along(x), function(j) {
    upper <- Inf
    if (!is.null(dots$omega) && is.finite(dots$omega)) {
      upper <- max(0, min(dots$omega - x[j], dots$omega - y[j]))
    }

    integrate(
      function(t) {
        exp(-delta[j] * t) *
          tpxy(x = x[j], y = y[j], t = t, tbl = tbl, model = model, ...) *
          hazard0(x[j] + t, tbl = tbl, model = model, ...)
      },
      lower = 0,
      upper = upper
    )$value
  })
}

#' Continuous contingent insurance: benefit on death of (y) if before (x)
#'
#' Computes \eqn{\overline{A}_{xy}^{\hspace{1mm}1} = \int_0^\infty v^t\,{}_tp_{xy}\mu_{y+t}\,dt}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' Abaryx1(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Abaryx1 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  y <- .check_nonneg_numeric(y, "y")
  i <- .check_i_ch12(i)
  delta <- .check_delta_ch12(i)

  dots <- list(...)
  xyi <- .recycle3_ch12(x, y, delta, "x", "y", "i")
  x <- xyi$a
  y <- xyi$b
  delta <- xyi$c

  sapply(seq_along(x), function(j) {
    upper <- Inf
    if (!is.null(dots$omega) && is.finite(dots$omega)) {
      upper <- max(0, min(dots$omega - x[j], dots$omega - y[j]))
    }

    integrate(
      function(t) {
        exp(-delta[j] * t) *
          tpxy(x = x[j], y = y[j], t = t, tbl = tbl, model = model, ...) *
          hazard0(y[j] + t, tbl = tbl, model = model, ...)
      },
      lower = 0,
      upper = upper
    )$value
  })
}

#' Continuous contingent insurance: benefit on death of (x) if after (y)
#'
#' Computes \eqn{\overline{A}_{xy}^{2} = \overline{A}_x - \overline{A}_{xy}^{1}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' Abarxy2(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Abarxy2 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  Abarx(x = x, i = i, tbl = tbl, model = model, ...) -
    Abarxy1(x = x, y = y, i = i, tbl = tbl, model = model, ...)
}

#' Continuous contingent insurance: benefit on death of (y) if after (x)
#'
#' Computes \eqn{\overline{A}_{xy}^{\hspace{1mm}2} = \overline{A}_y - \overline{A}_{xy}^{\hspace{1mm}1}}.
#'
#' @inheritParams ch12_multilife_cont_params
#' @return Numeric vector.
#' @examples
#' Abaryx2(40, 50, i = 0.05, model = "uniform", omega = 100)
#' @export
Abaryx2 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  Abarx(x = y, i = i, tbl = tbl, model = model, ...) -
    Abaryx1(x = x, y = y, i = i, tbl = tbl, model = model, ...)
}
