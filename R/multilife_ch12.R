# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.multilife_check_basis <- function(tbl = NULL, model = NULL) {
  if (is.null(tbl) && is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }

  if (!is.null(tbl) && !is.null(model)) {
    stop("Supply only one of tbl or model, not both.", call. = FALSE)
  }

  if (!is.null(tbl)) {
    .validate_life_table(tbl)
  }

  invisible(TRUE)
}

#' @noRd
.multilife_check_nonnegative <- function(x, name) {
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
.multilife_check_integerish <- function(x, name) {
  x <- .multilife_check_nonnegative(x, name)

  if (any(abs(x - round(x)) > 1e-10)) {
    stop(
      name,
      " must contain nonnegative integer-like values.",
      call. = FALSE
    )
  }

  as.integer(round(x))
}

#' @noRd
.multilife_check_interest <- function(i) {
  i <- as.numeric(i)

  if (length(i) == 0L) {
    stop("i must have positive length.", call. = FALSE)
  }

  if (any(!is.finite(i)) || any(i <= -1)) {
    stop(
      "i must contain finite values greater than -1.",
      call. = FALSE
    )
  }

  i
}

#' @noRd
.multilife_check_sum_controls <- function(k_max, tol) {
  k_max <- as.numeric(k_max)
  tol <- as.numeric(tol)

  if (length(k_max) != 1L ||
      !is.finite(k_max) ||
      k_max <= 0 ||
      abs(k_max - round(k_max)) > 1e-10) {
    stop("k_max must be a positive integer.", call. = FALSE)
  }

  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite number.", call. = FALSE)
  }

  list(
    k_max = as.integer(round(k_max)),
    tol = tol
  )
}

#' @noRd
.multilife_recycle <- function(..., .names = NULL) {
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

  for (j in seq_along(values)) {
    if (!lengths[j] %in% c(1L, common_length)) {
      stop(
        .names[j],
        " must have length 1 or the common length ",
        common_length,
        ".",
        call. = FALSE
      )
    }

    values[[j]] <- rep(values[[j]], length.out = common_length)
  }

  values
}

#' @noRd
.multilife_table_survival <- function(tbl, x, t) {
  x <- as.numeric(x)
  t <- as.numeric(t)

  values <- .multilife_recycle(
    x,
    t,
    .names = c("x", "t")
  )

  x <- values[[1]]
  t <- values[[2]]

  out <- numeric(length(x))

  # Survival for zero years is always one, including at the limiting age.
  zero_time <- t == 0
  out[zero_time] <- 1

  positive_time <- !zero_time

  if (any(positive_time)) {
    out[positive_time] <- npx(
      tbl,
      x = x[positive_time],
      n = t[positive_time]
    )
  }

  out[!is.finite(out)] <- 0
  pmin(pmax(out, 0), 1)
}

#' @noRd
.multilife_single_survival <- function(x, t, tbl = NULL, model = NULL, ...) {
  .multilife_check_basis(tbl, model)

  values <- .multilife_recycle(
    x,
    t,
    .names = c("x", "t")
  )

  x <- values[[1]]
  t <- values[[2]]

  if (!is.null(tbl)) {
    return(.multilife_table_survival(tbl, x, t))
  }

  out <- tpx(
    x = x,
    t = t,
    model = model,
    ...
  )

  pmin(pmax(out, 0), 1)
}

#' @noRd
.multilife_joint_table_horizon <- function(tbl, x, y) {
  maximum_age <- max(tbl$x, na.rm = TRUE)

  lx_x <- lx(tbl, x = x)
  lx_y <- lx(tbl, x = y)

  if (is.na(lx_x) || is.na(lx_y)) {
    if (x > maximum_age || y > maximum_age) {
      return(-1L)
    }

    stop(
      "One or both ages are not available in the life table.",
      call. = FALSE
    )
  }

  if (lx_x <= 0 || lx_y <= 0) {
    return(-1L)
  }

  max(
    floor(min(maximum_age - x, maximum_age - y)),
    0L
  )
}

#' @noRd
.multilife_sum_series <- function(fun, start, end, tol,
                                  warn_if_unconverged = FALSE) {
  if (end < start) {
    return(0)
  }

  times <- seq.int(start, end)
  values <- fun(times)

  values[!is.finite(values)] <- 0
  total <- sum(values)

  if (warn_if_unconverged &&
      length(values) > 0L &&
      abs(tail(values, 1L)) > tol) {
    warning(
      "Series may not have converged; consider increasing k_max.",
      call. = FALSE
    )
  }

  total
}

# -------------------------------------------------------------------------
# Joint-life and last-survivor survival probabilities
# -------------------------------------------------------------------------

#' Multi-life survival and failure probabilities
#'
#' Computes joint-life and last-survivor survival and failure probabilities
#' for two independent lives.
#'
#' \code{tpxy()} computes
#' \eqn{{}_t p_{xy} = {}_t p_x\,{}_t p_y}.
#'
#' \code{tqxy()} computes
#' \eqn{{}_t q_{xy} = 1 - {}_t p_{xy}}.
#'
#' \code{tpxybar()} computes
#' \eqn{{}_t p_{\overline{xy}}
#' = {}_t p_x + {}_t p_y - {}_t p_{xy}}.
#'
#' \code{tqxybar()} computes
#' \eqn{{}_t q_{\overline{xy}}
#' = 1 - {}_t p_{\overline{xy}}}.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param t Nonnegative duration. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the survival model.
#'
#' @return A numeric vector of probabilities.
#'
#' @examples
#' tpxy(40, 50, t = 10, model = "uniform", omega = 100)
#' tqxy(40, 50, t = 10, model = "uniform", omega = 100)
#' tpxybar(40, 50, t = 10, model = "uniform", omega = 100)
#' tqxybar(40, 50, t = 10, model = "uniform", omega = 100)
#'
#' @rdname multilife_survival
#' @export
tpxy <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  .multilife_check_basis(tbl, model)

  x <- .multilife_check_nonnegative(x, "x")
  y <- .multilife_check_nonnegative(y, "y")
  t <- .multilife_check_nonnegative(t, "t")

  values <- .multilife_recycle(
    x,
    y,
    t,
    .names = c("x", "y", "t")
  )

  x <- values[[1]]
  y <- values[[2]]
  t <- values[[3]]

  out <- .multilife_single_survival(
    x = x,
    t = t,
    tbl = tbl,
    model = model,
    ...
  ) * .multilife_single_survival(
    x = y,
    t = t,
    tbl = tbl,
    model = model,
    ...
  )

  pmin(pmax(out, 0), 1)
}

#' @rdname multilife_survival
#' @export
tqxy <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  1 - tpxy(
    x = x,
    y = y,
    t = t,
    tbl = tbl,
    model = model,
    ...
  )
}

#' @rdname multilife_survival
#' @export
tpxybar <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  .multilife_check_basis(tbl, model)

  x <- .multilife_check_nonnegative(x, "x")
  y <- .multilife_check_nonnegative(y, "y")
  t <- .multilife_check_nonnegative(t, "t")

  values <- .multilife_recycle(
    x,
    y,
    t,
    .names = c("x", "y", "t")
  )

  x <- values[[1]]
  y <- values[[2]]
  t <- values[[3]]

  px <- .multilife_single_survival(
    x = x,
    t = t,
    tbl = tbl,
    model = model,
    ...
  )

  py <- .multilife_single_survival(
    x = y,
    t = t,
    tbl = tbl,
    model = model,
    ...
  )

  out <- px + py - px * py

  pmin(pmax(out, 0), 1)
}

#' @rdname multilife_survival
#' @export
tqxybar <- function(x, y, t, tbl = NULL, model = NULL, ...) {
  1 - tpxybar(
    x = x,
    y = y,
    t = t,
    tbl = tbl,
    model = model,
    ...
  )
}

# -------------------------------------------------------------------------
# Pure endowments
# -------------------------------------------------------------------------

#' Multi-life pure endowments
#'
#' Computes pure endowment actuarial present values for joint-life and
#' last-survivor statuses.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param n Nonnegative integer term. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the survival model.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @examples
#' nExy(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#' nExybar(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#'
#' @rdname multilife_pure_endowment
#' @export
nExy <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  .multilife_check_basis(tbl, model)

  x <- .multilife_check_nonnegative(x, "x")
  y <- .multilife_check_nonnegative(y, "y")
  n <- .multilife_check_integerish(n, "n")
  i <- .multilife_check_interest(i)

  values <- .multilife_recycle(
    x,
    y,
    n,
    i,
    .names = c("x", "y", "n", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  n <- values[[3]]
  i <- values[[4]]

  discount(i, n) * tpxy(
    x = x,
    y = y,
    t = n,
    tbl = tbl,
    model = model,
    ...
  )
}

#' @rdname multilife_pure_endowment
#' @export
nExybar <- function(x, y, n, i, tbl = NULL, model = NULL, ...) {
  .multilife_check_basis(tbl, model)

  x <- .multilife_check_nonnegative(x, "x")
  y <- .multilife_check_nonnegative(y, "y")
  n <- .multilife_check_integerish(n, "n")
  i <- .multilife_check_interest(i)

  values <- .multilife_recycle(
    x,
    y,
    n,
    i,
    .names = c("x", "y", "n", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  n <- values[[3]]
  i <- values[[4]]

  discount(i, n) * tpxybar(
    x = x,
    y = y,
    t = n,
    tbl = tbl,
    model = model,
    ...
  )
}

# -------------------------------------------------------------------------
# Joint-life annuities
# -------------------------------------------------------------------------

#' Joint-life annuities
#'
#' Computes temporary and whole-life joint-life annuities-due and
#' annuities-immediate for two independent lives.
#'
#' \code{adotxyn()} computes an \eqn{n}-year temporary joint-life
#' annuity-due.
#'
#' \code{axyn()} computes an \eqn{n}-year temporary joint-life
#' annuity-immediate.
#'
#' \code{adotxy()} computes a whole-life joint-life annuity-due.
#'
#' \code{axy()} computes a whole-life joint-life annuity-immediate.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param n Nonnegative integer term. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the survival model.
#' @param k_max Maximum summation horizon used for non-terminating parametric
#'   survival models.
#' @param tol Positive convergence tolerance for whole-life summations.
#'
#' @return A numeric vector of annuity actuarial present values.
#'
#' @details
#' All calculations assume independence between the two future lifetimes.
#'
#' The temporary annuity-due is
#' \deqn{
#' \ddot{a}_{xy:\overline{n}|}
#' =
#' \sum_{k=0}^{n-1} v^k\,{}_kp_{xy}.
#' }
#'
#' The temporary annuity-immediate is
#' \deqn{
#' a_{xy:\overline{n}|}
#' =
#' \sum_{k=1}^{n} v^k\,{}_kp_{xy}.
#' }
#'
#' For life-table calculations, the whole-life sums terminate at the last
#' duration supported jointly by the two lives. For parametric models, the
#' sums are evaluated until convergence or until \code{k_max} is reached.
#'
#' @examples
#' adotxyn(
#'   x = 40,
#'   y = 50,
#'   n = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' axyn(
#'   x = 40,
#'   y = 50,
#'   n = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' adotxy(
#'   x = 40,
#'   y = 50,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' axy(
#'   x = 40,
#'   y = 50,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' @name joint_life_annuities
#' @rdname joint_life_annuities
#' @export
adotxyn <- function(x, y, n, i,
                    tbl = NULL, model = NULL, ...) {
  .multilife_check_basis(tbl, model)

  x <- .multilife_check_nonnegative(x, "x")
  y <- .multilife_check_nonnegative(y, "y")
  n <- .multilife_check_integerish(n, "n")
  i <- .multilife_check_interest(i)

  values <- .multilife_recycle(
    x,
    y,
    n,
    i,
    .names = c("x", "y", "n", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  n <- values[[3]]
  i <- values[[4]]

  vapply(
    seq_along(x),
    function(j) {
      if (n[j] == 0L) {
        return(0)
      }

      times <- 0:(n[j] - 1L)

      sum(
        discount(i[j], times) *
          tpxy(
            x = x[j],
            y = y[j],
            t = times,
            tbl = tbl,
            model = model,
            ...
          )
      )
    },
    numeric(1)
  )
}


#' @rdname joint_life_annuities
#' @export
axyn <- function(x, y, n, i,
                 tbl = NULL, model = NULL, ...) {
  adotxyn(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - 1 +
    nExy(
      x = x,
      y = y,
      n = n,
      i = i,
      tbl = tbl,
      model = model,
      ...
    )
}


#' @rdname joint_life_annuities
#' @export
adotxy <- function(x, y, i,
                   tbl = NULL, model = NULL, ...,
                   k_max = 5000L, tol = 1e-12) {
  .multilife_check_basis(tbl, model)

  x <- .multilife_check_nonnegative(x, "x")
  y <- .multilife_check_nonnegative(y, "y")
  i <- .multilife_check_interest(i)

  controls <- .multilife_check_sum_controls(k_max, tol)
  k_max <- controls$k_max
  tol <- controls$tol

  values <- .multilife_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  vapply(
    seq_along(x),
    function(j) {
      if (!is.null(tbl)) {
        upper <- .multilife_joint_table_horizon(
          tbl = tbl,
          x = x[j],
          y = y[j]
        )

        if (upper < 0L) {
          return(0)
        }

        return(
          .multilife_sum_series(
            fun = function(k) {
              discount(i[j], k) *
                tpxy(
                  x = x[j],
                  y = y[j],
                  t = k,
                  tbl = tbl
                )
            },
            start = 0L,
            end = upper,
            tol = tol
          )
        )
      }

      .multilife_sum_series(
        fun = function(k) {
          discount(i[j], k) *
            tpxy(
              x = x[j],
              y = y[j],
              t = k,
              model = model,
              ...
            )
        },
        start = 0L,
        end = k_max - 1L,
        tol = tol,
        warn_if_unconverged = TRUE
      )
    },
    numeric(1)
  )
}


#' @rdname joint_life_annuities
#' @export
axy <- function(x, y, i,
                tbl = NULL, model = NULL, ...,
                k_max = 5000L, tol = 1e-12) {
  .multilife_check_basis(tbl, model)

  x <- .multilife_check_nonnegative(x, "x")
  y <- .multilife_check_nonnegative(y, "y")
  i <- .multilife_check_interest(i)

  controls <- .multilife_check_sum_controls(k_max, tol)
  k_max <- controls$k_max
  tol <- controls$tol

  values <- .multilife_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  vapply(
    seq_along(x),
    function(j) {
      if (!is.null(tbl)) {
        upper <- .multilife_joint_table_horizon(
          tbl = tbl,
          x = x[j],
          y = y[j]
        )

        if (upper < 1L) {
          return(0)
        }

        return(
          .multilife_sum_series(
            fun = function(k) {
              discount(i[j], k) *
                tpxy(
                  x = x[j],
                  y = y[j],
                  t = k,
                  tbl = tbl
                )
            },
            start = 1L,
            end = upper,
            tol = tol
          )
        )
      }

      .multilife_sum_series(
        fun = function(k) {
          discount(i[j], k) *
            tpxy(
              x = x[j],
              y = y[j],
              t = k,
              model = model,
              ...
            )
        },
        start = 1L,
        end = k_max,
        tol = tol,
        warn_if_unconverged = TRUE
      )
    },
    numeric(1)
  )
}

# -------------------------------------------------------------------------
# Joint-life insurance
# -------------------------------------------------------------------------

#' Joint-life insurance functions
#'
#' Computes temporary and whole-life insurance values for two lives under
#' the joint-life status.
#'
#' \code{Axyn1()} computes an \code{n}-year joint-life term insurance.
#'
#' \code{Axyn()} computes an \code{n}-year joint-life endowment insurance.
#'
#' \code{Axy()} computes joint-life whole-life insurance.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param n Term in years. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the survival model or
#'   life-table functions.
#' @param k_max Maximum number of terms used for an infinite series.
#' @param tol Numerical tolerance used to assess convergence.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @name joint_life_insurance
#' @export
Axyn1 <- function(x, y, n, i,
                  tbl = NULL, model = NULL, ...) {
  .multilife_check_basis(tbl, model)

  i <- .multilife_check_interest(i)

  d <- i / (1 + i)

  1 -
    d * adotxyn(
      x = x,
      y = y,
      n = n,
      i = i,
      tbl = tbl,
      model = model,
      ...
    ) -
    nExy(
      x = x,
      y = y,
      n = n,
      i = i,
      tbl = tbl,
      model = model,
      ...
    )
}

#' @rdname joint_life_insurance
#' @export
Axyn <- function(x, y, n, i,
                 tbl = NULL, model = NULL, ...) {
  Axyn1(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + nExy(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}

#' @rdname joint_life_insurance
#' @export
Axy <- function(x, y, i,
                tbl = NULL, model = NULL, ...,
                k_max = 5000L, tol = 1e-12) {
  .multilife_check_basis(tbl, model)

  i <- .multilife_check_interest(i)
  d <- i / (1 + i)

  1 - d * adotxy(
    x = x,
    y = y,
    i = i,
    tbl = tbl,
    model = model,
    ...,
    k_max = k_max,
    tol = tol
  )
}

# -------------------------------------------------------------------------
# Last-survivor annuities
# -------------------------------------------------------------------------

#' Last-survivor annuity functions
#'
#' Computes temporary and whole-life annuities for two lives under the
#' last-survivor status.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param n Term in years. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the survival model or
#'   life-table functions.
#' @param k_max Maximum number of terms used for an infinite series.
#' @param tol Numerical tolerance used to assess convergence.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @name last_survivor_annuities
#' @export
adotxybarn <- function(x, y, n, i,
                       tbl = NULL, model = NULL, ...) {
  adotxn(
    x = x,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + adotxn(
    x = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - adotxyn(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}

#' @rdname last_survivor_annuities
#' @export
axybarn <- function(x, y, n, i,
                    tbl = NULL, model = NULL, ...) {
  axn(
    x = x,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + axn(
    x = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - axyn(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}

#' @rdname last_survivor_annuities
#' @export
adotxybar <- function(x, y, i,
                      tbl = NULL, model = NULL, ...,
                      k_max = 5000L, tol = 1e-12) {
  adotx(
    x = x,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + adotx(
    x = y,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - adotxy(
    x = x,
    y = y,
    i = i,
    tbl = tbl,
    model = model,
    ...,
    k_max = k_max,
    tol = tol
  )
}

#' @rdname last_survivor_annuities
#' @export
axybar <- function(x, y, i,
                   tbl = NULL, model = NULL, ...,
                   k_max = 5000L, tol = 1e-12) {
  ax(
    x = x,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + ax(
    x = y,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - axy(
    x = x,
    y = y,
    i = i,
    tbl = tbl,
    model = model,
    ...,
    k_max = k_max,
    tol = tol
  )
}

# -------------------------------------------------------------------------
# Last-survivor insurance
# -------------------------------------------------------------------------

#' Last-survivor insurance functions
#'
#' Computes temporary and whole-life insurance values for two lives under
#' the last-survivor status.
#'
#' \code{Axybarn1()} computes temporary insurance payable at the second
#' death within the term.
#'
#' \code{Axybarn()} computes last-survivor endowment insurance.
#'
#' \code{Axybar()} computes last-survivor whole-life insurance.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param n Term in years. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the survival model or
#'   life-table functions.
#' @param k_max Maximum number of terms used for an infinite series.
#' @param tol Numerical tolerance used to assess convergence.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @name last_survivor_insurance
#' @export
Axybarn1 <- function(x, y, n, i,
                     tbl = NULL, model = NULL, ...) {
  Axn1(
    x = x,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + Axn1(
    x = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - Axyn1(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}

#' @rdname last_survivor_insurance
#' @export
Axybarn <- function(x, y, n, i,
                    tbl = NULL, model = NULL, ...) {
  Axybarn1(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + nExybar(
    x = x,
    y = y,
    n = n,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}

#' @rdname last_survivor_insurance
#' @export
Axybar <- function(x, y, i,
                   tbl = NULL, model = NULL, ...,
                   k_max = 5000L, tol = 1e-12) {
  Ax(
    x = x,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) + Ax(
    x = y,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - Axy(
    x = x,
    y = y,
    i = i,
    tbl = tbl,
    model = model,
    ...,
    k_max = k_max,
    tol = tol
  )
}

# -------------------------------------------------------------------------
# Reversionary annuities
# -------------------------------------------------------------------------

#' Reversionary annuity functions
#'
#' Computes reversionary whole-life annuities payable to one life after the
#' death of the other life.
#'
#' \code{ax_y()} computes an annuity payable to the second life after the
#' death of the first life.
#'
#' \code{ay_x()} computes an annuity payable to the first life after the
#' death of the second life.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the survival model or
#'   life-table functions.
#' @param k_max Maximum number of terms used for an infinite series.
#' @param tol Numerical tolerance used to assess convergence.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @name reversionary_annuities
#' @export
ax_y <- function(x, y, i,
                 tbl = NULL, model = NULL, ...,
                 k_max = 5000L, tol = 1e-12) {
  ax(
    x = y,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - axy(
    x = x,
    y = y,
    i = i,
    tbl = tbl,
    model = model,
    ...,
    k_max = k_max,
    tol = tol
  )
}

#' @rdname reversionary_annuities
#' @export
ay_x <- function(x, y, i,
                 tbl = NULL, model = NULL, ...,
                 k_max = 5000L, tol = 1e-12) {
  ax(
    x = x,
    i = i,
    tbl = tbl,
    model = model,
    ...
  ) - axy(
    x = x,
    y = y,
    i = i,
    tbl = tbl,
    model = model,
    ...,
    k_max = k_max,
    tol = tol
  )
}
