# reserve_extended.R
#
# Retrospective reserves, deferred-insurance reserves, annuity reserves,
# continuous-style gain and loss calculations, and numerical helpers for
# Thiele's differential equation.

# -------------------------------------------------------------------------
# Internal validation and recycling helpers
# -------------------------------------------------------------------------

#' @noRd
.reserve_ext_check_finite <- function(x, name) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop(name, " must contain finite numeric values.", call. = FALSE)
  }

  x
}

#' @noRd
.reserve_ext_check_nonnegative <- function(x, name) {
  x <- .reserve_ext_check_finite(x, name)

  if (any(x < 0)) {
    stop(name, " must contain nonnegative finite values.", call. = FALSE)
  }

  x
}

#' @noRd
.reserve_ext_check_positive <- function(x, name) {
  x <- .reserve_ext_check_finite(x, name)

  if (any(x <= 0)) {
    stop(name, " must contain positive finite values.", call. = FALSE)
  }

  x
}

#' @noRd
.reserve_ext_check_integerish <- function(x, name, positive = FALSE) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  invalid_numeric <- any(!is.finite(x))
  invalid_sign <- if (positive) any(x <= 0) else any(x < 0)
  invalid_integer <- any(abs(x - round(x)) > 1e-10)

  if (invalid_numeric || invalid_sign || invalid_integer) {
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
.reserve_ext_check_probability <- function(x, name) {
  x <- .reserve_ext_check_finite(x, name)

  if (any(x < 0 | x > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }

  x
}

#' @noRd
.reserve_ext_check_interest <- function(i, name = "i") {
  i <- .reserve_ext_check_finite(i, name)

  if (any(i <= -1)) {
    stop(name, " must contain values greater than -1.", call. = FALSE)
  }

  i
}

#' @noRd
.reserve_ext_check_basis <- function(tbl, model) {
  if (is.null(tbl) && is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }

  if (!is.null(tbl) && !is.null(model)) {
    stop("Supply only one of tbl or model.", call. = FALSE)
  }

  invisible(TRUE)
}

#' @noRd
.reserve_ext_recycle <- function(..., .names = NULL) {
  values <- list(...)
  lengths <- vapply(values, length, integer(1))

  if (length(values) == 0L || any(lengths == 0L)) {
    stop("Inputs must have positive length.", call. = FALSE)
  }

  common_length <- max(lengths)

  if (is.null(.names)) {
    .names <- rep("Input", length(values))
  }

  if (length(.names) != length(values)) {
    stop("Internal error: incorrect number of argument names.", call. = FALSE)
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
  }

  lapply(values, rep, length.out = common_length)
}

#' @noRd
.reserve_ext_recycle_to <- function(x, size, name) {
  if (!length(x) %in% c(1L, size)) {
    stop(
      name,
      " must have length 1 or ",
      size,
      ".",
      call. = FALSE
    )
  }

  rep(x, length.out = size)
}

#' @noRd
.reserve_ext_check_t_le_n <- function(t, n, allow_equal = TRUE) {
  invalid <- if (allow_equal) t > n else t >= n

  if (any(invalid)) {
    relation <- if (allow_equal) "t <= n" else "t < n"
    stop("t must satisfy ", relation, ".", call. = FALSE)
  }

  invisible(TRUE)
}

#' @noRd
.reserve_ext_check_positive_denominator <- function(x, name) {
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop(name, " must be positive and finite.", call. = FALSE)
  }

  x
}


# Internal compatibility helpers
#
# These wrappers preserve internal helper names used by reserve and
# gain/loss functions in existing package source files. They are not
# exported and do not create documentation pages.

#' @noRd
.check_nonneg_numeric <- function(x, name) {
  .reserve_ext_check_nonnegative(x, name)
}

#' @noRd
.check_nonneg_integerish <- function(x, name) {
  .reserve_ext_check_integerish(
    x = x,
    name = name,
    positive = FALSE
  )
}

#' @noRd
.check_prob_ch10 <- function(x, name) {
  .reserve_ext_check_probability(x, name)
}

#' @noRd
.check_interest_ch10 <- function(i, name = "i") {
  .reserve_ext_check_interest(i, name)
}

#' @noRd
.common_length_ch10 <- function(...) {
  lengths <- vapply(list(...), length, integer(1))

  if (length(lengths) == 0L) {
    return(0L)
  }

  max(lengths)
}

#' @noRd
.recycle_common_ch10 <- function(...) {
  values <- list(...)
  argument_names <- names(values)

  if (is.null(argument_names) ||
      any(!nzchar(argument_names))) {
    argument_names <- rep("Input", length(values))
  }

  .reserve_ext_recycle(
    ...,
    .names = argument_names
  )
}

#' @noRd
.recycle2_ch10 <- function(a,
                           b,
                           a_name = "a",
                           b_name = "b") {
  values <- .reserve_ext_recycle(
    a,
    b,
    .names = c(a_name, b_name)
  )

  list(
    a = values[[1]],
    b = values[[2]]
  )
}

#' @noRd
.recycle3_ch10 <- function(a,
                           b,
                           c,
                           a_name = "a",
                           b_name = "b",
                           c_name = "c") {
  values <- .reserve_ext_recycle(
    a,
    b,
    c,
    .names = c(a_name, b_name, c_name)
  )

  list(
    a = values[[1]],
    b = values[[2]],
    c = values[[3]]
  )
}

#' @noRd
.check_t_le_n <- function(t, n, allow_equal = TRUE) {
  .reserve_ext_check_t_le_n(
    t = t,
    n = n,
    allow_equal = allow_equal
  )
}


# -------------------------------------------------------------------------
# Retrospective reserves
# -------------------------------------------------------------------------

#' Retrospective whole life reserve
#'
#' Computes the retrospective net level premium reserve
#' \deqn{{}_tV_x =
#' P_x\ddot{s}_{x:\overline{t}|}
#' -
#' \frac{A_{x:\overline{t}|}^{1}}{{}_tE_x}.}
#'
#' The mortality basis may be supplied through either a life table or a
#' parametric survival model.
#'
#' @param x Issue age. May be scalar or vector.
#' @param t Nonnegative integer duration. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the actuarial functions.
#' @param tbl Optional life table object. Supply by name.
#'
#' @return Numeric vector of retrospective reserve values.
#' @examples
#' tVx_ret(
#'   40,
#'   t = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#' @export
tVx_ret <- function(x, t, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  t <- .reserve_ext_check_integerish(t, "t")
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    t,
    i,
    .names = c("x", "t", "i")
  )

  x <- values[[1]]
  t <- values[[2]]
  i <- values[[3]]

  vapply(seq_along(x), function(j) {
    premium <- Px(
      x = x[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    accumulated_premiums <- sdotxn(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    accumulated_claims <- Axn1(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    ) / nEx(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    premium * accumulated_premiums - accumulated_claims
  }, numeric(1))
}

#' Retrospective endowment insurance reserve
#'
#' Computes the retrospective reserve for an endowment insurance at a
#' duration satisfying \eqn{0 \le t \le n}. At maturity, the reserve is 1.
#'
#' @inheritParams tVx_ret
#' @param n Nonnegative integer contract term. May be scalar or vector.
#'
#' @return Numeric vector of retrospective reserve values.
#' @examples
#' tVxn_ret(
#'   40,
#'   n = 20,
#'   t = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#' @export
tVxn_ret <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n")
  t <- .reserve_ext_check_integerish(t, "t")
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    t,
    i,
    .names = c("x", "n", "t", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  t <- values[[3]]
  i <- values[[4]]

  .reserve_ext_check_t_le_n(t, n)

  vapply(seq_along(x), function(j) {
    if (t[j] == n[j]) {
      return(1)
    }

    premium <- Pxn(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    accumulated_premiums <- sdotxn(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    accumulated_claims <- Axn1(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    ) / nEx(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    premium * accumulated_premiums - accumulated_claims
  }, numeric(1))
}

#' Retrospective term insurance reserve
#'
#' Computes the retrospective reserve for a term insurance at a duration
#' satisfying \eqn{0 \le t \le n}. At expiry, the reserve is 0.
#'
#' @inheritParams tVxn_ret
#'
#' @return Numeric vector of retrospective reserve values.
#' @examples
#' tVxn1_ret(
#'   40,
#'   n = 20,
#'   t = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#' @export
tVxn1_ret <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n")
  t <- .reserve_ext_check_integerish(t, "t")
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    t,
    i,
    .names = c("x", "n", "t", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  t <- values[[3]]
  i <- values[[4]]

  .reserve_ext_check_t_le_n(t, n)

  vapply(seq_along(x), function(j) {
    if (t[j] == n[j]) {
      return(0)
    }

    premium <- Pxn1(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    accumulated_premiums <- sdotxn(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    accumulated_claims <- Axn1(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    ) / nEx(
      x = x[j],
      n = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    premium * accumulated_premiums - accumulated_claims
  }, numeric(1))
}


# -------------------------------------------------------------------------
# Deferred insurance reserves
# -------------------------------------------------------------------------

#' Deferred insurance reserves
#'
#' Computes prospective reserves for deferred whole life insurance contracts.
#'
#' \code{tVnAx()} computes the reserve at duration \code{t} for an
#' \code{n}-year deferred whole life insurance funded by level annual
#' premiums during the deferral period.
#'
#' \code{htVnAx()} computes the corresponding reserve when premiums are
#' limited to the first \code{h} years, where \code{h <= n}.
#'
#' @param x Issue age.
#' @param n Deferral period in years.
#' @param t Duration in years.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model name.
#' @param ... Additional arguments passed to the underlying actuarial
#'   functions.
#'
#' @return A numeric vector of prospective reserve values.
#'
#' @examples
#' tVnAx(
#'   x = 40, n = 20, t = 10, i = 0.05,
#'   model = "uniform", omega = 100
#' )
#'
#' htVnAx(
#'   x = 40, n = 20, h = 10, t = 5, i = 0.05,
#'   model = "uniform", omega = 100
#' )
#'
#' @name deferred_insurance_reserves
#' @aliases tVnAx htVnAx
#' @export
tVnAx <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n")
  t <- .reserve_ext_check_integerish(t, "t")
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    t,
    i,
    .names = c("x", "n", "t", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  t <- values[[3]]
  i <- values[[4]]

  vapply(seq_along(x), function(j) {
    if (t[j] < n[j]) {
      future_benefit <- nAx(
        x = x[j] + t[j],
        n = n[j] - t[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )

      future_premiums <- PnAx(
        x = x[j],
        n = n[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      ) * adotxn(
        x = x[j] + t[j],
        n = n[j] - t[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )

      return(future_benefit - future_premiums)
    }

    Ax(
      x = x[j] + t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )
  }, numeric(1))
}


#' @rdname deferred_insurance_reserves
#' @param h Premium-paying period in years.
#' @export
htVnAx <- function(x, n, h, t, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n")
  h <- .reserve_ext_check_integerish(h, "h", positive = TRUE)
  t <- .reserve_ext_check_integerish(t, "t")
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    h,
    t,
    i,
    .names = c("x", "n", "h", "t", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  h <- values[[3]]
  t <- values[[4]]
  i <- values[[5]]

  if (any(h > n)) {
    stop("h must satisfy h <= n.", call. = FALSE)
  }

  vapply(seq_along(x), function(j) {
    if (t[j] < h[j]) {
      future_benefit <- nAx(
        x = x[j] + t[j],
        n = n[j] - t[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )

      future_premiums <- tPnAx(
        x = x[j],
        n = n[j],
        t = h[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      ) * adotxn(
        x = x[j] + t[j],
        n = h[j] - t[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )

      return(future_benefit - future_premiums)
    }

    if (t[j] < n[j]) {
      return(nAx(
        x = x[j] + t[j],
        n = n[j] - t[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      ))
    }

    Ax(
      x = x[j] + t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )
  }, numeric(1))
}

# -------------------------------------------------------------------------
# Annuity reserves
# -------------------------------------------------------------------------

#' Net premium for a deferred annuity-due
#'
#' Computes
#' \deqn{
#' P({}_{n|}\ddot{a}_x)
#' =
#' \frac{{}_{n|}\ddot{a}_x}
#' {\ddot{a}_{x:\overline{n}|}}.
#' }
#'
#' @param x Issue age. May be scalar or vector.
#' @param n Positive integer deferral period. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the actuarial functions.
#' @param tbl Optional life table object. Supply by name.
#'
#' @return Numeric vector of net annual premiums.
#' @examples
#' PnAdotx(
#'   40,
#'   n = 20,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#' @export
PnAdotx <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n", positive = TRUE)
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    i,
    .names = c("x", "n", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  i <- values[[3]]

  vapply(seq_along(x), function(j) {
    denominator <- adotxn(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    .reserve_ext_check_positive_denominator(
      denominator,
      "The premium annuity denominator"
    )

    nadotx(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    ) / denominator
  }, numeric(1))
}

#' Reserve for a deferred annuity-due
#'
#' Computes the reserve during the deferral period, where \eqn{0 \le t < n}.
#'
#' @inheritParams PnAdotx
#' @param t Nonnegative integer duration. May be scalar or vector.
#'
#' @return Numeric vector of reserve values.
#' @examples
#' tVnAdotx(
#'   40,
#'   n = 20,
#'   t = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#' @export
tVnAdotx <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n", positive = TRUE)
  t <- .reserve_ext_check_integerish(t, "t")
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    t,
    i,
    .names = c("x", "n", "t", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  t <- values[[3]]
  i <- values[[4]]

  .reserve_ext_check_t_le_n(t, n, allow_equal = FALSE)

  vapply(seq_along(x), function(j) {
    future_benefits <- nadotx(
      x = x[j] + t[j],
      n = n[j] - t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    future_premiums <- PnAdotx(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    ) * adotxn(
      x = x[j] + t[j],
      n = n[j] - t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    future_benefits - future_premiums
  }, numeric(1))
}

#' Net premium for a deferred annuity-immediate
#'
#' Computes
#' \deqn{
#' P({}_{n|}a_x)
#' =
#' \frac{{}_{n|}a_x}
#' {\ddot{a}_{x:\overline{n}|}}.
#' }
#'
#' @inheritParams PnAdotx
#'
#' @return Numeric vector of net annual premiums.
#' @examples
#' Pnax(
#'   40,
#'   n = 20,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#' @export
Pnax <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n", positive = TRUE)
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    i,
    .names = c("x", "n", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  i <- values[[3]]

  vapply(seq_along(x), function(j) {
    denominator <- adotxn(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    .reserve_ext_check_positive_denominator(
      denominator,
      "The premium annuity denominator"
    )

    nax(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    ) / denominator
  }, numeric(1))
}

#' Reserve for a deferred annuity-immediate
#'
#' Computes the reserve during the deferral period, where \eqn{0 \le t < n}.
#'
#' @inheritParams tVnAdotx
#'
#' @return Numeric vector of reserve values.
#' @examples
#' tVnax(
#'   40,
#'   n = 20,
#'   t = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#' @export
tVnax <- function(x, n, t, i, model = NULL, ..., tbl = NULL) {
  .reserve_ext_check_basis(tbl, model)

  x <- .reserve_ext_check_nonnegative(x, "x")
  n <- .reserve_ext_check_integerish(n, "n", positive = TRUE)
  t <- .reserve_ext_check_integerish(t, "t")
  i <- .reserve_ext_check_interest(i)

  values <- .reserve_ext_recycle(
    x,
    n,
    t,
    i,
    .names = c("x", "n", "t", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  t <- values[[3]]
  i <- values[[4]]

  .reserve_ext_check_t_le_n(t, n, allow_equal = FALSE)

  vapply(seq_along(x), function(j) {
    future_benefits <- nax(
      x = x[j] + t[j],
      n = n[j] - t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    future_premiums <- Pnax(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    ) * adotxn(
      x = x[j] + t[j],
      n = n[j] - t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    future_benefits - future_premiums
  }, numeric(1))
}

# -------------------------------------------------------------------------
# Continuous-style gain and loss functions
# -------------------------------------------------------------------------

#' Total gain for a continuous-style one-step recursion
#'
#' Computes the amount accumulated during a step, less the expected reserve
#' required at the end of the step.
#'
#' Reserves, premium rates, benefits, and forces of interest may be negative
#' when such values are meaningful for the application. The step length must
#' be positive, and the survival probability must lie in \eqn{[0,1]}.
#'
#' @param Vt Reserve at time `t`.
#' @param Vt1 Reserve at time `t + h`.
#' @param P Premium rate.
#' @param delta_actual Actual force of interest.
#' @param p_actual Actual survival probability over the step.
#' @param benefit Benefit paid at the start of the step.
#' @param h Positive step length.
#'
#' @return Numeric vector of gain values.
#' @examples
#' GT_cont(
#'   Vt = 10,
#'   Vt1 = 11,
#'   P = 1,
#'   delta_actual = 0.05,
#'   p_actual = 0.99
#' )
#' @export
GT_cont <- function(Vt,
                    Vt1,
                    P,
                    delta_actual,
                    p_actual,
                    benefit = 0,
                    h = 1) {
  Vt <- .reserve_ext_check_finite(Vt, "Vt")
  Vt1 <- .reserve_ext_check_finite(Vt1, "Vt1")
  P <- .reserve_ext_check_finite(P, "P")
  delta_actual <- .reserve_ext_check_finite(
    delta_actual,
    "delta_actual"
  )
  p_actual <- .reserve_ext_check_probability(p_actual, "p_actual")
  benefit <- .reserve_ext_check_finite(benefit, "benefit")
  h <- .reserve_ext_check_positive(h, "h")

  values <- .reserve_ext_recycle(
    Vt,
    Vt1,
    P,
    delta_actual,
    p_actual,
    benefit,
    h,
    .names = c(
      "Vt",
      "Vt1",
      "P",
      "delta_actual",
      "p_actual",
      "benefit",
      "h"
    )
  )

  Vt <- values[[1]]
  Vt1 <- values[[2]]
  P <- values[[3]]
  delta_actual <- values[[4]]
  p_actual <- values[[5]]
  benefit <- values[[6]]
  h <- values[[7]]

  (Vt - benefit + P * h) * exp(delta_actual * h) -
    p_actual * Vt1
}

#' Mortality gain for a continuous-style recursion
#'
#' Evaluates the one-step gain using the assumed force of interest and the
#' actual survival probability.
#'
#' @inheritParams GT_cont
#' @param delta_assumed Assumed force of interest.
#'
#' @return Numeric vector of mortality gain values.
#' @examples
#' GM_cont(
#'   Vt = 10,
#'   Vt1 = 11,
#'   P = 1,
#'   delta_assumed = 0.05,
#'   p_actual = 0.99
#' )
#' @export
GM_cont <- function(Vt,
                    Vt1,
                    P,
                    delta_assumed,
                    p_actual,
                    benefit = 0,
                    h = 1) {
  GT_cont(
    Vt = Vt,
    Vt1 = Vt1,
    P = P,
    delta_actual = delta_assumed,
    p_actual = p_actual,
    benefit = benefit,
    h = h
  )
}

#' Interest gain for a continuous-style recursion
#'
#' Evaluates the one-step gain using the actual force of interest and the
#' assumed survival probability.
#'
#' @inheritParams GT_cont
#' @param p_assumed Assumed survival probability over the step.
#'
#' @return Numeric vector of interest gain values.
#' @examples
#' GI_cont(
#'   Vt = 10,
#'   Vt1 = 11,
#'   P = 1,
#'   delta_actual = 0.05,
#'   p_assumed = 0.99
#' )
#' @export
GI_cont <- function(Vt,
                    Vt1,
                    P,
                    delta_actual,
                    p_assumed,
                    benefit = 0,
                    h = 1) {
  GT_cont(
    Vt = Vt,
    Vt1 = Vt1,
    P = P,
    delta_actual = delta_actual,
    p_actual = p_assumed,
    benefit = benefit,
    h = h
  )
}

# -------------------------------------------------------------------------
# Thiele equation helpers
# -------------------------------------------------------------------------

#' One backward numerical step for Thiele's equation
#'
#' Approximates the reserve at time `t` from a known reserve at time
#' `t + h`.
#'
#' The implemented step is
#' \deqn{
#' V_t =
#' \frac{
#' V_{t+h} - hP + h\mu B
#' }{
#' 1 + h(\delta + \mu)
#' }.
#' }
#'
#' @param V_next Reserve at time `t + h`.
#' @param P Premium rate.
#' @param delta Force of interest.
#' @param mu Nonnegative force of mortality at time `t`.
#' @param benefit Benefit amount.
#' @param h Positive step size.
#'
#' @return Numeric vector of reserve approximations.
#' @examples
#' thiele_backward_step(
#'   V_next = 1000,
#'   P = 26.96,
#'   delta = 0.058,
#'   mu = 0.002,
#'   benefit = 1000,
#'   h = 1
#' )
#' @export
thiele_backward_step <- function(V_next,
                                 P,
                                 delta,
                                 mu,
                                 benefit = 1,
                                 h = 1) {
  V_next <- .reserve_ext_check_finite(V_next, "V_next")
  P <- .reserve_ext_check_finite(P, "P")
  delta <- .reserve_ext_check_finite(delta, "delta")
  mu <- .reserve_ext_check_nonnegative(mu, "mu")
  benefit <- .reserve_ext_check_finite(benefit, "benefit")
  h <- .reserve_ext_check_positive(h, "h")

  values <- .reserve_ext_recycle(
    V_next,
    P,
    delta,
    mu,
    benefit,
    h,
    .names = c(
      "V_next",
      "P",
      "delta",
      "mu",
      "benefit",
      "h"
    )
  )

  V_next <- values[[1]]
  P <- values[[2]]
  delta <- values[[3]]
  mu <- values[[4]]
  benefit <- values[[5]]
  h <- values[[6]]

  denominator <- 1 + h * (delta + mu)

  if (any(abs(denominator) < .Machine$double.eps^0.5)) {
    stop(
      "The backward-step denominator is zero or numerically indistinguishable from zero.",
      call. = FALSE
    )
  }

  (V_next - h * P + h * mu * benefit) / denominator
}

#' Reserve derivative from Thiele's equation
#'
#' Computes
#' \deqn{
#' \frac{dV}{dt}
#' =
#' P + \delta V - \mu(B - V).
#' }
#'
#' @param V Reserve at time `t`.
#' @param P Premium rate.
#' @param delta Force of interest.
#' @param mu Nonnegative force of mortality.
#' @param benefit Benefit amount.
#'
#' @return Numeric vector of reserve derivatives.
#' @examples
#' thiele_dVdt(
#'   V = 900,
#'   P = 25,
#'   delta = 0.05,
#'   mu = 0.002,
#'   benefit = 1000
#' )
#' @export
thiele_dVdt <- function(V, P, delta, mu, benefit = 1) {
  V <- .reserve_ext_check_finite(V, "V")
  P <- .reserve_ext_check_finite(P, "P")
  delta <- .reserve_ext_check_finite(delta, "delta")
  mu <- .reserve_ext_check_nonnegative(mu, "mu")
  benefit <- .reserve_ext_check_finite(benefit, "benefit")

  values <- .reserve_ext_recycle(
    V,
    P,
    delta,
    mu,
    benefit,
    .names = c("V", "P", "delta", "mu", "benefit")
  )

  V <- values[[1]]
  P <- values[[2]]
  delta <- values[[3]]
  mu <- values[[4]]
  benefit <- values[[5]]

  P + delta * V - mu * (benefit - V)
}

#' Backward reserve path from a terminal value
#'
#' Starting from a terminal reserve at the final time, computes reserves
#' backward over a strictly increasing time grid using
#' `thiele_backward_step()`.
#'
#' @param times Finite numeric vector of strictly increasing times.
#' @param V_terminal Finite scalar reserve at the final time.
#' @param P Premium rate, scalar or vector with one value per time step.
#' @param delta Force of interest, scalar or vector with one value per step.
#' @param mu Nonnegative force of mortality, scalar or vector with one value
#'   per step.
#' @param benefit Benefit amount, scalar or vector with one value per step.
#'
#' @return Numeric vector of reserve values corresponding to `times`.
#' @examples
#' times <- seq(19, 20, by = 0.25)
#'
#' thiele_backward_path(
#'   times,
#'   V_terminal = 1000,
#'   P = 26.96,
#'   delta = 0.058,
#'   mu = 0.002,
#'   benefit = 1000
#' )
#' @export
thiele_backward_path <- function(times,
                                 V_terminal,
                                 P,
                                 delta,
                                 mu,
                                 benefit = 1) {
  times <- .reserve_ext_check_finite(times, "times")

  if (length(times) < 2L) {
    stop("times must have length at least 2.", call. = FALSE)
  }

  if (any(diff(times) <= 0)) {
    stop("times must be strictly increasing.", call. = FALSE)
  }

  V_terminal <- .reserve_ext_check_finite(
    V_terminal,
    "V_terminal"
  )

  if (length(V_terminal) != 1L) {
    stop("V_terminal must be a finite scalar.", call. = FALSE)
  }

  P <- .reserve_ext_check_finite(P, "P")
  delta <- .reserve_ext_check_finite(delta, "delta")
  mu <- .reserve_ext_check_nonnegative(mu, "mu")
  benefit <- .reserve_ext_check_finite(benefit, "benefit")

  n_steps <- length(times) - 1L

  P <- .reserve_ext_recycle_to(P, n_steps, "P")
  delta <- .reserve_ext_recycle_to(delta, n_steps, "delta")
  mu <- .reserve_ext_recycle_to(mu, n_steps, "mu")
  benefit <- .reserve_ext_recycle_to(
    benefit,
    n_steps,
    "benefit"
  )

  reserves <- numeric(length(times))
  reserves[length(times)] <- V_terminal

  for (k in seq.int(n_steps, 1L, by = -1L)) {
    h <- times[k + 1L] - times[k]

    reserves[k] <- thiele_backward_step(
      V_next = reserves[k + 1L],
      P = P[k],
      delta = delta[k],
      mu = mu[k],
      benefit = benefit[k],
      h = h
    )
  }

  reserves
}
