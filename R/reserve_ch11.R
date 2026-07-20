# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.reserve11_check_basis <- function(tbl = NULL, model = NULL) {
  if (is.null(tbl) && is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }

  if (!is.null(tbl) && !is.null(model)) {
    stop("Supply either tbl or model, not both.", call. = FALSE)
  }

  if (!is.null(tbl)) {
    .validate_life_table(tbl)
  }

  invisible(TRUE)
}

#' @noRd
.reserve11_check_nonnegative <- function(x, name) {
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
.reserve11_check_integerish <- function(x, name, positive = FALSE) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  invalid_sign <- if (positive) {
    any(x <= 0)
  } else {
    any(x < 0)
  }

  if (any(!is.finite(x)) ||
      invalid_sign ||
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
.reserve11_check_interest <- function(i, name = "i") {
  i <- as.numeric(i)

  if (length(i) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }

  if (any(!is.finite(i)) || any(i <= -1)) {
    stop(
      name,
      " must contain finite values greater than -1.",
      call. = FALSE
    )
  }

  i
}

#' @noRd
.reserve11_check_probability <- function(x, name) {
  x <- as.numeric(x)

  if (length(x) == 0L ||
      any(!is.finite(x)) ||
      any(x < 0) ||
      any(x > 1)) {
    stop(
      name,
      " must contain values in [0, 1].",
      call. = FALSE
    )
  }

  x
}

#' @noRd
.reserve11_recycle <- function(..., .names = NULL) {
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
.reserve11_check_t_less_n <- function(t, n) {
  if (any(t >= n)) {
    stop("t must satisfy t < n.", call. = FALSE)
  }

  invisible(TRUE)
}

#' @noRd
.reserve11_scalar_or_matrix <- function(x, column_names) {
  x <- as.matrix(x)
  colnames(x) <- column_names

  if (nrow(x) == 1L) {
    return(stats::setNames(as.numeric(x[1, ]), column_names))
  }

  x
}

# -------------------------------------------------------------------------
# Full preliminary term modified premiums and reserves
# -------------------------------------------------------------------------

#' Full preliminary term modified premiums and reserves
#'
#' Computes modified premiums and reserves under the full preliminary term
#' method for whole-life insurance.
#'
#' \code{alphaF()} computes the first-year modified premium
#' \eqn{\alpha^F = vq_x}.
#'
#' \code{betaF()} computes the renewal modified premium
#' \eqn{\beta^F = P_{x+1}}.
#'
#' \code{tVFx()} computes the full preliminary term reserve. The reserve is
#' zero at durations 0 and 1. For \eqn{t > 1}, it equals the net level premium
#' reserve at duration \eqn{t - 1} for a policy issued at age \eqn{x + 1}.
#'
#' @param x Issue age. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the actuarial functions.
#' @param t Nonnegative integer duration. May be scalar or vector.
#'
#' @return A numeric vector of modified premiums or reserves.
#'
#' @details
#' Under the full preliminary term method, the first policy year is treated
#' as one-year term insurance. The first-year modified premium is therefore
#' the actuarial present value of one-year term insurance at age \eqn{x}.
#'
#' Renewal premiums are based on a whole-life policy issued one year later,
#' at age \eqn{x + 1}. Accordingly, reserves after the first policy year are
#' obtained from the corresponding net level premium reserve for that
#' deferred issue age.
#'
#' @examples
#' alphaF(
#'   x = 40,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' betaF(
#'   x = 40,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' tVFx(
#'   x = 40,
#'   t = 5,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' @name full_preliminary_term
#' @rdname full_preliminary_term
#' @export
alphaF <- function(x, i, tbl = NULL, model = NULL, ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  i <- .reserve11_check_interest(i)

  values <- .reserve11_recycle(
    x,
    i,
    .names = c("x", "i")
  )

  x <- values[[1]]
  i <- values[[2]]

  Axn1(
    x = x,
    n = 1,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}


#' @rdname full_preliminary_term
#' @export
betaF <- function(x, i, tbl = NULL, model = NULL, ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  i <- .reserve11_check_interest(i)

  values <- .reserve11_recycle(
    x,
    i,
    .names = c("x", "i")
  )

  x <- values[[1]]
  i <- values[[2]]

  Px(
    x = x + 1,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}


#' @rdname full_preliminary_term
#' @export
tVFx <- function(x, t, i, tbl = NULL, model = NULL, ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  t <- .reserve11_check_integerish(t, "t")
  i <- .reserve11_check_interest(i)

  values <- .reserve11_recycle(
    x,
    t,
    i,
    .names = c("x", "t", "i")
  )

  x <- values[[1]]
  t <- values[[2]]
  i <- values[[3]]

  vapply(
    seq_along(x),
    function(j) {
      if (t[j] <= 1L) {
        return(0)
      }

      tVx(
        x = x[j] + 1,
        t = t[j] - 1L,
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )
    },
    numeric(1)
  )
}


# -------------------------------------------------------------------------
# Fractional-duration whole life reserves
# -------------------------------------------------------------------------

#' Fractional-duration whole life reserves
#'
#' Computes fractional-duration whole life reserves using linear
#' interpolation between the reserve immediately after the premium at
#' duration \eqn{t} and the reserve at duration \eqn{t+1}.
#'
#' \code{tsVx()} computes the reserve at fractional duration
#' \eqn{t+s}, where \eqn{0 \le s \le 1}.
#'
#' \code{meanVx()} computes the reserve at the midpoint of the policy
#' year (\eqn{s=0.5}).
#'
#' @param x Issue age. May be scalar or vector.
#' @param t Nonnegative integer duration. May be scalar or vector.
#' @param s Fractional duration in \eqn{[0,1]}. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the actuarial functions.
#'
#' @return A numeric vector of reserve values.
#'
#' @details
#' The fractional reserve is computed using
#'
#' \deqn{
#' {}_{t+s}V_x
#' =
#' ({}_tV_x + P_x)(1-s)
#' +
#' {}_{t+1}V_x s.
#' }
#'
#' The function \code{meanVx()} is a convenience wrapper corresponding to
#' \eqn{s=0.5}.
#'
#' @examples
#' tsVx(
#'   40,
#'   t = 10,
#'   s = 0.5,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' meanVx(
#'   40,
#'   t = 10,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' @name fractional_duration_reserves
#' @rdname fractional_duration_reserves
#' @export
tsVx <- function(x, t, s, i, tbl = NULL, model = NULL, ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  t <- .reserve11_check_integerish(t, "t")
  s <- .reserve11_check_probability(s, "s")
  i <- .reserve11_check_interest(i)

  values <- .reserve11_recycle(
    x,
    t,
    s,
    i,
    .names = c("x", "t", "s", "i")
  )

  x <- values[[1]]
  t <- values[[2]]
  s <- values[[3]]
  i <- values[[4]]

  vapply(
    seq_along(x),
    function(j) {
      reserve_t <- tVx(
        x = x[j],
        t = t[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )

      reserve_t1 <- tVx(
        x = x[j],
        t = t[j] + 1L,
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )

      premium <- Px(
        x = x[j],
        i = i[j],
        tbl = tbl,
        model = model,
        ...
      )

      (reserve_t + premium) * (1 - s[j]) +
        reserve_t1 * s[j]
    },
    numeric(1)
  )
}

#' @rdname fractional_duration_reserves
#' @export
meanVx <- function(x, t, i, tbl = NULL, model = NULL, ...) {
  tsVx(
    x = x,
    t = t,
    s = 0.5,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )
}

# -------------------------------------------------------------------------
# Fractional-duration term and endowment reserves
# -------------------------------------------------------------------------

#' Fractional-duration term and endowment reserves
#'
#' Computes fractional-duration reserves for term and endowment insurance
#' using linear interpolation:
#'
#' \deqn{
#' {}_{t+s}V =
#' ({}_tV + P)(1-s) + {}_{t+1}V s.
#' }
#'
#' \code{tsVxn()} applies the formula to endowment insurance.
#'
#' \code{tsVxn1()} applies the formula to term insurance.
#'
#' @param x Issue age. May be scalar or vector.
#' @param n Positive integer term. May be scalar or vector.
#' @param t Nonnegative integer duration satisfying \eqn{t<n}.
#' @param s Fractional duration in \eqn{[0,1]}.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the actuarial functions.
#'
#' @return A numeric vector of reserve values.
#'
#' @examples
#' tsVxn(
#'   40, n = 20, t = 10, s = 0.5, i = 0.05,
#'   model = "uniform", omega = 100
#' )
#'
#' tsVxn1(
#'   40, n = 20, t = 10, s = 0.5, i = 0.05,
#'   model = "uniform", omega = 100
#' )
#'
#' @export
tsVxn <- function(x, n, t, s, i,
                  tbl = NULL, model = NULL, ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  n <- .reserve11_check_integerish(n, "n", positive = TRUE)
  t <- .reserve11_check_integerish(t, "t")
  s <- .reserve11_check_probability(s, "s")
  i <- .reserve11_check_interest(i)

  values <- .reserve11_recycle(
    x,
    n,
    t,
    s,
    i,
    .names = c("x", "n", "t", "s", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  t <- values[[3]]
  s <- values[[4]]
  i <- values[[5]]

  .reserve11_check_t_less_n(t, n)

  vapply(seq_along(x), function(j) {
    reserve_t <- tVxn(
      x = x[j],
      n = n[j],
      t = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    reserve_t1 <- tVxn(
      x = x[j],
      n = n[j],
      t = t[j] + 1L,
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    premium <- Pxn(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    (reserve_t + premium) * (1 - s[j]) +
      reserve_t1 * s[j]
  }, numeric(1))
}

#' @rdname tsVxn
#' @export
tsVxn1 <- function(x, n, t, s, i,
                   tbl = NULL, model = NULL, ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  n <- .reserve11_check_integerish(n, "n", positive = TRUE)
  t <- .reserve11_check_integerish(t, "t")
  s <- .reserve11_check_probability(s, "s")
  i <- .reserve11_check_interest(i)

  values <- .reserve11_recycle(
    x,
    n,
    t,
    s,
    i,
    .names = c("x", "n", "t", "s", "i")
  )

  x <- values[[1]]
  n <- values[[2]]
  t <- values[[3]]
  s <- values[[4]]
  i <- values[[5]]

  .reserve11_check_t_less_n(t, n)

  vapply(seq_along(x), function(j) {
    reserve_t <- tVxn1(
      x = x[j],
      n = n[j],
      t = t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    reserve_t1 <- tVxn1(
      x = x[j],
      n = n[j],
      t = t[j] + 1L,
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    premium <- Pxn1(
      x = x[j],
      n = n[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    (reserve_t + premium) * (1 - s[j]) +
      reserve_t1 * s[j]
  }, numeric(1))
}

# -------------------------------------------------------------------------
# Gross premium and expense reserves
# -------------------------------------------------------------------------

#' Whole life gross premium and expense reserves
#'
#' Computes prospective gross premium and expense reserves for fully
#' discrete whole life insurance after issue.
#'
#' \code{tVGx()} computes the gross premium reserve:
#'
#' \deqn{
#' {}_tV_x^G =
#' (b+s)A_{x+t}
#' -
#' \left[(1-r)G-e\right]\ddot{a}_{x+t}.
#' }
#'
#' \code{tVEx()} computes the corresponding expense reserve as the gross
#' premium reserve minus the net benefit reserve.
#'
#' @param x Issue age. May be scalar or vector.
#' @param t Nonnegative integer duration. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param G Gross annual premium. May be scalar or vector.
#' @param benefit Insurance benefit amount.
#' @param renewal_premium_pct Renewal percent-of-premium expense in
#'   \eqn{[0,1]}.
#' @param renewal_policy_exp Renewal per-policy expense.
#' @param settlement_exp Settlement expense paid at death.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional parameters passed to the actuarial functions.
#'
#' @return A numeric vector of reserve values.
#'
#' @examples
#' tVGx(
#'   x = 40, t = 10, i = 0.05, G = 0.03,
#'   benefit = 1, renewal_premium_pct = 0.10,
#'   renewal_policy_exp = 0.002,
#'   settlement_exp = 0.02,
#'   model = "uniform", omega = 100
#' )
#'
#' @export
tVGx <- function(x, t, i, G,
                 benefit = 1,
                 renewal_premium_pct = 0,
                 renewal_policy_exp = 0,
                 settlement_exp = 0,
                 tbl = NULL,
                 model = NULL,
                 ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  t <- .reserve11_check_integerish(t, "t")
  i <- .reserve11_check_interest(i)
  G <- .reserve11_check_nonnegative(G, "G")
  benefit <- .reserve11_check_nonnegative(benefit, "benefit")
  renewal_premium_pct <- .reserve11_check_probability(
    renewal_premium_pct,
    "renewal_premium_pct"
  )
  renewal_policy_exp <- .reserve11_check_nonnegative(
    renewal_policy_exp,
    "renewal_policy_exp"
  )
  settlement_exp <- .reserve11_check_nonnegative(
    settlement_exp,
    "settlement_exp"
  )

  values <- .reserve11_recycle(
    x,
    t,
    i,
    G,
    benefit,
    renewal_premium_pct,
    renewal_policy_exp,
    settlement_exp,
    .names = c(
      "x",
      "t",
      "i",
      "G",
      "benefit",
      "renewal_premium_pct",
      "renewal_policy_exp",
      "settlement_exp"
    )
  )

  x <- values[[1]]
  t <- values[[2]]
  i <- values[[3]]
  G <- values[[4]]
  benefit <- values[[5]]
  renewal_premium_pct <- values[[6]]
  renewal_policy_exp <- values[[7]]
  settlement_exp <- values[[8]]

  vapply(seq_along(x), function(j) {
    insurance_value <- Ax(
      x = x[j] + t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    annuity_value <- adotx(
      x = x[j] + t[j],
      i = i[j],
      tbl = tbl,
      model = model,
      ...
    )

    (benefit[j] + settlement_exp[j]) * insurance_value -
      (
        (1 - renewal_premium_pct[j]) * G[j] -
          renewal_policy_exp[j]
      ) * annuity_value
  }, numeric(1))
}

#' @rdname tVGx
#' @export
tVEx <- function(x, t, i, G,
                 benefit = 1,
                 renewal_premium_pct = 0,
                 renewal_policy_exp = 0,
                 settlement_exp = 0,
                 tbl = NULL,
                 model = NULL,
                 ...) {
  .reserve11_check_basis(tbl, model)

  x <- .reserve11_check_nonnegative(x, "x")
  t <- .reserve11_check_integerish(t, "t")
  i <- .reserve11_check_interest(i)
  G <- .reserve11_check_nonnegative(G, "G")
  benefit <- .reserve11_check_nonnegative(benefit, "benefit")
  renewal_premium_pct <- .reserve11_check_probability(
    renewal_premium_pct,
    "renewal_premium_pct"
  )
  renewal_policy_exp <- .reserve11_check_nonnegative(
    renewal_policy_exp,
    "renewal_policy_exp"
  )
  settlement_exp <- .reserve11_check_nonnegative(
    settlement_exp,
    "settlement_exp"
  )

  values <- .reserve11_recycle(
    x,
    t,
    i,
    G,
    benefit,
    renewal_premium_pct,
    renewal_policy_exp,
    settlement_exp,
    .names = c(
      "x",
      "t",
      "i",
      "G",
      "benefit",
      "renewal_premium_pct",
      "renewal_policy_exp",
      "settlement_exp"
    )
  )

  x <- values[[1]]
  t <- values[[2]]
  i <- values[[3]]
  G <- values[[4]]
  benefit <- values[[5]]
  renewal_premium_pct <- values[[6]]
  renewal_policy_exp <- values[[7]]
  settlement_exp <- values[[8]]

  gross_reserve <- tVGx(
    x = x,
    t = t,
    i = i,
    G = G,
    benefit = benefit,
    renewal_premium_pct = renewal_premium_pct,
    renewal_policy_exp = renewal_policy_exp,
    settlement_exp = settlement_exp,
    tbl = tbl,
    model = model,
    ...
  )

  net_reserve <- tVx(
    x = x,
    t = t,
    i = i,
    tbl = tbl,
    model = model,
    ...
  )

  gross_reserve - benefit * net_reserve
}

# -------------------------------------------------------------------------
# Gross gain analysis
# -------------------------------------------------------------------------

#' Total gross gain for a discrete insurance contract
#'
#' Computes total gain during one policy year under gross premiums,
#' gross reserves, actual mortality, actual interest, and actual expenses.
#'
#' @param VtG Gross reserve at duration \code{t}.
#' @param Vt1G Gross reserve at duration \code{t + 1}.
#' @param G Gross premium.
#' @param i_actual Actual annual effective interest rate.
#' @param q_actual Actual mortality probability.
#' @param r_actual Actual percent-of-premium expense rate.
#' @param e_actual Actual per-policy expense.
#' @param s_actual Actual settlement expense.
#' @param b Benefit amount.
#'
#' @return A numeric vector of total gains.
#'
#' @examples
#' GTg_disc(
#'   VtG = 0.10,
#'   Vt1G = 0.12,
#'   G = 0.02,
#'   i_actual = 0.05,
#'   q_actual = 0.01,
#'   r_actual = 0.03,
#'   e_actual = 0,
#'   s_actual = 0.01,
#'   b = 1
#' )
#'
#' @export
GTg_disc <- function(VtG, Vt1G, G,
                     i_actual, q_actual,
                     r_actual = 0,
                     e_actual = 0,
                     s_actual = 0,
                     b = 1) {
  VtG <- .reserve11_check_nonnegative(VtG, "VtG")
  Vt1G <- .reserve11_check_nonnegative(Vt1G, "Vt1G")
  G <- .reserve11_check_nonnegative(G, "G")
  i_actual <- .reserve11_check_interest(i_actual, "i_actual")
  q_actual <- .reserve11_check_probability(q_actual, "q_actual")
  r_actual <- .reserve11_check_probability(r_actual, "r_actual")
  e_actual <- .reserve11_check_nonnegative(e_actual, "e_actual")
  s_actual <- .reserve11_check_nonnegative(s_actual, "s_actual")
  b <- .reserve11_check_nonnegative(b, "b")

  values <- .reserve11_recycle(
    VtG,
    Vt1G,
    G,
    i_actual,
    q_actual,
    r_actual,
    e_actual,
    s_actual,
    b,
    .names = c(
      "VtG",
      "Vt1G",
      "G",
      "i_actual",
      "q_actual",
      "r_actual",
      "e_actual",
      "s_actual",
      "b"
    )
  )

  VtG <- values[[1]]
  Vt1G <- values[[2]]
  G <- values[[3]]
  i_actual <- values[[4]]
  q_actual <- values[[5]]
  r_actual <- values[[6]]
  e_actual <- values[[7]]
  s_actual <- values[[8]]
  b <- values[[9]]

  (VtG + G * (1 - r_actual) - e_actual) *
    (1 + i_actual) -
    (
      (b + s_actual) * q_actual +
        (1 - q_actual) * Vt1G
    )
}

#' Ordered decomposition of gross gain
#'
#' Decomposes gross gain into interest, mortality, and expense components
#' in a user-specified order.
#'
#' For scalar input, the function returns a named numeric vector. For
#' vectorized input, it returns a numeric matrix with one row per calculation.
#'
#' @param VtG Gross reserve at duration \code{t}.
#' @param Vt1G Gross reserve at duration \code{t + 1}.
#' @param G Gross premium.
#' @param i_assumed Assumed annual effective interest rate.
#' @param q_assumed Assumed mortality probability.
#' @param r_assumed Assumed percent-of-premium expense rate.
#' @param e_assumed Assumed per-policy expense.
#' @param s_assumed Assumed settlement expense.
#' @param i_actual Actual annual effective interest rate.
#' @param q_actual Actual mortality probability.
#' @param r_actual Actual percent-of-premium expense rate.
#' @param e_actual Actual per-policy expense.
#' @param s_actual Actual settlement expense.
#' @param b Benefit amount.
#' @param order Character vector containing \code{"interest"},
#'   \code{"mortality"}, and \code{"expense"} exactly once.
#'
#' @return For scalar input, a named numeric vector. For vectorized input,
#'   a numeric matrix with columns \code{total_gain}, \code{interest},
#'   \code{mortality}, \code{expense}, and \code{check}.
#'
#' @examples
#' decompGg_disc(
#'   VtG = 3950.73,
#'   Vt1G = 4607.07,
#'   G = 685,
#'   i_assumed = 0.06,
#'   q_assumed = 0.00592,
#'   r_assumed = 0.05,
#'   e_assumed = 0,
#'   s_assumed = 300,
#'   i_actual = 0.065,
#'   q_actual = 0.005,
#'   r_actual = 0.06,
#'   e_actual = 0,
#'   s_actual = 100,
#'   b = 50000
#' )
#'
#' @export
decompGg_disc <- function(
    VtG,
    Vt1G,
    G,
    i_assumed,
    q_assumed,
    r_assumed = 0,
    e_assumed = 0,
    s_assumed = 0,
    i_actual,
    q_actual,
    r_actual = 0,
    e_actual = 0,
    s_actual = 0,
    b = 1,
    order = c("interest", "mortality", "expense")) {

  allowed <- c("interest", "mortality", "expense")

  if (!is.character(order) ||
      length(order) != 3L ||
      length(unique(order)) != 3L ||
      !setequal(order, allowed)) {
    stop(
      "order must contain 'interest', 'mortality', and 'expense' exactly once.",
      call. = FALSE
    )
  }

  VtG <- .reserve11_check_nonnegative(VtG, "VtG")
  Vt1G <- .reserve11_check_nonnegative(Vt1G, "Vt1G")
  G <- .reserve11_check_nonnegative(G, "G")

  i_assumed <- .reserve11_check_interest(i_assumed, "i_assumed")
  q_assumed <- .reserve11_check_probability(q_assumed, "q_assumed")
  r_assumed <- .reserve11_check_probability(r_assumed, "r_assumed")
  e_assumed <- .reserve11_check_nonnegative(e_assumed, "e_assumed")
  s_assumed <- .reserve11_check_nonnegative(s_assumed, "s_assumed")

  i_actual <- .reserve11_check_interest(i_actual, "i_actual")
  q_actual <- .reserve11_check_probability(q_actual, "q_actual")
  r_actual <- .reserve11_check_probability(r_actual, "r_actual")
  e_actual <- .reserve11_check_nonnegative(e_actual, "e_actual")
  s_actual <- .reserve11_check_nonnegative(s_actual, "s_actual")

  b <- .reserve11_check_nonnegative(b, "b")

  values <- .reserve11_recycle(
    VtG,
    Vt1G,
    G,
    i_assumed,
    q_assumed,
    r_assumed,
    e_assumed,
    s_assumed,
    i_actual,
    q_actual,
    r_actual,
    e_actual,
    s_actual,
    b,
    .names = c(
      "VtG",
      "Vt1G",
      "G",
      "i_assumed",
      "q_assumed",
      "r_assumed",
      "e_assumed",
      "s_assumed",
      "i_actual",
      "q_actual",
      "r_actual",
      "e_actual",
      "s_actual",
      "b"
    )
  )

  VtG <- values[[1]]
  Vt1G <- values[[2]]
  G <- values[[3]]
  i_assumed <- values[[4]]
  q_assumed <- values[[5]]
  r_assumed <- values[[6]]
  e_assumed <- values[[7]]
  s_assumed <- values[[8]]
  i_actual <- values[[9]]
  q_actual <- values[[10]]
  r_actual <- values[[11]]
  e_actual <- values[[12]]
  s_actual <- values[[13]]
  b <- values[[14]]

  calculate_one <- function(j) {
    assumed_state <- list(
      i = i_assumed[j],
      q = q_assumed[j],
      r = r_assumed[j],
      e = e_assumed[j],
      s = s_assumed[j]
    )

    actual_state <- list(
      i = i_actual[j],
      q = q_actual[j],
      r = r_actual[j],
      e = e_actual[j],
      s = s_actual[j]
    )

    evaluate_gain <- function(state) {
      as.numeric(GTg_disc(
        VtG = VtG[j],
        Vt1G = Vt1G[j],
        G = G[j],
        i_actual = state$i,
        q_actual = state$q,
        r_actual = state$r,
        e_actual = state$e,
        s_actual = state$s,
        b = b[j]
      ))
    }

    components <- c(
      interest = 0,
      mortality = 0,
      expense = 0
    )

    current_state <- assumed_state
    running_gain <- 0
    total_gain <- evaluate_gain(actual_state)

    # Compute the first two components sequentially.
    for (k in seq_len(2L)) {
      component <- order[k]

      if (component == "interest") {
        current_state$i <- actual_state$i
      } else if (component == "mortality") {
        current_state$q <- actual_state$q
      } else {
        current_state$r <- actual_state$r
        current_state$e <- actual_state$e
        current_state$s <- actual_state$s
      }

      new_gain <- evaluate_gain(current_state)

      components[component] <- new_gain - running_gain
      running_gain <- new_gain
    }

    # Use the final component as the balancing item so that the
    # decomposition sums exactly to the directly calculated total gain.
    final_component <- order[3L]

    components[final_component] <-
      total_gain - sum(components[order[1:2]])

    c(
      total_gain = total_gain,
      interest = unname(components["interest"]),
      mortality = unname(components["mortality"]),
      expense = unname(components["expense"]),
      check = sum(components)
    )
  }

  rows <- lapply(seq_along(VtG), calculate_one)

  if (length(rows) == 1L) {
    return(rows[[1L]])
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
