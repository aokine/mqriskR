#' Reserve functions for Chapter 11
#'
#' Chapter 11 functions for modified reserves, fractional-duration reserves,
#' gross premium reserves, expense reserves, and gross gain analysis.
#'
#' @name reserve_ch11
NULL

# -------------------------------------------------------------------------
# Full preliminary term (FPT) modified reserves
# -------------------------------------------------------------------------

#' Full preliminary term first-year modified premium
#'
#' Computes
#' \eqn{\alpha^F = v q_x = A_{x:\overline{1}|}^{1}}.
#'
#' @param x Issue age.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' alphaF(40, i = 0.05, model = "uniform", omega = 100)
#' @export
alphaF <- function(x, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  Axn1(x = x, n = 1, i = i, tbl = tbl, model = model, ...)
}

#' Full preliminary term renewal modified premium
#'
#' Computes the FPT renewal premium
#' \eqn{\beta^F = P_{x+1}}
#' for whole life insurance.
#'
#' @inheritParams alphaF
#'
#' @return Numeric vector.
#' @examples
#' betaF(40, i = 0.05, model = "uniform", omega = 100)
#' @export
betaF <- function(x, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  Px(x = x + 1, i = i, tbl = tbl, model = model, ...)
}

#' Full preliminary term reserve for whole life insurance
#'
#' Computes the FPT reserve for a whole life insurance.
#' For whole life insurance,
#' \eqn{{}_1V^F = 0}
#' and for \eqn{t \ge 1},
#' \eqn{{}_tV^F = {}_{t-1}V_{x+1}^{NLP}}.
#'
#' @param x Issue age.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tVFx(40, t = 5, i = 0.05, model = "uniform", omega = 100)
#' @export
tVFx <- function(x, t, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")

  xt <- .recycle2_ch10(x, t, "x", "t")
  x <- xt$a
  t <- xt$b

  out <- numeric(length(x))

  for (j in seq_along(x)) {
    if (t[j] <= 1L) {
      out[j] <- 0
    } else {
      out[j] <- tVx(
        x = x[j] + 1,
        t = t[j] - 1,
        i = i,
        tbl = tbl,
        model = model,
        ...
      )
    }
  }

  out
}

# -------------------------------------------------------------------------
# Fractional-duration reserves
# -------------------------------------------------------------------------

#' Fractional-duration whole life reserve
#'
#' Computes the interpolated reserve
#' \eqn{{}_{t+s}V = ({}_tV + P_x)(1-s) + {}_{t+1}V \cdot s}
#' for \eqn{0 \le s \le 1}.
#'
#' @param x Issue age.
#' @param t Integer contract duration.
#' @param s Fractional part of duration in [0, 1].
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tsVx(40, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100)
#' @export
tsVx <- function(x, t, s, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")
  s <- as.numeric(s)

  if (length(s) == 0L || any(!is.finite(s)) || any(s < 0) || any(s > 1)) {
    stop("s must contain values in [0, 1].", call. = FALSE)
  }

  xts <- .recycle3_ch10(x, t, s, "x", "t", "s")
  x <- xts$a
  t <- xts$b
  s <- xts$c

  out <- numeric(length(x))

  for (j in seq_along(x)) {
    Vt <- tVx(x = x[j], t = t[j], i = i, tbl = tbl, model = model, ...)
    Vt1 <- tVx(x = x[j], t = t[j] + 1, i = i, tbl = tbl, model = model, ...)
    P <- Px(x = x[j], i = i, tbl = tbl, model = model, ...)

    out[j] <- (Vt + P) * (1 - s[j]) + Vt1 * s[j]
  }

  out
}

#' Mean reserve for whole life insurance
#'
#' Computes the mean reserve
#' \eqn{{}_{t+1/2}V}.
#'
#' @inheritParams tsVx
#'
#' @return Numeric vector.
#' @examples
#' meanVx(40, t = 10, i = 0.05, model = "uniform", omega = 100)
#' @export
meanVx <- function(x, t, i, tbl = NULL, model = NULL, ...) {
  tsVx(x = x, t = t, s = 0.5, i = i, tbl = tbl, model = model, ...)
}

#' Fractional-duration endowment reserve
#'
#' Computes the interpolated reserve
#' \eqn{{}_{t+s}V = ({}_tV + P)(1-s) + {}_{t+1}V \cdot s}
#' for an n-year endowment insurance.
#'
#' @param x Issue age.
#' @param n Term in years.
#' @param t Integer duration with t < n.
#' @param s Fractional part in [0, 1].
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tsVxn(40, n = 20, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100)
#' @export
tsVxn <- function(x, n, t, s, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")
  s <- as.numeric(s)

  if (length(s) == 0L || any(!is.finite(s)) || any(s < 0) || any(s > 1)) {
    stop("s must contain values in [0, 1].", call. = FALSE)
  }

  xnts <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnts$a
  n <- xnts$b
  t <- xnts$c

  if (!length(s) %in% c(1L, length(x))) {
    stop("s must have length 1 or common length.", call. = FALSE)
  }
  s <- rep(s, length.out = length(x))

  if (any(t >= n)) {
    stop("t must satisfy t < n.", call. = FALSE)
  }
  if (any(t + 1 > n)) {
    stop("t + 1 must satisfy t + 1 <= n.", call. = FALSE)
  }

  out <- numeric(length(x))

  for (j in seq_along(x)) {
    Vt <- tVxn(x = x[j], n = n[j], t = t[j], i = i, tbl = tbl, model = model, ...)
    Vt1 <- tVxn(x = x[j], n = n[j], t = t[j] + 1, i = i, tbl = tbl, model = model, ...)
    P <- Pxn(x = x[j], n = n[j], i = i, tbl = tbl, model = model, ...)

    out[j] <- (Vt + P) * (1 - s[j]) + Vt1 * s[j]
  }

  out
}

#' Fractional-duration term reserve
#'
#' Computes the interpolated reserve
#' \eqn{{}_{t+s}V = ({}_tV + P)(1-s) + {}_{t+1}V \cdot s}
#' for an n-year term insurance.
#'
#' @inheritParams tsVxn
#'
#' @return Numeric vector.
#' @examples
#' tsVxn1(40, n = 20, t = 10, s = 0.5, i = 0.05, model = "uniform", omega = 100)
#' @export
tsVxn1 <- function(x, n, t, s, i, tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  n <- .check_nonneg_integerish(n, "n")
  t <- .check_nonneg_integerish(t, "t")
  s <- as.numeric(s)

  if (length(s) == 0L || any(!is.finite(s)) || any(s < 0) || any(s > 1)) {
    stop("s must contain values in [0, 1].", call. = FALSE)
  }

  xnts <- .recycle3_ch10(x, n, t, "x", "n", "t")
  x <- xnts$a
  n <- xnts$b
  t <- xnts$c

  if (!length(s) %in% c(1L, length(x))) {
    stop("s must have length 1 or common length.", call. = FALSE)
  }
  s <- rep(s, length.out = length(x))

  if (any(t >= n)) {
    stop("t must satisfy t < n.", call. = FALSE)
  }
  if (any(t + 1 > n)) {
    stop("t + 1 must satisfy t + 1 <= n.", call. = FALSE)
  }

  out <- numeric(length(x))

  for (j in seq_along(x)) {
    Vt <- tVxn1(x = x[j], n = n[j], t = t[j], i = i, tbl = tbl, model = model, ...)
    Vt1 <- tVxn1(x = x[j], n = n[j], t = t[j] + 1, i = i, tbl = tbl, model = model, ...)
    P <- Pxn1(x = x[j], n = n[j], i = i, tbl = tbl, model = model, ...)

    out[j] <- (Vt + P) * (1 - s[j]) + Vt1 * s[j]
  }

  out
}

# -------------------------------------------------------------------------
# Gross premium and expense reserves
# -------------------------------------------------------------------------

#' Whole life gross premium reserve
#'
#' Computes the Chapter 11 prospective gross premium reserve for a
#' fully discrete whole life insurance with annual premiums and
#' renewal expenses.
#'
#' This function is intended for durations after issue, where future
#' expenses are modeled through renewal premium expenses, renewal
#' per-policy expenses, and settlement expense.
#'
#' @param x Issue age.
#' @param t Duration.
#' @param i Effective annual interest rate.
#' @param G Gross annual premium.
#' @param benefit Insurance amount. Default 1.
#' @param renewal_premium_pct Renewal percent-of-premium expense.
#' @param renewal_policy_exp Renewal per-policy expense.
#' @param settlement_exp Settlement expense at death.
#' @param tbl Optional life table object.
#' @param model Optional parametric survival model.
#' @param ... Additional model parameters.
#'
#' @return Numeric vector.
#' @examples
#' tVGx(
#'   x = 40, t = 10, i = 0.05, G = 0.03,
#'   benefit = 1, renewal_premium_pct = 0.10,
#'   renewal_policy_exp = 0.002, settlement_exp = 0.02,
#'   model = "uniform", omega = 100
#' )
#' @export
tVGx <- function(x, t, i, G,
                 benefit = 1,
                 renewal_premium_pct = 0,
                 renewal_policy_exp = 0,
                 settlement_exp = 0,
                 tbl = NULL, model = NULL, ...) {
  x <- .check_nonneg_numeric(x, "x")
  t <- .check_nonneg_integerish(t, "t")
  G <- .check_nonneg_numeric(G, "G")
  benefit <- .check_nonneg_numeric(benefit, "benefit")
  renewal_premium_pct <- as.numeric(renewal_premium_pct)
  renewal_policy_exp <- .check_nonneg_numeric(renewal_policy_exp, "renewal_policy_exp")
  settlement_exp <- .check_nonneg_numeric(settlement_exp, "settlement_exp")

  if (any(!is.finite(renewal_premium_pct)) ||
      any(renewal_premium_pct < 0) ||
      any(renewal_premium_pct > 1)) {
    stop("renewal_premium_pct must lie in [0, 1].", call. = FALSE)
  }

  vals <- .recycle_common_ch10(
    x, t, G, benefit, renewal_premium_pct, renewal_policy_exp, settlement_exp
  )
  x <- vals[[1]]
  t <- vals[[2]]
  G <- vals[[3]]
  benefit <- vals[[4]]
  renewal_premium_pct <- vals[[5]]
  renewal_policy_exp <- vals[[6]]
  settlement_exp <- vals[[7]]

  out <- numeric(length(x))

  for (j in seq_along(x)) {
    A <- Ax(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)
    a <- adotx(x = x[j] + t[j], i = i, tbl = tbl, model = model, ...)

    out[j] <- (benefit[j] + settlement_exp[j]) * A -
      ((1 - renewal_premium_pct[j]) * G[j] - renewal_policy_exp[j]) * a
  }

  out
}

#' Whole life expense reserve
#'
#' Computes the Chapter 11 expense reserve for a fully discrete whole
#' life insurance with renewal expenses.
#'
#' @inheritParams tVGx
#'
#' @return Numeric vector.
#' @examples
#' tVEx(
#'   x = 40, t = 10, i = 0.05, G = 0.03,
#'   benefit = 1, renewal_premium_pct = 0.10,
#'   renewal_policy_exp = 0.002, settlement_exp = 0.02,
#'   model = "uniform", omega = 100
#' )
#' @export
tVEx <- function(x, t, i, G,
                 benefit = 1,
                 renewal_premium_pct = 0,
                 renewal_policy_exp = 0,
                 settlement_exp = 0,
                 tbl = NULL, model = NULL, ...) {
  tVGx(
    x = x, t = t, i = i, G = G,
    benefit = benefit,
    renewal_premium_pct = renewal_premium_pct,
    renewal_policy_exp = renewal_policy_exp,
    settlement_exp = settlement_exp,
    tbl = tbl, model = model, ...
  ) - benefit * tVx(
    x = x, t = t, i = i,
    tbl = tbl, model = model, ...
  )
}

# -------------------------------------------------------------------------
# Gross gain analysis
# -------------------------------------------------------------------------

#' Total gross gain for a discrete insurance contract
#'
#' Computes the Chapter 11 total gain under gross premiums and gross reserves.
#'
#' @param VtG Gross reserve at duration t.
#' @param Vt1G Gross reserve at duration t+1.
#' @param G Gross premium.
#' @param i_actual Actual annual effective interest rate.
#' @param q_actual Actual mortality rate.
#' @param r_actual Actual percent-of-premium expense rate.
#' @param e_actual Actual per-policy expense.
#' @param s_actual Actual settlement expense.
#' @param b Benefit amount. Default 1.
#'
#' @return Numeric vector.
#' @examples
#' GTg_disc(
#'   VtG = 0.10, Vt1G = 0.12, G = 0.02,
#'   i_actual = 0.05, q_actual = 0.01,
#'   r_actual = 0.03, e_actual = 0, s_actual = 0.01, b = 1
#' )
#' @export
GTg_disc <- function(VtG, Vt1G, G,
                     i_actual, q_actual,
                     r_actual = 0, e_actual = 0, s_actual = 0,
                     b = 1) {
  vals <- .recycle_common_ch10(
    .check_nonneg_numeric(VtG, "VtG"),
    .check_nonneg_numeric(Vt1G, "Vt1G"),
    .check_nonneg_numeric(G, "G"),
    as.numeric(i_actual),
    as.numeric(q_actual),
    as.numeric(r_actual),
    .check_nonneg_numeric(e_actual, "e_actual"),
    .check_nonneg_numeric(s_actual, "s_actual"),
    .check_nonneg_numeric(b, "b")
  )

  VtG <- vals[[1]]
  Vt1G <- vals[[2]]
  G <- vals[[3]]
  i_actual <- vals[[4]]
  q_actual <- vals[[5]]
  r_actual <- vals[[6]]
  e_actual <- vals[[7]]
  s_actual <- vals[[8]]
  b <- vals[[9]]

  if (any(!is.finite(i_actual)) || any(i_actual <= -1)) {
    stop("i_actual must be finite and greater than -1.", call. = FALSE)
  }
  if (any(!is.finite(q_actual)) || any(q_actual < 0) || any(q_actual > 1)) {
    stop("q_actual must lie in [0, 1].", call. = FALSE)
  }
  if (any(!is.finite(r_actual)) || any(r_actual < 0) || any(r_actual > 1)) {
    stop("r_actual must lie in [0, 1].", call. = FALSE)
  }

  (VtG + G * (1 - r_actual) - e_actual) * (1 + i_actual) -
    ((b + s_actual) * q_actual + (1 - q_actual) * Vt1G)
}


#' Ordered gross gain decomposition
#'
#' Decomposes total gross gain into interest, mortality, and expense
#' components in a user-specified order.
#'
#' @param VtG Gross reserve at duration t.
#' @param Vt1G Gross reserve at duration t+1.
#' @param G Gross premium.
#' @param i_assumed Assumed annual effective interest rate.
#' @param q_assumed Assumed mortality rate.
#' @param r_assumed Assumed percent-of-premium expense rate.
#' @param e_assumed Assumed per-policy expense.
#' @param s_assumed Assumed settlement expense.
#' @param i_actual Actual annual effective interest rate.
#' @param q_actual Actual mortality rate.
#' @param r_actual Actual percent-of-premium expense rate.
#' @param e_actual Actual per-policy expense.
#' @param s_actual Actual settlement expense.
#' @param b Benefit amount. Default 1.
#' @param order Character vector giving the order of decomposition.
#'
#' @return Named numeric vector.
#' @examples
#' decompGg_disc(
#'   VtG = 3950.73, Vt1G = 4607.07, G = 685,
#'   i_assumed = 0.06, q_assumed = 0.00592,
#'   r_assumed = 0.05, e_assumed = 0, s_assumed = 300,
#'   i_actual = 0.065, q_actual = 0.005,
#'   r_actual = 0.06, e_actual = 0, s_actual = 100,
#'   b = 50000,
#'   order = c("interest", "mortality", "expense")
#' )
#' @export
decompGg_disc <- function(VtG, Vt1G, G,
                          i_assumed, q_assumed,
                          r_assumed = 0, e_assumed = 0, s_assumed = 0,
                          i_actual, q_actual,
                          r_actual = 0, e_actual = 0, s_actual = 0,
                          b = 1,
                          order = c("interest", "mortality", "expense")) {
  order <- match.arg(order, several.ok = TRUE)

  if (length(order) != 3L || length(unique(order)) != 3L) {
    stop("order must contain 'interest', 'mortality', and 'expense' exactly once.",
         call. = FALSE)
  }

  state_assumed <- list(
    i = i_assumed, q = q_assumed, r = r_assumed, e = e_assumed, s = s_assumed
  )
  state_actual <- list(
    i = i_actual, q = q_actual, r = r_actual, e = e_actual, s = s_actual
  )

  eval_gain <- function(st) {
    GTg_disc(
      VtG = VtG,
      Vt1G = Vt1G,
      G = G,
      i_actual = st$i,
      q_actual = st$q,
      r_actual = st$r,
      e_actual = st$e,
      s_actual = st$s,
      b = b
    )
  }

  comps <- c(interest = 0, mortality = 0, expense = 0)

  current <- state_assumed
  running_gain <- 0
  total <- eval_gain(state_actual)

  # First two gains are computed sequentially
  for (k in 1:2) {
    nm <- order[k]

    if (nm == "interest") {
      current$i <- state_actual$i
    } else if (nm == "mortality") {
      current$q <- state_actual$q
    } else if (nm == "expense") {
      current$r <- state_actual$r
      current$e <- state_actual$e
      current$s <- state_actual$s
    }

    new_gain <- eval_gain(current)
    comps[nm] <- new_gain - running_gain
    running_gain <- new_gain
  }

  # Last gain is the balancing item
  comps[order[3]] <- total - sum(comps[order[1:2]])

  c(total_gain = total, comps, check = sum(comps))
}
