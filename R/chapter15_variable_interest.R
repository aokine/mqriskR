# =========================================================
# Variable-interest models
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.variable_interest_check_numeric <- function(x, name, positive_length = TRUE) {
  if (!is.numeric(x)) {
    stop(name, " must be numeric.", call. = FALSE)
  }
  if (positive_length && length(x) == 0L) {
    stop(name, " must have positive length.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain finite values.", call. = FALSE)
  }
  as.numeric(x)
}

#' @noRd
.variable_interest_check_probability <- function(x, name) {
  x <- .variable_interest_check_numeric(x, name)
  if (any(x < 0 | x > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }
  x
}

#' @noRd
.variable_interest_check_rate <- function(x, name) {
  x <- .variable_interest_check_numeric(x, name)
  if (any(x <= -1)) {
    stop(name, " must contain values greater than -1.", call. = FALSE)
  }
  x
}

#' @noRd
.variable_interest_check_nonnegative <- function(x, name) {
  x <- .variable_interest_check_numeric(x, name)
  if (any(x < 0)) {
    stop(name, " must contain nonnegative values.", call. = FALSE)
  }
  x
}

#' @noRd
.variable_interest_check_scalar <- function(
    x,
    name,
    lower = -Inf,
    strict_lower = FALSE) {
  x <- .variable_interest_check_numeric(x, name)

  if (length(x) != 1L) {
    stop(name, " must be a numeric scalar.", call. = FALSE)
  }

  invalid <- if (strict_lower) x <= lower else x < lower
  if (invalid) {
    comparator <- if (strict_lower) "greater than" else "at least"
    stop(name, " must be ", comparator, " ", lower, ".", call. = FALSE)
  }

  x
}

#' @noRd
.variable_interest_check_integer <- function(
    x,
    name,
    positive = FALSE) {
  x <- .variable_interest_check_scalar(
    x,
    name,
    lower = 0,
    strict_lower = positive
  )

  if (abs(x - round(x)) > 1e-10) {
    descriptor <- if (positive) "positive" else "nonnegative"
    stop(name, " must be a ", descriptor, " integer.", call. = FALSE)
  }

  as.integer(round(x))
}

#' @noRd
.variable_interest_check_same_length <- function(x, y, x_name, y_name) {
  if (length(x) != length(y)) {
    stop(x_name, " and ", y_name, " must have the same length.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.variable_interest_survival_start <- function(qx) {
  c(1, cumprod(1 - qx))[seq_along(qx)]
}

#' @noRd
.variable_interest_survival_end <- function(qx) {
  cumprod(1 - qx)
}

#' @noRd
.variable_interest_benefit_scalar <- function(benefit) {
  benefit <- .variable_interest_check_nonnegative(benefit, "benefit")
  if (length(benefit) != 1L) {
    stop("benefit must be a numeric scalar.", call. = FALSE)
  }
  benefit
}

# -------------------------------------------------------------------------
# Variable-interest discount factors
# -------------------------------------------------------------------------

#' Discount factors under variable annual interest rates
#'
#' Computes cumulative discount factors for a sequence of annual effective
#' interest rates:
#' \deqn{
#' v_t = \prod_{k=1}^{t}(1+i_k)^{-1},
#' \qquad t=1,\ldots,n.
#' }
#'
#' @param i Numeric vector of annual effective interest rates. Each value must
#'   be greater than \code{-1}.
#'
#' @return A numeric vector of cumulative discount factors with the same
#'   length as \code{i}.
#'
#' @examples
#' vt_var(c(0.06, 0.07, 0.08))
#' vt_var(c(-0.01, 0.02, 0.03))
#'
#' @export
vt_var <- function(i) {
  i <- .variable_interest_check_rate(i, "i")
  cumprod(1 / (1 + i))
}


# -------------------------------------------------------------------------
# APVs under a sequence of annual effective interest rates
# -------------------------------------------------------------------------

#' Actuarial present values under variable annual interest rates
#'
#' Computes life-contingent actuarial present values using a specified
#' sequence of annual effective interest rates.
#'
#' The vectors \code{qx} and \code{i} represent one valuation scenario and
#' must have the same positive length.
#'
#' \code{nEx_var()} computes a pure endowment.
#'
#' \code{Axn1_var()} computes term insurance payable at the end of the year
#' of death.
#'
#' \code{Axn_var()} computes endowment insurance.
#'
#' \code{axn_var()} computes a temporary annuity-immediate or
#' annuity-due.
#'
#' @param qx Numeric vector of one-year mortality probabilities.
#' @param i Numeric vector of annual effective interest rates. Each value
#'   must be greater than \code{-1}.
#' @param benefit Nonnegative scalar benefit or annuity payment amount.
#' @param type Character string equal to \code{"immediate"} or
#'   \code{"due"}.
#'
#' @return A numeric scalar.
#'
#' @details
#' Each year's payment is discounted using the cumulative product of the
#' annual effective interest rates supplied in \code{i}.
#'
#' The pure endowment is
#' \deqn{
#' {}_nE
#' =
#' {}_np_x\,v_n,
#' }
#' where \eqn{v_n} is the cumulative discount factor implied by the sequence
#' of annual effective interest rates.
#'
#' Term insurance is obtained by discounting each possible death benefit
#' using the cumulative discount factor applicable to its payment year.
#'
#' Endowment insurance equals the sum of the corresponding term insurance
#' and pure endowment.
#'
#' Temporary annuities discount each payment using the cumulative discount
#' factors derived from the interest-rate sequence.
#'
#' @examples
#' qx <- c(0.03, 0.04, 0.05, 0.06, 0.07)
#' rates <- c(0.06, 0.07, 0.08, 0.09, 0.10)
#'
#' nEx_var(qx, rates, benefit = 1000)
#' Axn1_var(qx, rates)
#' Axn_var(qx, rates)
#' axn_var(qx, rates, type = "due")
#'
#' @name variable_interest_apvs
#' @rdname variable_interest_apvs
#' @export
nEx_var <- function(qx, i, benefit = 1) {
  qx <- .variable_interest_check_probability(qx, "qx")
  i <- .variable_interest_check_rate(i, "i")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, i, "qx", "i")

  benefit * prod(1 - qx) * tail(vt_var(i), 1L)
}

#' @rdname variable_interest_apvs
#' @export
Axn1_var <- function(qx, i, benefit = 1) {
  qx <- .variable_interest_check_probability(qx, "qx")
  i <- .variable_interest_check_rate(i, "i")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, i, "qx", "i")

  discount <- vt_var(i)
  survival_start <- .variable_interest_survival_start(qx)

  benefit * sum(discount * survival_start * qx)
}

#' @rdname variable_interest_apvs
#' @export
Axn_var <- function(qx, i, benefit = 1) {
  qx <- .variable_interest_check_probability(qx, "qx")
  i <- .variable_interest_check_rate(i, "i")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, i, "qx", "i")

  Axn1_var(qx = qx, i = i, benefit = benefit) +
    nEx_var(qx = qx, i = i, benefit = benefit)
}

#' @rdname variable_interest_apvs
#' @export
axn_var <- function(qx, i, type = c("immediate", "due"), benefit = 1) {
  type <- match.arg(type)

  qx <- .variable_interest_check_probability(qx, "qx")
  i <- .variable_interest_check_rate(i, "i")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, i, "qx", "i")

  discount <- vt_var(i)

  if (type == "immediate") {
    survival_end <- .variable_interest_survival_end(qx)
    return(benefit * sum(discount * survival_end))
  }

  survival_start <- .variable_interest_survival_start(qx)
  discount_due <- c(1, head(discount, -1L))

  benefit * sum(discount_due * survival_start)
}

# -------------------------------------------------------------------------
# APVs under spot rates
# -------------------------------------------------------------------------

#' Actuarial present values under spot rates
#'
#' Computes life-contingent actuarial present values using annual effective
#' spot rates by maturity.
#'
#' If \eqn{z_t} denotes the annual effective spot rate for maturity
#' \eqn{t}, the corresponding discount factor is
#' \deqn{(1+z_t)^{-t}.}
#'
#' \code{nEx_spot()} computes a pure endowment.
#'
#' \code{Axn1_spot()} computes term insurance payable at the end of the year
#' of death.
#'
#' \code{Axn_spot()} computes endowment insurance.
#'
#' \code{axn_spot()} computes a temporary annuity-immediate or
#' annuity-due.
#'
#' @param qx Numeric vector of one-year mortality probabilities.
#' @param z Numeric vector of annual effective spot rates for maturities
#'   \code{1, \dots, n}. Each value must be greater than \code{-1}.
#' @param benefit Nonnegative scalar benefit or annuity payment amount.
#' @param type Character string equal to \code{"immediate"} or
#'   \code{"due"}.
#'
#' @return A numeric scalar.
#'
#' @details
#' Each payment is discounted using the spot rate corresponding to its
#' maturity rather than a single level interest rate.
#'
#' The pure endowment is
#' \deqn{
#' {}_nE
#' =
#' {}_np_x(1+z_n)^{-n}.
#' }
#'
#' Term insurance is obtained by discounting each death benefit using the
#' spot rate corresponding to its payment year.
#'
#' Endowment insurance equals the sum of the corresponding term insurance
#' and pure endowment.
#'
#' Temporary annuities discount each payment using the spot rate for its
#' payment time.
#'
#' @examples
#' qx <- c(0.02, 0.03, 0.04, 0.05, 0.06)
#' spot <- c(0.03, 0.04, 0.05, 0.06, 0.07)
#'
#' nEx_spot(qx, spot, benefit = 1000)
#' Axn1_spot(qx, spot)
#' Axn_spot(qx, spot)
#' axn_spot(qx, spot, type = "due")
#'
#' @name spot_interest_apvs
#' @rdname spot_interest_apvs
#' @export
nEx_spot <- function(qx, z, benefit = 1) {
  qx <- .variable_interest_check_probability(qx, "qx")
  z <- .variable_interest_check_rate(z, "z")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, z, "qx", "z")

  n <- length(z)
  benefit * prod(1 - qx) * (1 + z[n])^(-n)
}

#' @rdname spot_interest_apvs
#' @export
Axn1_spot <- function(qx, z, benefit = 1) {
  qx <- .variable_interest_check_probability(qx, "qx")
  z <- .variable_interest_check_rate(z, "z")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, z, "qx", "z")

  times <- seq_along(z)
  discount <- (1 + z)^(-times)
  survival_start <- .variable_interest_survival_start(qx)

  benefit * sum(discount * survival_start * qx)
}

#' @rdname spot_interest_apvs
#' @export
Axn_spot <- function(qx, z, benefit = 1) {
  qx <- .variable_interest_check_probability(qx, "qx")
  z <- .variable_interest_check_rate(z, "z")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, z, "qx", "z")

  Axn1_spot(qx = qx, z = z, benefit = benefit) +
    nEx_spot(qx = qx, z = z, benefit = benefit)
}

#' @rdname spot_interest_apvs
#' @export
axn_spot <- function(qx, z, type = c("immediate", "due"), benefit = 1) {
  type <- match.arg(type)

  qx <- .variable_interest_check_probability(qx, "qx")
  z <- .variable_interest_check_rate(z, "z")
  benefit <- .variable_interest_benefit_scalar(benefit)

  .variable_interest_check_same_length(qx, z, "qx", "z")

  n <- length(z)

  if (type == "immediate") {
    discount <- (1 + z)^(-seq_len(n))
    survival_end <- .variable_interest_survival_end(qx)
    return(benefit * sum(discount * survival_end))
  }

  survival_start <- .variable_interest_survival_start(qx)
  discount_due <- c(
    1,
    (1 + head(z, -1L))^(-seq_len(n - 1L))
  )

  benefit * sum(discount_due * survival_start)
}

# -------------------------------------------------------------------------
# Spot-rate valuation and bootstrapping
# -------------------------------------------------------------------------

#' Present value of deterministic cash flows using spot rates
#'
#' Discounts each cash flow using the spot rate matched to its payment time.
#' Time-zero cash flows are not discounted.
#'
#' @param amounts Numeric vector of cash-flow amounts.
#' @param times Numeric vector of payment times in years.
#' @param spot Numeric vector of spot rates matched elementwise to
#'   \code{times}. The value corresponding to a time-zero cash flow is ignored.
#' @param compounding Character string equal to \code{"annual"} or
#'   \code{"semiannual"}.
#'
#' @return A numeric scalar.
#'
#' @examples
#' pv_spot_cashflows(
#'   amounts = c(200000, 50000, 50000, 100000),
#'   times = c(0, 0.5, 1, 2),
#'   spot = c(0, 0.02440, 0.02601, 0.02936),
#'   compounding = "semiannual"
#' )
#'
#' @export
pv_spot_cashflows <- function(
    amounts,
    times,
    spot,
    compounding = c("annual", "semiannual")) {
  compounding <- match.arg(compounding)

  amounts <- .variable_interest_check_numeric(amounts, "amounts")
  times <- .variable_interest_check_nonnegative(times, "times")
  spot <- .variable_interest_check_numeric(spot, "spot")

  if (!(length(amounts) == length(times) &&
        length(times) == length(spot))) {
    stop(
      "amounts, times, and spot must have the same length.",
      call. = FALSE
    )
  }

  positive_time <- times > 0
  if (compounding == "annual" &&
      any(1 + spot[positive_time] <= 0)) {
    stop(
      "Annual spot rates at positive times must be greater than -1.",
      call. = FALSE
    )
  }
  if (compounding == "semiannual" &&
      any(1 + spot[positive_time] / 2 <= 0)) {
    stop(
      "Semiannual nominal spot rates at positive times must be greater than -2.",
      call. = FALSE
    )
  }

  discount <- rep(1, length(times))

  if (compounding == "annual") {
    discount[positive_time] <-
      (1 + spot[positive_time])^(-times[positive_time])
  } else {
    discount[positive_time] <-
      (1 + spot[positive_time] / 2)^(-2 * times[positive_time])
  }

  sum(amounts * discount)
}

#' Bootstrap semiannual nominal spot rates
#'
#' Bootstraps nominal annual spot rates convertible semiannually from par
#' coupon yields at consecutive half-year maturities.
#'
#' @param maturity Numeric vector of positive maturities in years, in strictly
#'   increasing order. Maturities must be consecutive multiples of
#'   \code{0.5}.
#' @param coupon_yield Numeric vector of nominal annual par coupon yields
#'   convertible semiannually. Values must be greater than \code{-2}.
#' @param par Positive scalar par value.
#'
#' @return A numeric vector of nominal annual spot rates convertible
#'   semiannually.
#'
#' @examples
#' maturity <- c(0.5, 1.0, 1.5, 2.0)
#' coupon_yield <- c(0.0244, 0.0260, 0.0276, 0.0293)
#' z_from_coupon_semi(maturity, coupon_yield)
#'
#' @export
z_from_coupon_semi <- function(maturity, coupon_yield, par = 1000) {
  maturity <- .variable_interest_check_numeric(maturity, "maturity")
  coupon_yield <- .variable_interest_check_numeric(
    coupon_yield,
    "coupon_yield"
  )
  par <- .variable_interest_check_scalar(
    par,
    "par",
    lower = 0,
    strict_lower = TRUE
  )

  .variable_interest_check_same_length(
    maturity,
    coupon_yield,
    "maturity",
    "coupon_yield"
  )

  if (any(maturity <= 0) || is.unsorted(maturity, strictly = TRUE)) {
    stop(
      "maturity must be positive and strictly increasing.",
      call. = FALSE
    )
  }
  if (any(coupon_yield <= -2)) {
    stop(
      "coupon_yield must contain values greater than -2.",
      call. = FALSE
    )
  }

  periods <- round(2 * maturity)
  if (any(abs(2 * maturity - periods) > 1e-10)) {
    stop(
      "maturity must contain multiples of 0.5 years.",
      call. = FALSE
    )
  }
  if (!identical(as.integer(periods), seq_len(length(periods)))) {
    stop(
      "maturity must contain consecutive half-year maturities beginning at 0.5.",
      call. = FALSE
    )
  }

  spot <- numeric(length(maturity))

  for (j in seq_along(maturity)) {
    coupon <- par * coupon_yield[j] / 2
    m <- as.integer(periods[j])

    if (m == 1L) {
      spot[j] <- coupon_yield[j]
      next
    }

    prior_periods <- seq_len(m - 1L)
    pv_prior <- sum(
      coupon / (1 + spot[prior_periods] / 2)^prior_periods
    )
    final_cash_flow <- par + coupon
    remaining_price <- par - pv_prior

    if (remaining_price <= 0) {
      stop(
        "The supplied coupon yields do not produce a positive final discount factor.",
        call. = FALSE
      )
    }

    accumulation <- final_cash_flow / remaining_price
    spot[j] <- 2 * (accumulation^(1 / m) - 1)
  }

  spot
}

#' Bootstrap annual effective spot rates
#'
#' Bootstraps annual effective spot rates from par coupon yields at consecutive
#' integer maturities.
#'
#' @param maturity Numeric vector of positive integer maturities in strictly
#'   increasing order. Maturities must be consecutive and begin at 1.
#' @param coupon_yield Numeric vector of annual effective par coupon yields.
#'   Values must be greater than \code{-1}.
#' @param par Positive scalar par value.
#'
#' @return A numeric vector of annual effective spot rates.
#'
#' @examples
#' maturity <- 1:4
#' coupon_yield <- c(0.02, 0.04, 0.06, 0.08)
#' z_from_coupon_annual(maturity, coupon_yield)
#'
#' @export
z_from_coupon_annual <- function(maturity, coupon_yield, par = 1000) {
  maturity <- .variable_interest_check_numeric(maturity, "maturity")
  coupon_yield <- .variable_interest_check_rate(
    coupon_yield,
    "coupon_yield"
  )
  par <- .variable_interest_check_scalar(
    par,
    "par",
    lower = 0,
    strict_lower = TRUE
  )

  .variable_interest_check_same_length(
    maturity,
    coupon_yield,
    "maturity",
    "coupon_yield"
  )

  if (any(maturity <= 0) || is.unsorted(maturity, strictly = TRUE)) {
    stop(
      "maturity must be positive and strictly increasing.",
      call. = FALSE
    )
  }
  if (any(abs(maturity - round(maturity)) > 1e-10)) {
    stop("maturity must contain integer years.", call. = FALSE)
  }

  maturity <- as.integer(round(maturity))
  if (!identical(maturity, seq_len(length(maturity)))) {
    stop(
      "maturity must contain consecutive integer maturities beginning at 1.",
      call. = FALSE
    )
  }

  spot <- numeric(length(maturity))

  for (j in seq_along(maturity)) {
    coupon <- par * coupon_yield[j]
    m <- maturity[j]

    if (m == 1L) {
      spot[j] <- coupon_yield[j]
      next
    }

    prior_periods <- seq_len(m - 1L)
    pv_prior <- sum(coupon / (1 + spot[prior_periods])^prior_periods)
    final_cash_flow <- par + coupon
    remaining_price <- par - pv_prior

    if (remaining_price <= 0) {
      stop(
        "The supplied coupon yields do not produce a positive final discount factor.",
        call. = FALSE
      )
    }

    accumulation <- final_cash_flow / remaining_price
    spot[j] <- accumulation^(1 / m) - 1
  }

  spot
}

# -------------------------------------------------------------------------
# Forward and spot rate conversions
# -------------------------------------------------------------------------

#' Forward rate implied by spot rates
#'
#' Computes the annual effective forward rate for the interval from time
#' \code{n} to time \code{n + k}:
#' \deqn{
#' (1+z_{n+k})^{n+k}
#' =
#' (1+z_n)^n(1+f_{n,k})^k.
#' }
#'
#' @param z Numeric vector of annual effective spot rates for maturities
#'   \code{1, ..., length(z)}. Each value must be greater than \code{-1}.
#' @param n Nonnegative integer forward-start time.
#' @param k Positive integer forward period.
#'
#' @return A numeric scalar.
#'
#' @examples
#' z <- c(0.03, 0.04, 0.05, 0.06, 0.07)
#' fnk_from_z(z, n = 1, k = 4)
#' fnk_from_z(z, n = 2, k = 2)
#'
#' @export
fnk_from_z <- function(z, n, k) {
  z <- .variable_interest_check_rate(z, "z")
  n <- .variable_interest_check_integer(n, "n")
  k <- .variable_interest_check_integer(k, "k", positive = TRUE)

  if (n + k > length(z)) {
    stop("Spot rates are required through maturity n + k.", call. = FALSE)
  }

  z_n <- if (n == 0L) 0 else z[n]
  z_nk <- z[n + k]

  ((1 + z_nk)^(n + k) / (1 + z_n)^n)^(1 / k) - 1
}

#' Matrix of forward rates implied by spot rates
#'
#' Constructs a matrix of annual effective forward rates. Rows correspond to
#' forward-start times \eqn{n=1,\ldots,m-1}; columns correspond to forward
#' periods \eqn{k=1,\ldots,m-1}. Entries requiring maturities beyond
#' \eqn{m} are returned as \code{NA}.
#'
#' @param z Numeric vector of at least two annual effective spot rates. Each
#'   value must be greater than \code{-1}.
#'
#' @return A numeric matrix.
#'
#' @examples
#' forward_matrix_from_z(c(0.03, 0.04, 0.05, 0.06, 0.07))
#'
#' @export
forward_matrix_from_z <- function(z) {
  z <- .variable_interest_check_rate(z, "z")

  if (length(z) < 2L) {
    stop("z must contain at least two spot rates.", call. = FALSE)
  }

  m <- length(z)
  out <- matrix(
    NA_real_,
    nrow = m - 1L,
    ncol = m - 1L,
    dimnames = list(
      paste0("n=", seq_len(m - 1L)),
      paste0("k=", seq_len(m - 1L))
    )
  )

  for (n in seq_len(m - 1L)) {
    for (k in seq_len(m - n)) {
      out[n, k] <- fnk_from_z(z, n = n, k = k)
    }
  }

  out
}

#' Spot rates from one-year forward rates
#'
#' Converts annual effective one-year forward rates
#' \eqn{f_{0,1},f_{1,1},\ldots,f_{n-1,1}} into annual effective spot rates:
#' \deqn{
#' (1+z_n)^n
#' =
#' \prod_{j=0}^{n-1}(1+f_{j,1}).
#' }
#'
#' @param fn1 Numeric vector of annual effective one-year forward rates. Each
#'   value must be greater than \code{-1}.
#'
#' @return A numeric vector of annual effective spot rates.
#'
#' @examples
#' z_from_fn1(c(0.04, 0.05, 0.06, 0.07, 0.08))
#'
#' @export
z_from_fn1 <- function(fn1) {
  fn1 <- .variable_interest_check_rate(fn1, "fn1")

  accumulation <- cumprod(1 + fn1)
  accumulation^(1 / seq_along(fn1)) - 1
}
