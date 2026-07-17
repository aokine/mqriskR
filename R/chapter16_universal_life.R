# =========================================================
# Universal life insurance functions
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' @noRd
.ul_check_numeric <- function(x, name, allow_infinite = FALSE) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop(name, " must be a numeric vector with positive length.", call. = FALSE)
  }

  invalid <- if (allow_infinite) {
    is.na(x) | is.nan(x)
  } else {
    !is.finite(x)
  }

  if (any(invalid)) {
    descriptor <- if (allow_infinite) {
      "must not contain missing or NaN values."
    } else {
      "must contain finite values."
    }
    stop(name, " ", descriptor, call. = FALSE)
  }

  as.numeric(x)
}

#' @noRd
.ul_check_nonnegative <- function(x, name, allow_infinite = FALSE) {
  x <- .ul_check_numeric(x, name, allow_infinite = allow_infinite)

  if (any(x < 0)) {
    stop(name, " must contain nonnegative values.", call. = FALSE)
  }

  x
}

#' @noRd
.ul_check_probability <- function(x, name) {
  x <- .ul_check_numeric(x, name)

  if (any(x < 0 | x > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }

  x
}

#' @noRd
.ul_check_rate <- function(x, name) {
  x <- .ul_check_numeric(x, name)

  if (any(x <= -1)) {
    stop(name, " must contain values greater than -1.", call. = FALSE)
  }

  x
}

#' @noRd
.ul_check_scalar <- function(x, name, nonnegative = FALSE, positive = FALSE) {
  if (nonnegative || positive) {
    x <- .ul_check_nonnegative(x, name)
  } else {
    x <- .ul_check_numeric(x, name)
  }

  if (length(x) != 1L) {
    stop(name, " must be a numeric scalar.", call. = FALSE)
  }

  if (positive && x <= 0) {
    stop(name, " must be positive.", call. = FALSE)
  }

  x
}

#' @noRd
.ul_check_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }

  x
}

#' @noRd
.ul_recycle <- function(..., .names = NULL) {
  values <- list(...)
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

# -------------------------------------------------------------------------
# Cost of insurance
# -------------------------------------------------------------------------

#' Cost of insurance for Type B universal life
#'
#' Computes the one-period cost of insurance for a Type B universal life
#' contract:
#' \deqn{
#' \mathrm{COI}_t = \frac{B_t q_{x+t-1}}{1+i_t^q}.
#' }
#'
#' The arguments may be scalars or vectors. Scalar arguments are recycled to
#' the common length.
#'
#' @param B Face amount.
#' @param qx Mortality probability for the period.
#' @param iq Interest rate used in the cost-of-insurance calculation. Values
#'   must be greater than \code{-1}.
#'
#' @return A numeric vector of cost-of-insurance charges.
#'
#' @examples
#' coi_ul_typeB(B = 100000, qx = 0.00076, iq = 0.03)
#' coi_ul_typeB(
#'   B = 100000,
#'   qx = c(0.00076, 0.00081),
#'   iq = c(0.03, 0.035)
#' )
#'
#' @export
coi_ul_typeB <- function(B, qx, iq) {
  B <- .ul_check_nonnegative(B, "B")
  qx <- .ul_check_probability(qx, "qx")
  iq <- .ul_check_rate(iq, "iq")

  values <- .ul_recycle(B, qx, iq, .names = c("B", "qx", "iq"))

  B <- values[[1]]
  qx <- values[[2]]
  iq <- values[[3]]

  B * qx / (1 + iq)
}

# -------------------------------------------------------------------------
# Universal life account-value paths
# -------------------------------------------------------------------------

#' Type B universal life account-value path
#'
#' Computes a year-by-year Type B universal life account-value path. Premiums
#' and expenses are applied at the beginning of each period, followed by the
#' cost-of-insurance charge and credited interest.
#'
#' @param G Premium amount by period.
#' @param r Percent-of-premium expense rate by period. Values must lie in
#'   \code{[0, 1]}.
#' @param e Fixed expense by period.
#' @param qx Mortality probability by period.
#' @param ic Credited annual effective interest rate by period. Values must be
#'   greater than \code{-1}.
#' @param B Face amount by period.
#' @param iq Interest rate used in the cost-of-insurance calculation. Defaults
#'   to \code{ic}; values must be greater than \code{-1}.
#' @param AV0 Nonnegative scalar initial account value.
#'
#' @return A data frame containing the policy duration, premium,
#'   net contribution, cost-of-insurance charge, and account value.
#'
#' @examples
#' qx <- c(0.00076, 0.00081, 0.00085, 0.00090, 0.00095)
#' r <- c(0.75, rep(0.10, 4))
#' e <- c(100, rep(20, 4))
#'
#' AV_path_ul_typeB(
#'   G = 5000,
#'   r = r,
#'   e = e,
#'   qx = qx,
#'   ic = 0.03,
#'   B = 100000
#' )
#'
#' @export
AV_path_ul_typeB <- function(G, r, e, qx, ic, B, iq = ic, AV0 = 0) {
  G <- .ul_check_nonnegative(G, "G")
  r <- .ul_check_probability(r, "r")
  e <- .ul_check_nonnegative(e, "e")
  qx <- .ul_check_probability(qx, "qx")
  ic <- .ul_check_rate(ic, "ic")
  B <- .ul_check_nonnegative(B, "B")
  iq <- .ul_check_rate(iq, "iq")
  AV0 <- .ul_check_scalar(AV0, "AV0", nonnegative = TRUE)

  values <- .ul_recycle(
    G, r, e, qx, ic, B, iq,
    .names = c("G", "r", "e", "qx", "ic", "B", "iq")
  )

  G <- values[[1]]
  r <- values[[2]]
  e <- values[[3]]
  qx <- values[[4]]
  ic <- values[[5]]
  B <- values[[6]]
  iq <- values[[7]]

  n <- length(qx)
  account_value <- numeric(n + 1L)
  net_contribution <- numeric(n)
  cost_of_insurance <- numeric(n)
  account_value[1L] <- AV0

  for (k in seq_len(n)) {
    net_contribution[k] <- G[k] * (1 - r[k]) - e[k]
    cost_of_insurance[k] <- coi_ul_typeB(B[k], qx[k], iq[k])
    account_value[k + 1L] <-
      (account_value[k] + net_contribution[k] - cost_of_insurance[k]) *
      (1 + ic[k])
  }

  data.frame(
    t = 0:n,
    premium = c(NA_real_, G),
    net_contribution = c(NA_real_, net_contribution),
    COI = c(NA_real_, cost_of_insurance),
    AV = account_value
  )
}

#' Type A universal life account-value path
#'
#' Computes a year-by-year Type A universal life account-value path. The death
#' benefit is fixed, so the net amount at risk depends on the ending account
#' value and the roll-forward is solved explicitly each period.
#'
#' @inheritParams AV_path_ul_typeB
#'
#' @return A data frame containing the policy duration, premium, and account
#'   value.
#'
#' @examples
#' qx <- c(0.00076, 0.00081, 0.00085, 0.00090, 0.00095)
#' r <- c(0.75, rep(0.10, 4))
#' e <- c(100, rep(20, 4))
#'
#' AV_path_ul_typeA(
#'   G = 5000,
#'   r = r,
#'   e = e,
#'   qx = qx,
#'   ic = 0.03,
#'   B = 100000
#' )
#'
#' @export
AV_path_ul_typeA <- function(G, r, e, qx, ic, B, iq = ic, AV0 = 0) {
  G <- .ul_check_nonnegative(G, "G")
  r <- .ul_check_probability(r, "r")
  e <- .ul_check_nonnegative(e, "e")
  qx <- .ul_check_probability(qx, "qx")
  ic <- .ul_check_rate(ic, "ic")
  B <- .ul_check_nonnegative(B, "B")
  iq <- .ul_check_rate(iq, "iq")
  AV0 <- .ul_check_scalar(AV0, "AV0", nonnegative = TRUE)

  values <- .ul_recycle(
    G, r, e, qx, ic, B, iq,
    .names = c("G", "r", "e", "qx", "ic", "B", "iq")
  )

  G <- values[[1]]
  r <- values[[2]]
  e <- values[[3]]
  qx <- values[[4]]
  ic <- values[[5]]
  B <- values[[6]]
  iq <- values[[7]]

  n <- length(qx)
  account_value <- numeric(n + 1L)
  account_value[1L] <- AV0

  for (k in seq_len(n)) {
    mortality_factor <- qx[k] * (1 + ic[k]) / (1 + iq[k])
    denominator <- 1 - mortality_factor

    if (denominator <= 0) {
      stop("The Type A roll-forward denominator must be positive.", call. = FALSE)
    }

    numerator <-
      (account_value[k] + G[k] * (1 - r[k]) - e[k]) * (1 + ic[k]) -
      B[k] * mortality_factor

    account_value[k + 1L] <- numerator / denominator
  }

  data.frame(
    t = 0:n,
    premium = c(NA_real_, G),
    AV = account_value
  )
}

# -------------------------------------------------------------------------
# Indexed universal life credited-rate helpers
# -------------------------------------------------------------------------

#' Point-to-point index growth rates
#'
#' Computes consecutive point-to-point growth rates from index values.
#'
#' @param index Numeric vector of strictly positive index values.
#'
#' @return A numeric vector with length one less than \code{index}.
#'
#' @examples
#' iP_eiul(c(1000, 1050, 1200, 1100))
#'
#' @export
iP_eiul <- function(index) {
  index <- .ul_check_numeric(index, "index")

  if (length(index) < 2L) {
    stop("index must contain at least two values.", call. = FALSE)
  }
  if (any(index <= 0)) {
    stop("index values must be positive.", call. = FALSE)
  }

  index[-1L] / index[-length(index)] - 1
}

#' Monthly-average index growth rate
#'
#' Computes the monthly-average growth rate from an initial index value and
#' twelve monthly closing values.
#'
#' @param index Numeric vector of length 13 containing a strictly positive
#'   initial index value followed by twelve nonnegative monthly closing values.
#'
#' @return A numeric scalar.
#'
#' @examples
#' index <- c(
#'   1000, 1020, 1100, 1150, 1080, 1040, 960,
#'   1030, 1000, 1070, 1150, 1200, 1150
#' )
#' iMA_eiul(index)
#'
#' @export
iMA_eiul <- function(index) {
  index <- .ul_check_nonnegative(index, "index")

  if (length(index) != 13L) {
    stop(
      "index must contain an initial value and 12 monthly closing values.",
      call. = FALSE
    )
  }
  if (index[1L] <= 0) {
    stop("The initial index value must be positive.", call. = FALSE)
  }

  mean(index[-1L]) / index[1L] - 1
}

#' Credited rates from index growth rates
#'
#' Applies a participation rate, floor, cap, and optional margin to raw index
#' growth rates.
#'
#' @param i_raw Numeric vector of raw index growth rates.
#' @param part Nonnegative scalar participation rate.
#' @param floor Scalar minimum credited rate.
#' @param cap Scalar maximum credited rate. The default is \code{Inf}.
#' @param margin Nonnegative scalar index margin.
#' @param margin_after_participation Logical scalar. If \code{TRUE}, the margin
#'   is subtracted after applying participation; otherwise it is subtracted
#'   before applying participation.
#'
#' @return A numeric vector of credited rates.
#'
#' @examples
#' raw <- iP_eiul(c(1000, 1050, 1200, 1100))
#' i_credit_eiul(raw, part = 1.10, floor = 0.01, cap = 0.10)
#' i_credit_eiul(raw)
#'
#' @export
i_credit_eiul <- function(i_raw, part = 1, floor = 0, cap = Inf,
                          margin = 0, margin_after_participation = TRUE) {
  i_raw <- .ul_check_numeric(i_raw, "i_raw")
  part <- .ul_check_scalar(part, "part", nonnegative = TRUE)
  floor <- .ul_check_scalar(floor, "floor")
  cap <- .ul_check_numeric(cap, "cap", allow_infinite = TRUE)
  margin <- .ul_check_scalar(margin, "margin", nonnegative = TRUE)
  margin_after_participation <- .ul_check_logical_scalar(
    margin_after_participation,
    "margin_after_participation"
  )

  if (length(cap) != 1L) {
    stop("cap must be a numeric scalar.", call. = FALSE)
  }
  if (floor > cap) {
    stop("floor cannot exceed cap.", call. = FALSE)
  }

  adjusted <- if (margin_after_participation) {
    part * i_raw - margin
  } else {
    part * (i_raw - margin)
  }

  pmin(cap, pmax(floor, adjusted))
}

# -------------------------------------------------------------------------
# Persistency under mortality and withdrawal
# -------------------------------------------------------------------------

#' Universal life persistency probabilities
#'
#' \code{pxtau_ul()} computes one-year persistency probabilities under
#' mortality and withdrawal.
#'
#' \code{tpxtau_ul()} computes cumulative persistency through the end of each
#' policy year.
#'
#' @param qd Mortality probabilities.
#' @param qw Withdrawal probabilities.
#' @param year_end_withdrawal Logical scalar. If \code{TRUE}, withdrawal is
#'   modeled at year-end and persistency is
#'   \eqn{(1-q^{(d)})(1-q^{(w)})}. Otherwise persistency is
#'   \eqn{1-q^{(d)}-q^{(w)}}.
#'
#' @return A numeric vector.
#'
#' @examples
#' qd <- c(0.001, 0.002, 0.003)
#' qw <- c(0.02, 0.02, 0.03)
#'
#' pxtau_ul(qd, qw)
#' tpxtau_ul(qd, qw)
#'
#' @rdname ul_persistency
#' @export
pxtau_ul <- function(qd, qw, year_end_withdrawal = TRUE) {
  qd <- .ul_check_probability(qd, "qd")
  qw <- .ul_check_probability(qw, "qw")
  year_end_withdrawal <- .ul_check_logical_scalar(
    year_end_withdrawal,
    "year_end_withdrawal"
  )

  values <- .ul_recycle(qd, qw, .names = c("qd", "qw"))
  qd <- values[[1]]
  qw <- values[[2]]

  if (year_end_withdrawal) {
    return((1 - qd) * (1 - qw))
  }

  persistency <- 1 - qd - qw

  if (any(persistency < 0)) {
    stop("qd + qw must not exceed 1.", call. = FALSE)
  }

  persistency
}

#' @rdname ul_persistency
#' @export
tpxtau_ul <- function(qd, qw, year_end_withdrawal = TRUE) {
  cumprod(
    pxtau_ul(
      qd = qd,
      qw = qw,
      year_end_withdrawal = year_end_withdrawal
    )
  )
}

# -------------------------------------------------------------------------
# Universal life reserving helpers
# -------------------------------------------------------------------------

#' Guaranteed maturity fund roll-forward
#'
#' Computes a one-period guaranteed maturity fund roll-forward. Arguments may
#' be scalars or vectors and follow common-length recycling.
#'
#' @param GMF_prev Prior guaranteed maturity fund.
#' @param GMP Guaranteed maturity premium.
#' @param r Percent-of-premium expense rate in \code{[0, 1]}.
#' @param policy_charge Guaranteed policy charge.
#' @param i Guaranteed annual effective interest rate greater than \code{-1}.
#'
#' @return A numeric vector.
#'
#' @examples
#' GMF_rollforward_ul(140.40, 14.49, 0.04, 11.80, 0.03)
#'
#' @export
GMF_rollforward_ul <- function(GMF_prev, GMP, r, policy_charge, i) {
  GMF_prev <- .ul_check_nonnegative(GMF_prev, "GMF_prev")
  GMP <- .ul_check_nonnegative(GMP, "GMP")
  r <- .ul_check_probability(r, "r")
  policy_charge <- .ul_check_nonnegative(policy_charge, "policy_charge")
  i <- .ul_check_rate(i, "i")

  values <- .ul_recycle(
    GMF_prev, GMP, r, policy_charge, i,
    .names = c("GMF_prev", "GMP", "r", "policy_charge", "i")
  )

  GMF_prev <- values[[1]]
  GMP <- values[[2]]
  r <- values[[3]]
  policy_charge <- values[[4]]
  i <- values[[5]]

  (GMF_prev + GMP * (1 - r) - policy_charge) * (1 + i)
}

#' Account-value to guaranteed-fund ratio
#'
#' Computes the ratio of account value to guaranteed maturity fund, capped at
#' one.
#'
#' @param AV Nonnegative account value.
#' @param GMF Positive guaranteed maturity fund.
#'
#' @return A numeric vector.
#'
#' @examples
#' rt_ul(AV = 4918.20, GMF = 14678.57)
#'
#' @export
rt_ul <- function(AV, GMF) {
  AV <- .ul_check_nonnegative(AV, "AV")
  GMF <- .ul_check_nonnegative(GMF, "GMF")

  values <- .ul_recycle(AV, GMF, .names = c("AV", "GMF"))
  AV <- values[[1]]
  GMF <- values[[2]]

  if (any(GMF <= 0)) {
    stop("GMF must contain positive values.", call. = FALSE)
  }

  pmin(AV / GMF, 1)
}

#' Pre-floor CRVM reserve
#'
#' Computes the pre-floor reserve by multiplying the funding ratio by the
#' difference between the present value of future benefits and future
#' premiums.
#'
#' @param r Funding ratio. Values must lie in \code{[0, 1]}.
#' @param pvfb_minus_pvfp Numeric difference between the present value of
#'   future benefits and future premiums.
#'
#' @return A numeric vector.
#'
#' @examples
#' Vprefloor_crvm_ul(r = 0.33506, pvfb_minus_pvfp = 70)
#'
#' @export
Vprefloor_crvm_ul <- function(r, pvfb_minus_pvfp) {
  r <- .ul_check_probability(r, "r")
  pvfb_minus_pvfp <- .ul_check_numeric(pvfb_minus_pvfp, "pvfb_minus_pvfp")

  values <- .ul_recycle(
    r, pvfb_minus_pvfp,
    .names = c("r", "pvfb_minus_pvfp")
  )

  values[[1]] * values[[2]]
}

#' AG 38 prefunding ratio
#'
#' Computes the excess-payment-to-required-net-single-premium ratio, capped at
#' one.
#'
#' @param excess_payment Nonnegative excess payment or shadow-fund amount.
#' @param nsp_required Positive net single premium required to fully fund the
#'   guarantee.
#'
#' @return A numeric vector.
#'
#' @examples
#' ag38_prefunding_ratio(60000, 100000)
#'
#' @export
ag38_prefunding_ratio <- function(excess_payment, nsp_required) {
  excess_payment <- .ul_check_nonnegative(excess_payment, "excess_payment")
  nsp_required <- .ul_check_nonnegative(nsp_required, "nsp_required")

  values <- .ul_recycle(
    excess_payment, nsp_required,
    .names = c("excess_payment", "nsp_required")
  )

  excess_payment <- values[[1]]
  nsp_required <- values[[2]]

  if (any(nsp_required <= 0)) {
    stop("nsp_required must contain positive values.", call. = FALSE)
  }

  pmin(excess_payment / nsp_required, 1)
}

#' AG 38 reserve calculation
#'
#' Computes the prefunding ratio, net additional amount, reduced deficiency
#' reserve, intermediate reserve, and increased basic reserve.
#'
#' Scalar inputs preserve the original named-list output. Vectorized inputs
#' return a data frame with one row per calculation.
#'
#' @param basic_reserve Nonnegative basic reserve.
#' @param deficiency_reserve Nonnegative deficiency reserve.
#' @param excess_payment Nonnegative excess payment or shadow-fund amount.
#' @param nsp_required Positive net single premium required to fully fund the
#'   guarantee.
#' @param valuation_nsp Nonnegative valuation net single premium.
#' @param surrender_charge Nonnegative surrender charge.
#'
#' @return For scalar inputs, a named list. For vectorized inputs, a data frame
#'   containing the same calculated quantities.
#'
#' @examples
#' ag38_reserve_ul(
#'   basic_reserve = 10000,
#'   deficiency_reserve = 0,
#'   excess_payment = 60000,
#'   nsp_required = 100000,
#'   valuation_nsp = 150000,
#'   surrender_charge = 5000
#' )
#'
#' @export
ag38_reserve_ul <- function(basic_reserve,
                            deficiency_reserve = 0,
                            excess_payment,
                            nsp_required,
                            valuation_nsp,
                            surrender_charge = 0) {
  basic_reserve <- .ul_check_nonnegative(basic_reserve, "basic_reserve")
  deficiency_reserve <- .ul_check_nonnegative(
    deficiency_reserve,
    "deficiency_reserve"
  )
  excess_payment <- .ul_check_nonnegative(excess_payment, "excess_payment")
  nsp_required <- .ul_check_nonnegative(nsp_required, "nsp_required")
  valuation_nsp <- .ul_check_nonnegative(valuation_nsp, "valuation_nsp")
  surrender_charge <- .ul_check_nonnegative(
    surrender_charge,
    "surrender_charge"
  )

  values <- .ul_recycle(
    basic_reserve,
    deficiency_reserve,
    excess_payment,
    nsp_required,
    valuation_nsp,
    surrender_charge,
    .names = c(
      "basic_reserve",
      "deficiency_reserve",
      "excess_payment",
      "nsp_required",
      "valuation_nsp",
      "surrender_charge"
    )
  )

  basic_reserve <- values[[1]]
  deficiency_reserve <- values[[2]]
  excess_payment <- values[[3]]
  nsp_required <- values[[4]]
  valuation_nsp <- values[[5]]
  surrender_charge <- values[[6]]

  if (any(nsp_required <= 0)) {
    stop("nsp_required must contain positive values.", call. = FALSE)
  }

  prefunding_ratio <- ag38_prefunding_ratio(excess_payment, nsp_required)
  base_plus_deficiency <- basic_reserve + deficiency_reserve

  net_amount_additional <-
    prefunding_ratio * (valuation_nsp - base_plus_deficiency)

  reduced_deficiency_reserve <-
    (1 - prefunding_ratio) * deficiency_reserve

  step8_reserve <-
    pmin(
      valuation_nsp,
      base_plus_deficiency + net_amount_additional
    ) - surrender_charge

  increased_basic_reserve <-
    step8_reserve - reduced_deficiency_reserve

  out <- data.frame(
    prefunding_ratio = prefunding_ratio,
    net_amount_additional = net_amount_additional,
    reduced_deficiency_reserve = reduced_deficiency_reserve,
    step8_reserve = step8_reserve,
    increased_basic_reserve = increased_basic_reserve,
    final_reserve = increased_basic_reserve
  )

  if (nrow(out) == 1L) {
    return(as.list(out[1L, , drop = FALSE]))
  }

  out
}
