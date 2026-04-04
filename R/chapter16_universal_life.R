# =========================================================
# Chapter 16 helpers for mqriskR
# Universal Life Insurance
# =========================================================

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.validate_prob_vector_ch16 <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite numeric vector.", name), call. = FALSE)
  }
  if (any(x < 0 | x > 1)) {
    stop(sprintf("%s must contain values in [0, 1].", name), call. = FALSE)
  }
  invisible(x)
}

.validate_nonneg_numeric_ch16 <- function(x, name) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite numeric vector.", name), call. = FALSE)
  }
  if (any(x < 0)) {
    stop(sprintf("%s must be nonnegative.", name), call. = FALSE)
  }
  invisible(x)
}

.recycle_ch16 <- function(x, n, name) {
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(sprintf("%s must have length 1 or %d.", name, n), call. = FALSE)
  }
  x
}

# -------------------------------------------------------------------------
# Cost of insurance
# -------------------------------------------------------------------------

#' Cost of insurance for Type B universal life
#'
#' Computes the one-period cost of insurance for a Type B universal life policy
#' using Equation (16.5a):
#' \deqn{
#' COI_t = \frac{B q_{x+t-1}}{1+i^q}.
#' }
#'
#' @param B Face amount.
#' @param qx Mortality rate for the period.
#' @param iq Interest rate used in the cost-of-insurance calculation.
#'
#' @return Numeric vector.
#'
#' @examples
#' coi_ul_typeB(B = 100000, qx = 0.00076, iq = 0.03)
#'
#' @export
coi_ul_typeB <- function(B, qx, iq) {
  .validate_nonneg_numeric_ch16(B, "B")
  .validate_prob_vector_ch16(qx, "qx")
  .validate_nonneg_numeric_ch16(iq, "iq")

  vals_n <- max(length(B), length(qx), length(iq))
  B <- .recycle_ch16(B, vals_n, "B")
  qx <- .recycle_ch16(qx, vals_n, "qx")
  iq <- .recycle_ch16(iq, vals_n, "iq")

  B * qx / (1 + iq)
}

# -------------------------------------------------------------------------
# Universal life account value roll-forward
# -------------------------------------------------------------------------

#' Account-value path for Type B universal life
#'
#' Computes the year-by-year account value roll-forward for Type B universal
#' life using Equations (16.4) and (16.5a).
#'
#' @param G Premium vector \eqn{G_t}.
#' @param r Percent-of-premium expense vector \eqn{r_t}.
#' @param e Fixed expense vector \eqn{e_t}.
#' @param qx Mortality vector \eqn{q_{x+t-1}}.
#' @param ic Credited interest rate vector \eqn{i^c}.
#' @param B Face amount.
#' @param iq Interest rate vector \eqn{i^q} used in cost of insurance.
#'   Defaults to \code{ic}.
#' @param AV0 Initial account value. Defaults to 0.
#'
#' @return A data frame with columns \code{t}, \code{premium},
#'   \code{net_contribution}, \code{COI}, and \code{AV}.
#'
#' @examples
#' qx <- c(.00076, .00081, .00085, .00090, .00095)
#' r <- c(.75, .10, .10, .10, .10)
#' e <- c(100, 20, 20, 20, 20)
#' G <- rep(5000, 5)
#'
#' AV_path_ul_typeB(G = G, r = r, e = e, qx = qx, ic = 0.03, B = 100000)
#'
#' @export
AV_path_ul_typeB <- function(G, r, e, qx, ic, B, iq = ic, AV0 = 0) {
  .validate_nonneg_numeric_ch16(G, "G")
  .validate_nonneg_numeric_ch16(r, "r")
  .validate_nonneg_numeric_ch16(e, "e")
  .validate_prob_vector_ch16(qx, "qx")
  .validate_nonneg_numeric_ch16(ic, "ic")
  .validate_nonneg_numeric_ch16(iq, "iq")
  .validate_nonneg_numeric_ch16(B, "B")
  .validate_nonneg_numeric_ch16(AV0, "AV0")

  n <- length(qx)
  G  <- .recycle_ch16(G,  n, "G")
  r  <- .recycle_ch16(r,  n, "r")
  e  <- .recycle_ch16(e,  n, "e")
  ic <- .recycle_ch16(ic, n, "ic")
  iq <- .recycle_ch16(iq, n, "iq")

  AV <- numeric(n + 1L)
  net_contribution <- numeric(n)
  COI <- numeric(n)
  AV[1L] <- AV0

  for (k in seq_len(n)) {
    net_contribution[k] <- G[k] * (1 - r[k]) - e[k]
    COI[k] <- coi_ul_typeB(B = B, qx = qx[k], iq = iq[k])
    AV[k + 1L] <- (AV[k] + net_contribution[k] - COI[k]) * (1 + ic[k])
  }

  data.frame(
    t = 0:n,
    premium = c(NA_real_, G),
    net_contribution = c(NA_real_, net_contribution),
    COI = c(NA_real_, COI),
    AV = AV
  )
}

#' Account-value path for Type A universal life
#'
#' Computes the year-by-year account value roll-forward for Type A universal
#' life using the explicit form in Equation (16.8), or Equation (16.9) when
#' \eqn{i^q = i^c}.
#'
#' @param G Premium vector \eqn{G_t}.
#' @param r Percent-of-premium expense vector \eqn{r_t}.
#' @param e Fixed expense vector \eqn{e_t}.
#' @param qx Mortality vector \eqn{q_{x+t-1}}.
#' @param ic Credited interest rate vector \eqn{i^c}.
#' @param B Fixed death benefit face amount.
#' @param iq Interest rate vector \eqn{i^q} used in cost of insurance.
#'   Defaults to \code{ic}.
#' @param AV0 Initial account value. Defaults to 0.
#'
#' @return A data frame with columns \code{t}, \code{premium}, \code{AV}.
#'
#' @examples
#' qx <- c(.00076, .00081, .00085, .00090, .00095)
#' r <- c(.75, .10, .10, .10, .10)
#' e <- c(100, 20, 20, 20, 20)
#' G <- rep(5000, 5)
#'
#' AV_path_ul_typeA(G = G, r = r, e = e, qx = qx, ic = 0.03, B = 100000)
#'
#' @export
AV_path_ul_typeA <- function(G, r, e, qx, ic, B, iq = ic, AV0 = 0) {
  .validate_nonneg_numeric_ch16(G, "G")
  .validate_nonneg_numeric_ch16(r, "r")
  .validate_nonneg_numeric_ch16(e, "e")
  .validate_prob_vector_ch16(qx, "qx")
  .validate_nonneg_numeric_ch16(ic, "ic")
  .validate_nonneg_numeric_ch16(iq, "iq")
  .validate_nonneg_numeric_ch16(B, "B")
  .validate_nonneg_numeric_ch16(AV0, "AV0")

  n <- length(qx)
  G  <- .recycle_ch16(G,  n, "G")
  r  <- .recycle_ch16(r,  n, "r")
  e  <- .recycle_ch16(e,  n, "e")
  ic <- .recycle_ch16(ic, n, "ic")
  iq <- .recycle_ch16(iq, n, "iq")

  AV <- numeric(n + 1L)
  AV[1L] <- AV0

  for (k in seq_len(n)) {
    num <- (AV[k] + G[k] * (1 - r[k]) - e[k]) * (1 + ic[k]) -
      B * qx[k] * (1 + ic[k]) / (1 + iq[k])

    den <- 1 - qx[k] * (1 + ic[k]) / (1 + iq[k])

    if (den <= 0) {
      stop("Type A roll-forward denominator is nonpositive.", call. = FALSE)
    }

    AV[k + 1L] <- num / den
  }

  data.frame(
    t = 0:n,
    premium = c(NA_real_, G),
    AV = AV
  )
}

# -------------------------------------------------------------------------
# EIUL credited-rate helpers
# -------------------------------------------------------------------------

#' Point-to-point index growth rates
#'
#' Computes annual point-to-point index growth rates from successive index
#' values, matching Equation (16.12).
#'
#' @param index Numeric vector of index closing values.
#'
#' @return Numeric vector of growth rates.
#'
#' @examples
#' iP_eiul(c(1000, 1050, 1200, 1100))
#'
#' @export
iP_eiul <- function(index) {
  .validate_nonneg_numeric_ch16(index, "index")

  if (length(index) < 2L) {
    stop("index must contain at least two values.", call. = FALSE)
  }

  index[-1L] / index[-length(index)] - 1
}

#' Monthly-average index growth rate
#'
#' Computes the monthly-average index growth rate from an initial index value
#' and 12 monthly closing values, matching Equation (16.13).
#'
#' @param index Numeric vector of length 13 containing the initial index value
#'   followed by the 12 monthly closing values.
#'
#' @return Numeric scalar.
#'
#' @examples
#' idx <- c(1000, 1020, 1100, 1150, 1080, 1040, 960, 1030, 1000, 1070, 1150, 1200, 1150)
#' iMA_eiul(idx)
#'
#' @export
iMA_eiul <- function(index) {
  .validate_nonneg_numeric_ch16(index, "index")

  if (length(index) != 13L) {
    stop("index must contain 13 values: initial plus 12 monthly closings.", call. = FALSE)
  }

  mean(index[-1L]) / index[1L] - 1
}

#' Credited rates from raw index growth rates
#'
#' Applies participation, floor, cap, and optional margin to raw index growth
#' rates for an indexed universal life contract.
#'
#' @param i_raw Numeric vector of raw index growth rates.
#' @param part Participation rate.
#' @param floor Minimum credited rate.
#' @param cap Maximum credited rate.
#' @param margin Index margin. Defaults to 0.
#' @param margin_after_participation Logical; if \code{TRUE}, subtract the
#'   margin after applying the participation rate.
#'
#' @return Numeric vector of credited rates.
#'
#' @examples
#' raw <- iP_eiul(c(1000, 1050, 1200, 1100))
#' i_credit_eiul(raw, part = 1.10, floor = 0.01, cap = 0.10)
#'
#' @export
i_credit_eiul <- function(i_raw, part = 1, floor = 0, cap = Inf,
                          margin = 0, margin_after_participation = TRUE) {
  if (!is.numeric(i_raw) || any(!is.finite(i_raw))) {
    stop("i_raw must be a finite numeric vector.", call. = FALSE)
  }
  .validate_nonneg_numeric_ch16(part, "part")
  .validate_nonneg_numeric_ch16(floor, "floor")
  .validate_nonneg_numeric_ch16(cap, "cap")
  .validate_nonneg_numeric_ch16(margin, "margin")

  if (length(part) != 1L || length(floor) != 1L || length(cap) != 1L || length(margin) != 1L) {
    stop("part, floor, cap, and margin must be numeric scalars.", call. = FALSE)
  }
  if (floor > cap) {
    stop("floor cannot exceed cap.", call. = FALSE)
  }

  adj <- if (margin_after_participation) {
    part * i_raw - margin
  } else {
    part * (i_raw - margin)
  }

  pmin(cap, pmax(floor, adj))
}

# -------------------------------------------------------------------------
# Persistency / survival under mortality and withdrawal
# -------------------------------------------------------------------------

#' One-year persistency rates for universal life
#'
#' Computes the one-year persistency rates \eqn{p_x^{(\tau)}} under either
#' Equation (16.14) or Equation (16.15).
#'
#' @param qd Mortality probabilities.
#' @param qw Withdrawal probabilities.
#' @param year_end_withdrawal Logical; if \code{TRUE}, use
#'   \eqn{(1-q^{(d)})(1-q^{(w)})}. Otherwise use \eqn{1-q^{(d)}-q^{(w)}}.
#'
#' @return Numeric vector.
#'
#' @examples
#' qd <- c(.001, .002, .003)
#' qw <- c(.02, .02, .03)
#' pxtau_ul(qd, qw)
#'
#' @export
pxtau_ul <- function(qd, qw, year_end_withdrawal = TRUE) {
  .validate_prob_vector_ch16(qd, "qd")
  .validate_prob_vector_ch16(qw, "qw")

  if (length(qd) != length(qw)) {
    stop("qd and qw must have the same length.", call. = FALSE)
  }

  if (year_end_withdrawal) {
    (1 - qd) * (1 - qw)
  } else {
    out <- 1 - qd - qw
    if (any(out < 0)) {
      stop("Some values of 1 - qd - qw are negative.", call. = FALSE)
    }
    out
  }
}

#' Cumulative persistency to the end of each policy year
#'
#' Computes \eqn{{}_tp_x^{(\tau)}} from the one-year persistency rates.
#'
#' @param qd Mortality probabilities.
#' @param qw Withdrawal probabilities.
#' @param year_end_withdrawal Logical; if \code{TRUE}, use Equation (16.15).
#'
#' @return Numeric vector.
#'
#' @examples
#' qd <- c(.001, .002, .003)
#' qw <- c(.02, .02, .03)
#' tpxtau_ul(qd, qw)
#'
#' @export
tpxtau_ul <- function(qd, qw, year_end_withdrawal = TRUE) {
  cumprod(pxtau_ul(qd = qd, qw = qw, year_end_withdrawal = year_end_withdrawal))
}

# -------------------------------------------------------------------------
# UL reserving helpers
# -------------------------------------------------------------------------

#' Guaranteed maturity fund roll-forward
#'
#' Computes the one-period guaranteed maturity fund roll-forward used in
#' Example 16.9.
#'
#' @param GMF_prev Prior guaranteed maturity fund.
#' @param GMP Guaranteed maturity premium.
#' @param r Expense factor applied to GMP.
#' @param policy_charge Guaranteed policy charge.
#' @param i Guaranteed interest rate.
#'
#' @return Numeric scalar.
#'
#' @examples
#' GMF_rollforward_ul(140.40, 14.49, 0.04, 11.80, 0.03)
#'
#' @export
GMF_rollforward_ul <- function(GMF_prev, GMP, r, policy_charge, i) {
  .validate_nonneg_numeric_ch16(GMF_prev, "GMF_prev")
  .validate_nonneg_numeric_ch16(GMP, "GMP")
  .validate_nonneg_numeric_ch16(r, "r")
  .validate_nonneg_numeric_ch16(policy_charge, "policy_charge")
  .validate_nonneg_numeric_ch16(i, "i")

  if (any(lengths(list(GMF_prev, GMP, r, policy_charge, i)) != 1L)) {
    stop("All arguments must be numeric scalars.", call. = FALSE)
  }

  (GMF_prev + GMP * (1 - r) - policy_charge) * (1 + i)
}

#' Ratio \eqn{r_t = AV_t / GMF_t} capped at 1
#'
#' Computes the ratio in Equation (16.16).
#'
#' @param AV Account value.
#' @param GMF Guaranteed maturity fund.
#'
#' @return Numeric scalar.
#'
#' @examples
#' rt_ul(AV = 4918.20, GMF = 14678.57)
#'
#' @export
rt_ul <- function(AV, GMF) {
  .validate_nonneg_numeric_ch16(AV, "AV")
  .validate_nonneg_numeric_ch16(GMF, "GMF")

  if (any(lengths(list(AV, GMF)) != 1L)) {
    stop("AV and GMF must be numeric scalars.", call. = FALSE)
  }
  if (GMF == 0) {
    stop("GMF must be positive.", call. = FALSE)
  }

  min(AV / GMF, 1)
}

#' Pre-floor CRVM reserve for universal life
#'
#' Computes the pre-floor CRVM reserve from Equation (16.17).
#'
#' @param r Ratio \eqn{r_t}.
#' @param pvfb_minus_pvfp Difference \eqn{(PVFB)_t - (PVFP)_t}.
#'
#' @return Numeric scalar.
#'
#' @examples
#' Vprefloor_crvm_ul(r = 0.33506, pvfb_minus_pvfp = 70)
#'
#' @export
Vprefloor_crvm_ul <- function(r, pvfb_minus_pvfp) {
  .validate_nonneg_numeric_ch16(r, "r")
  .validate_nonneg_numeric_ch16(pvfb_minus_pvfp, "pvfb_minus_pvfp")

  if (length(r) != 1L || length(pvfb_minus_pvfp) != 1L) {
    stop("r and pvfb_minus_pvfp must be numeric scalars.", call. = FALSE)
  }

  r * pvfb_minus_pvfp
}

#' AG 38 prefunding ratio
#'
#' Computes the prefunding ratio in Equation (16.18), capped at 1.
#'
#' @param excess_payment Excess payment or shadow-fund amount.
#' @param nsp_required Net single premium required to fully fund the guarantee.
#'
#' @return Numeric scalar.
#'
#' @examples
#' ag38_prefunding_ratio(60000, 100000)
#'
#' @export
ag38_prefunding_ratio <- function(excess_payment, nsp_required) {
  .validate_nonneg_numeric_ch16(excess_payment, "excess_payment")
  .validate_nonneg_numeric_ch16(nsp_required, "nsp_required")

  if (length(excess_payment) != 1L || length(nsp_required) != 1L) {
    stop("excess_payment and nsp_required must be numeric scalars.", call. = FALSE)
  }
  if (nsp_required == 0) {
    stop("nsp_required must be positive.", call. = FALSE)
  }

  min(excess_payment / nsp_required, 1)
}

#' AG 38 reserve calculation
#'
#' Computes the main quantities in the Chapter 16 AG 38 reserve calculation,
#' including the prefunding ratio, reduced deficiency reserve, Step (8)
#' reserve, and final increased basic reserve.
#'
#' @param basic_reserve Basic reserve.
#' @param deficiency_reserve Deficiency reserve.
#' @param excess_payment Excess payment or shadow-fund amount.
#' @param nsp_required Net single premium required to fully fund the guarantee.
#' @param valuation_nsp Valuation net single premium.
#' @param surrender_charge Applicable surrender charge.
#'
#' @return A named list.
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
  .validate_nonneg_numeric_ch16(basic_reserve, "basic_reserve")
  .validate_nonneg_numeric_ch16(deficiency_reserve, "deficiency_reserve")
  .validate_nonneg_numeric_ch16(excess_payment, "excess_payment")
  .validate_nonneg_numeric_ch16(nsp_required, "nsp_required")
  .validate_nonneg_numeric_ch16(valuation_nsp, "valuation_nsp")
  .validate_nonneg_numeric_ch16(surrender_charge, "surrender_charge")

  args <- list(
    basic_reserve = basic_reserve,
    deficiency_reserve = deficiency_reserve,
    excess_payment = excess_payment,
    nsp_required = nsp_required,
    valuation_nsp = valuation_nsp,
    surrender_charge = surrender_charge
  )
  if (any(vapply(args, length, integer(1)) != 1L)) {
    stop("All arguments must be numeric scalars.", call. = FALSE)
  }

  prefunding_ratio <- ag38_prefunding_ratio(
    excess_payment = excess_payment,
    nsp_required = nsp_required
  )

  base_plus_def <- basic_reserve + deficiency_reserve
  net_amount_additional <- prefunding_ratio * (valuation_nsp - base_plus_def)
  reduced_deficiency_reserve <- (1 - prefunding_ratio) * deficiency_reserve

  step8_reserve <- min(
    valuation_nsp,
    base_plus_def + net_amount_additional
  ) - surrender_charge

  increased_basic_reserve <- step8_reserve - reduced_deficiency_reserve

  list(
    prefunding_ratio = prefunding_ratio,
    net_amount_additional = net_amount_additional,
    reduced_deficiency_reserve = reduced_deficiency_reserve,
    step8_reserve = step8_reserve,
    increased_basic_reserve = increased_basic_reserve,
    final_reserve = increased_basic_reserve
  )
}
