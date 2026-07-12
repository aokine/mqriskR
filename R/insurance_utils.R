#' Insurance utilities
#'
#' Helper functions for insurance models.
#'
#' These utilities focus on:
#' \itemize{
#'   \item doubled-force calculations for second moments,
#'   \item UDD approximation multipliers,
#'   \item UDD approximations for continuous and m-thly insurance functions.
#' }
#'
#' General interest conversion functions such as \code{interest_convert()}
#' and \code{discount()} are defined elsewhere in the package and are reused here.
#'
#' @name insurance_utils
#' @keywords internal
NULL

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

.check_numeric_finite_ins_util <- function(x, name) {
  x <- as.numeric(x)

  if (length(x) == 0L || any(!is.finite(x))) {
    stop(name, " must contain finite numeric values.", call. = FALSE)
  }

  x
}

.check_i_ins_util <- function(i) {
  i <- .check_numeric_finite_ins_util(i, "i")

  if (any(i <= -1)) {
    stop("i must be greater than -1.", call. = FALSE)
  }

  i
}

.check_m_ins_util <- function(m) {
  m <- as.numeric(m)

  if (length(m) != 1L || !is.finite(m) || m <= 0 ||
      abs(m - round(m)) > 1e-10) {
    stop("m must be a positive integer.", call. = FALSE)
  }

  as.integer(round(m))
}

# ------------------------------------------------------------
# Doubled-force helpers for second moments
# ------------------------------------------------------------

#' Effective annual interest at doubled force
#'
#' If \eqn{\delta' = 2\delta}, then
#' \eqn{i' = (1+i)^2 - 1}.
#'
#' @param i Numeric vector of effective annual interest rates.
#'
#' @return Numeric vector of effective annual rates corresponding to doubled force.
#' @export
double_force_i <- function(i) {
  i <- .check_i_ins_util(i)
  (1 + i)^2 - 1
}

#' Doubled force of interest
#'
#' Computes \eqn{\delta' = 2\delta}.
#'
#' @param delta Numeric vector of forces of interest.
#'
#' @return Numeric vector of doubled forces of interest.
#' @export
double_force_delta <- function(delta) {
  delta <- .check_numeric_finite_ins_util(delta, "delta")
  2 * delta
}

# ------------------------------------------------------------
# UDD multipliers
# ------------------------------------------------------------

#' UDD multiplier for continuous insurance approximations
#'
#' Under UDD,
#' \eqn{\bar{A}_x = (i/\delta) A_x}
#' and similarly for term and deferred insurance.
#'
#' @param i Numeric vector of effective annual interest rates.
#'
#' @return Numeric vector equal to \eqn{i/\delta}.
#' @export
udd_continuous_multiplier <- function(i) {
  i <- .check_i_ins_util(i)

  delta <- log(1 + i)
  ifelse(abs(delta) < 1e-12, 1, i / delta)
}

#' UDD multiplier for m-thly insurance approximations
#'
#' Under UDD,
#' \eqn{A_x^{(m)} = (i / i^{(m)}) A_x}
#' and similarly for term and deferred insurance.
#'
#' @param i Numeric vector of effective annual interest rates.
#' @param m Positive integer payment frequency.
#'
#' @return Numeric vector equal to \eqn{i / i^{(m)}}.
#' @export
udd_mthly_multiplier <- function(i, m) {
  i <- .check_i_ins_util(i)
  m <- .check_m_ins_util(m)

  im <- m * ((1 + i)^(1 / m) - 1)
  ifelse(abs(im) < 1e-12, 1, i / im)
}

# ------------------------------------------------------------
# UDD approximations: continuous insurance
# ------------------------------------------------------------

#' UDD approximation of continuous whole life insurance
#'
#' Computes \eqn{\bar{A}_x = (i/\delta)A_x}.
#'
#' @param Ax Discrete whole life insurance APV.
#' @param i Effective annual interest rate.
#'
#' @return Continuous whole life insurance APV under UDD.
#' @export
Abarx_udd <- function(Ax, i) {
  Ax <- .check_numeric_finite_ins_util(Ax, "Ax")
  udd_continuous_multiplier(i) * Ax
}

#' UDD approximation of continuous term insurance
#'
#' Computes \eqn{\bar{A}_{x:\overline{n}|}^{1} = (i/\delta)A_{x:\overline{n}|}^{1}}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param i Effective annual interest rate.
#'
#' @return Continuous term insurance APV under UDD.
#' @export
Abarxn1_udd <- function(Axn1, i) {
  Axn1 <- .check_numeric_finite_ins_util(Axn1, "Axn1")
  udd_continuous_multiplier(i) * Axn1
}

#' UDD approximation of continuous deferred insurance
#'
#' Computes \eqn{{}_{n\mid}\bar{A}_x = (i/\delta)\,{}_{n\mid}A_x}.
#'
#' @param nAx Discrete deferred insurance APV.
#' @param i Effective annual interest rate.
#'
#' @return Continuous deferred insurance APV under UDD.
#' @export
nAbarx_udd <- function(nAx, i) {
  nAx <- .check_numeric_finite_ins_util(nAx, "nAx")
  udd_continuous_multiplier(i) * nAx
}

#' UDD approximation of continuous endowment insurance
#'
#' Computes
#' \eqn{\bar{A}_{x:\overline{n}|} = (i/\delta)A_{x:\overline{n}|}^{1} + {}_nE_x}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param nEx Pure endowment APV, \eqn{{}_nE_x}.
#' @param i Effective annual interest rate.
#'
#' @return Continuous endowment insurance APV under UDD.
#' @export
Abarxn_udd <- function(Axn1, nEx, i) {
  Axn1 <- .check_numeric_finite_ins_util(Axn1, "Axn1")
  nEx <- .check_numeric_finite_ins_util(nEx, "nEx")

  Abarxn1_udd(Axn1, i) + nEx
}

# ------------------------------------------------------------
# UDD approximations: m-thly insurance
# ------------------------------------------------------------

#' UDD approximation of m-thly whole life insurance
#'
#' Computes \eqn{A_x^{(m)} = (i/i^{(m)})A_x}.
#'
#' @param Ax Discrete whole life insurance APV.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#'
#' @return m-thly whole life insurance APV under UDD.
#' @export
Ax_m_udd <- function(Ax, i, m) {
  Ax <- .check_numeric_finite_ins_util(Ax, "Ax")
  udd_mthly_multiplier(i, m) * Ax
}

#' UDD approximation of m-thly term insurance
#'
#' Computes \eqn{A_{x:\overline{n}|}^{1(m)} = (i/i^{(m)})A_{x:\overline{n}|}^{1}}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#'
#' @return m-thly term insurance APV under UDD.
#' @export
Axn1_m_udd <- function(Axn1, i, m) {
  Axn1 <- .check_numeric_finite_ins_util(Axn1, "Axn1")
  udd_mthly_multiplier(i, m) * Axn1
}

#' UDD approximation of m-thly deferred insurance
#'
#' Computes \eqn{{}_{n\mid}A_x^{(m)} = (i/i^{(m)})\,{}_{n\mid}A_x}.
#'
#' @param nAx Discrete deferred insurance APV.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#'
#' @return m-thly deferred insurance APV under UDD.
#' @export
nAx_m_udd <- function(nAx, i, m) {
  nAx <- .check_numeric_finite_ins_util(nAx, "nAx")
  udd_mthly_multiplier(i, m) * nAx
}

#' UDD approximation of m-thly endowment insurance
#'
#' Computes
#' \eqn{A_{x:\overline{n}|}^{(m)} = (i/i^{(m)})A_{x:\overline{n}|}^{1} + {}_nE_x}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param nEx Pure endowment APV.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#'
#' @return m-thly endowment insurance APV under UDD.
#' @export
Axn_m_udd <- function(Axn1, nEx, i, m) {
  Axn1 <- .check_numeric_finite_ins_util(Axn1, "Axn1")
  nEx <- .check_numeric_finite_ins_util(nEx, "nEx")

  Axn1_m_udd(Axn1, i, m) + nEx
}
