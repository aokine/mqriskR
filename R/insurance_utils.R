#' Insurance utilities (Chapter 7)
#'
#' Helper functions for Chapter 7 insurance models.
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
  i <- as.numeric(i)
  if (any(!is.finite(i))) stop("i must be finite.", call. = FALSE)
  if (any(i <= -1)) stop("i must be > -1.", call. = FALSE)
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
  delta <- as.numeric(delta)
  if (any(!is.finite(delta))) stop("delta must be finite.", call. = FALSE)
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
  i <- as.numeric(i)
  if (any(!is.finite(i))) stop("i must be finite.", call. = FALSE)
  if (any(i <= -1)) stop("i must be > -1.", call. = FALSE)

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
  i <- as.numeric(i)
  if (any(!is.finite(i))) stop("i must be finite.", call. = FALSE)
  if (any(i <= -1)) stop("i must be > -1.", call. = FALSE)
  if (length(m) != 1 || !is.finite(m) || m <= 0 || m != as.integer(m)) {
    stop("m must be a positive integer.", call. = FALSE)
  }

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
  Ax <- as.numeric(Ax)
  if (any(!is.finite(Ax))) stop("Ax must be finite.", call. = FALSE)
  udd_continuous_multiplier(i) * Ax
}

#' UDD approximation of continuous term insurance
#'
#' Computes \eqn{\bar{A}_{x:\angl{n}}^{1} = (i/\delta)A_{x:\angl{n}}^{1}}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param i Effective annual interest rate.
#'
#' @return Continuous term insurance APV under UDD.
#' @export
Abarxn1_udd <- function(Axn1, i) {
  Axn1 <- as.numeric(Axn1)
  if (any(!is.finite(Axn1))) stop("Axn1 must be finite.", call. = FALSE)
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
  nAx <- as.numeric(nAx)
  if (any(!is.finite(nAx))) stop("nAx must be finite.", call. = FALSE)
  udd_continuous_multiplier(i) * nAx
}

#' UDD approximation of continuous endowment insurance
#'
#' Computes
#' \eqn{\bar{A}_{x:\angl{n}} = (i/\delta)A_{x:\angl{n}}^{1} + {}_nE_x}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param nEx Pure endowment APV, \eqn{{}_nE_x}.
#' @param i Effective annual interest rate.
#'
#' @return Continuous endowment insurance APV under UDD.
#' @export
Abarxn_udd <- function(Axn1, nEx, i) {
  Axn1 <- as.numeric(Axn1)
  nEx <- as.numeric(nEx)

  if (any(!is.finite(Axn1))) stop("Axn1 must be finite.", call. = FALSE)
  if (any(!is.finite(nEx))) stop("nEx must be finite.", call. = FALSE)

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
  Ax <- as.numeric(Ax)
  if (any(!is.finite(Ax))) stop("Ax must be finite.", call. = FALSE)
  udd_mthly_multiplier(i, m) * Ax
}

#' UDD approximation of m-thly term insurance
#'
#' Computes \eqn{A_{x:\angl{n}}^{1(m)} = (i/i^{(m)})A_{x:\angl{n}}^{1}}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#'
#' @return m-thly term insurance APV under UDD.
#' @export
Axn1_m_udd <- function(Axn1, i, m) {
  Axn1 <- as.numeric(Axn1)
  if (any(!is.finite(Axn1))) stop("Axn1 must be finite.", call. = FALSE)
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
  nAx <- as.numeric(nAx)
  if (any(!is.finite(nAx))) stop("nAx must be finite.", call. = FALSE)
  udd_mthly_multiplier(i, m) * nAx
}

#' UDD approximation of m-thly endowment insurance
#'
#' Computes
#' \eqn{A_{x:\angl{n}}^{(m)} = (i/i^{(m)})A_{x:\angl{n}}^{1} + {}_nE_x}.
#'
#' @param Axn1 Discrete term insurance APV.
#' @param nEx Pure endowment APV.
#' @param i Effective annual interest rate.
#' @param m Positive integer payment frequency.
#'
#' @return m-thly endowment insurance APV under UDD.
#' @export
Axn_m_udd <- function(Axn1, nEx, i, m) {
  Axn1 <- as.numeric(Axn1)
  nEx <- as.numeric(nEx)

  if (any(!is.finite(Axn1))) stop("Axn1 must be finite.", call. = FALSE)
  if (any(!is.finite(nEx))) stop("nEx must be finite.", call. = FALSE)

  Axn1_m_udd(Axn1, i, m) + nEx
}
