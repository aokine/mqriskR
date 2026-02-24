#' Internal: safe numerical integration helper
#'
#' @keywords internal
.safe_integrate <- function(f, lower, upper, ..., rel.tol = 1e-10) {
  out <- stats::integrate(f, lower = lower, upper = upper, ..., rel.tol = rel.tol)
  unname(out$value)
}

#' Internal: choose a practical upper bound for integrals to infinity
#'
#' Strategy: find t such that S0(t) is very small, or fall back to a cap.
#'
#' @keywords internal
.find_upper_t <- function(x, model, ..., eps = 1e-10, t_max = 500) {
  # we want S0(x + t) / S0(x) ~ 0, equivalently S0(x+t) small
  # Start with a modest grid and expand until survival is tiny.
  grid <- c(10, 20, 50, 100, 200, 300, 500)
  for (t in grid) {
    s <- S0(x + t, model = model, ...)
    if (is.finite(s) && s < eps) return(t)
  }
  t_max
}
