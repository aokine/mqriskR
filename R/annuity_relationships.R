#' Annuity-insurance relationships
#'
#' Identities linking annual and continuous annuity functions to the
#' corresponding insurance functions.
#'
#' @param x Age.
#' @param n Term or deferral period in years.
#' @param i Effective annual interest rate.
#' @param tbl Optional life table object for discrete identities.
#' @param model Optional survival model name.
#' @param ... Additional model parameters.
#' @param k_max Maximum summation horizon for non-terminating models.
#' @param tol Truncation tolerance for non-terminating models.
#'
#' @return
#' Numeric vector containing the annuity value computed from the corresponding
#' annuity-insurance identity.
#'
#' @name annuity_relationships
NULL

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_rel_i <- function(i) {
  i <- as.numeric(i)

  if (length(i) == 0L || any(!is.finite(i)) || any(i <= -1)) {
    stop("i must contain finite values greater than -1.", call. = FALSE)
  }

  i
}

.check_rel_n <- function(n) {
  n <- as.numeric(n)

  if (length(n) != 1L || !is.finite(n) || n < 0) {
    stop("n must be a single nonnegative number.", call. = FALSE)
  }

  n
}

.check_rel_integer_n <- function(n) {
  n <- .check_rel_n(n)

  if (abs(n - round(n)) > 1e-10) {
    stop("For annual annuity relationships, n must be an integer.", call. = FALSE)
  }

  as.integer(round(n))
}

.recycle_rel_vectors <- function(...) {
  args <- lapply(list(...), as.numeric)
  lens <- vapply(args, length, integer(1))

  if (any(lens == 0L)) {
    stop("All arguments must have positive length.", call. = FALSE)
  }

  if (any(vapply(args, function(z) any(!is.finite(z)), logical(1)))) {
    stop("All arguments must contain finite numeric values.", call. = FALSE)
  }

  common_len <- max(lens)

  if (any(!(lens %in% c(1L, common_len)))) {
    stop("Arguments must have compatible lengths.", call. = FALSE)
  }

  lapply(args, function(z) if (length(z) == 1L) rep(z, common_len) else z)
}

.check_rel_x <- function(x) {
  if (any(x < 0)) {
    stop("x must be nonnegative.", call. = FALSE)
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Whole life identities
# -------------------------------------------------------------------------

#' Whole life annuity-immediate from insurance identity
#'
#' Computes \eqn{a_x = (v - A_x)/d}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_ax <- function(x, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_rel_vectors(x, .check_rel_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    d <- i[j] / (1 + i[j])
    v <- 1 / (1 + i[j])

    (v - Ax(x[j], i[j], tbl = tbl, model = model, ...)) / d
  })
}

#' Whole life annuity-due from insurance identity
#'
#' Computes \eqn{\ddot{a}_x = (1 - A_x)/d}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_adotx <- function(x, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_rel_vectors(x, .check_rel_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    d <- i[j] / (1 + i[j])

    (1 - Ax(x[j], i[j], tbl = tbl, model = model, ...)) / d
  })
}

#' Continuous whole life annuity from insurance identity
#'
#' Computes \eqn{\bar{a}_x = (1 - \bar{A}_x)/\delta}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_abarx <- function(x, i, model, ...) {
  args <- .recycle_rel_vectors(x, .check_rel_i(i))
  x <- args[[1L]]
  i <- args[[2L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    delta <- log(1 + i[j])

    (1 - Abarx(x[j], i[j], model = model, ...)) / delta
  })
}

# -------------------------------------------------------------------------
# Temporary identities
# -------------------------------------------------------------------------

#' Temporary annuity-immediate from insurance identity
#'
#' Computes
#' \eqn{a_{x:\overline{n}|} = (1 - A_{x:\overline{n}|})/d - 1 + {}_nE_x}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_axn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_rel_vectors(x, n, .check_rel_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_integer_n(n[j])
    d <- i[j] / (1 + i[j])

    (1 - Axn(x[j], nn, i[j], tbl = tbl, model = model, ...)) / d - 1 +
      nEx(x[j], nn, i[j], tbl = tbl, model = model, ...)
  })
}

#' Temporary annuity-due from insurance identity
#'
#' Computes \eqn{\ddot{a}_{x:\overline{n}|} =
#' (1 - A_{x:\overline{n}|})/d}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_adotxn <- function(x, n, i, model = NULL, ..., tbl = NULL) {
  args <- .recycle_rel_vectors(x, n, .check_rel_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_integer_n(n[j])
    d <- i[j] / (1 + i[j])

    (1 - Axn(x[j], nn, i[j], tbl = tbl, model = model, ...)) / d
  })
}

#' Continuous temporary annuity from insurance identity
#'
#' Computes \eqn{\bar{a}_{x:\overline{n}|} =
#' (1 - \bar{A}_{x:\overline{n}|})/\delta}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_abarxn <- function(x, n, i, model, ...) {
  args <- .recycle_rel_vectors(x, n, .check_rel_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])
    delta <- log(1 + i[j])

    (1 - Abarxn(x[j], nn, i[j], model = model, ...)) / delta
  })
}

# -------------------------------------------------------------------------
# Deferred identities
# -------------------------------------------------------------------------

#' Deferred annuity-immediate from pure endowment identity
#'
#' Computes \eqn{{}_{n|}a_x = {}_nE_x a_{x+n}}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_nax <- function(x, n, i, model = NULL, ..., tbl = NULL,
                                 k_max = 5000, tol = 1e-12) {
  args <- .recycle_rel_vectors(x, n, .check_rel_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_integer_n(n[j])

    nEx(x[j], nn, i[j], tbl = tbl, model = model, ...) *
      ax(
        x[j] + nn,
        i[j],
        tbl = tbl,
        model = model,
        ...,
        k_max = k_max,
        tol = tol
      )
  })
}

#' Deferred annuity-due from pure endowment identity
#'
#' Computes \eqn{{}_{n|}\ddot{a}_x = {}_nE_x \ddot{a}_{x+n}}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_nadotx <- function(x, n, i, model = NULL, ..., tbl = NULL,
                                    k_max = 5000, tol = 1e-12) {
  args <- .recycle_rel_vectors(x, n, .check_rel_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_integer_n(n[j])

    nEx(x[j], nn, i[j], tbl = tbl, model = model, ...) *
      adotx(
        x[j] + nn,
        i[j],
        tbl = tbl,
        model = model,
        ...,
        k_max = k_max,
        tol = tol
      )
  })
}

#' Deferred continuous annuity from pure endowment identity
#'
#' Computes \eqn{{}_{n|}\bar{a}_x = {}_nE_x \bar{a}_{x+n}}.
#'
#' @rdname annuity_relationships
#' @export
annuity_identity_nabarx <- function(x, n, i, model, ..., tol = 1e-10) {
  args <- .recycle_rel_vectors(x, n, .check_rel_i(i))
  x <- args[[1L]]
  n <- args[[2L]]
  i <- args[[3L]]

  .check_rel_x(x)

  sapply(seq_along(x), function(j) {
    nn <- .check_rel_n(n[j])

    nEx(x[j], nn, i[j], model = model, ...) *
      abarx(x[j] + nn, i[j], model = model, ..., tol = tol)
  })
}
