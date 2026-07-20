# -------------------------------------------------------------------------
# Internal validation and recycling helpers
# -------------------------------------------------------------------------

#' @noRd
.multilife_ext_check_basis <- function(tbl = NULL,
                                       model = NULL,
                                       require_continuous = TRUE) {
  if (is.null(tbl) && is.null(model)) {
    stop("Supply either tbl or model.", call. = FALSE)
  }

  if (!is.null(tbl) && !is.null(model)) {
    stop("Supply only one of tbl or model, not both.", call. = FALSE)
  }

  if (!is.null(tbl)) {
    .validate_life_table(tbl)

    if (isTRUE(require_continuous)) {
      stop(
        paste0(
          "Continuous multi-life calculations require a parametric model. ",
          "An annual life table does not determine within-year survival ",
          "or forces of mortality without an interpolation assumption."
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(model)) {
    if (length(model) != 1L ||
        !is.character(model) ||
        is.na(model) ||
        !nzchar(model)) {
      stop("model must be a single nonempty character string.",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}


#' @noRd
.multilife_ext_check_nonnegative <- function(x, name) {
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
.multilife_ext_check_interest <- function(i) {
  i <- as.numeric(i)

  if (length(i) == 0L) {
    stop("i must have positive length.", call. = FALSE)
  }

  if (any(!is.finite(i)) || any(i <= -1)) {
    stop(
      "i must contain finite values greater than -1.",
      call. = FALSE
    )
  }

  i
}


#' @noRd
.multilife_ext_recycle <- function(..., .names = NULL) {
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
.multilife_ext_upper <- function(x,
                                 y,
                                 n = NULL,
                                 model,
                                 ...) {
  if (tolower(model) == "uniform") {
    dots <- list(...)

    if (is.null(dots$omega) ||
        length(dots$omega) != 1L ||
        !is.finite(dots$omega)) {
      stop(
        "For model = 'uniform', omega must be supplied as a finite scalar.",
        call. = FALSE
      )
    }

    upper <- min(dots$omega - x, dots$omega - y)

    if (!is.null(n)) {
      upper <- min(upper, n)
    }

    return(max(upper, 0))
  }

  if (!is.null(n)) {
    return(n)
  }

  Inf
}


#' @noRd
.multilife_ext_integrate <- function(f, upper) {
  if (!is.finite(upper) && upper != Inf) {
    stop("The integration limit must be finite or Inf.", call. = FALSE)
  }

  if (upper <= 0) {
    return(0)
  }

  value <- stats::integrate(
    f,
    lower = 0,
    upper = upper,
    rel.tol = 1e-9,
    subdivisions = 1000L,
    stop.on.error = TRUE
  )$value

  if (!is.finite(value)) {
    stop("The continuous multi-life integral did not return a finite value.",
         call. = FALSE)
  }

  value
}


#' @noRd
.multilife_ext_joint_survival <- function(t, x, y, model, ...) {
  tpxy(
    x = x,
    y = y,
    t = t,
    model = model,
    ...
  )
}


#' @noRd
.multilife_ext_single_survival <- function(t, x, model, ...) {
  tpx(
    t = t,
    x = x,
    model = model,
    ...
  )
}


#' @noRd
.multilife_ext_hazard <- function(age, model, ...) {
  hazard0(
    age,
    model = model,
    ...
  )
}


# -------------------------------------------------------------------------
# Contingent probabilities
# -------------------------------------------------------------------------

#' Contingent multi-life probabilities
#'
#' Computes probabilities associated with the order of death of two
#' independent lives over a specified term.
#'
#' \code{tqxy1()} computes the probability that the first life dies before
#' the second life within \eqn{n} years.
#'
#' \code{tqyx1()} computes the probability that the second life dies before
#' the first life within \eqn{n} years.
#'
#' \code{tqxy2()} computes the probability that the first life dies after
#' the second life but within \eqn{n} years.
#'
#' \code{tqyx2()} computes the probability that the second life dies after
#' the first life but within \eqn{n} years.
#'
#' These functions require a parametric survival model because an annual
#' life table does not uniquely determine the within-year order of death
#' without an additional interpolation assumption.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param n Term in years. May be scalar or vector.
#' @param tbl Optional life table object retained for backward compatibility.
#'   Continuous calculations currently require \code{model}.
#' @param model Parametric survival model.
#' @param ... Additional parameters passed to the survival and hazard
#'   functions.
#'
#' @return A numeric vector of contingent probabilities.
#'
#' @details
#' All calculations assume the two future lifetimes are independent.
#'
#' The probabilities are obtained by integrating the joint survival
#' function together with the appropriate force of mortality.
#'
#' The four functions partition the probability that one of the two lives
#' dies within the specified term according to the order of death.
#'
#' @examples
#' tqxy1(
#'   40, 50,
#'   n = 10,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' tqyx1(
#'   40, 50,
#'   n = 10,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' tqxy2(
#'   40, 50,
#'   n = 10,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' tqyx2(
#'   40, 50,
#'   n = 10,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' @name multilife_contingent_probabilities
#' @rdname multilife_contingent_probabilities
#' @export
tqxy1 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  n <- .multilife_ext_check_nonnegative(n, "n")

  values <- .multilife_ext_recycle(
    x,
    y,
    n,
    .names = c("x", "y", "n")
  )

  x <- values[[1]]
  y <- values[[2]]
  n <- values[[3]]

  vapply(
    seq_along(x),
    function(j) {
      upper <- .multilife_ext_upper(
        x = x[j],
        y = y[j],
        n = n[j],
        model = model,
        ...
      )

      .multilife_ext_integrate(
        function(t) {
          value <- .multilife_ext_joint_survival(
            t = t,
            x = x[j],
            y = y[j],
            model = model,
            ...
          ) *
            .multilife_ext_hazard(
              age = x[j] + t,
              model = model,
              ...
            )

          value[!is.finite(value)] <- 0
          value
        },
        upper = upper
      )
    },
    numeric(1)
  )
}

#' @rdname multilife_contingent_probabilities
#' @export
tqyx1 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  n <- .multilife_ext_check_nonnegative(n, "n")

  values <- .multilife_ext_recycle(
    x,
    y,
    n,
    .names = c("x", "y", "n")
  )

  x <- values[[1]]
  y <- values[[2]]
  n <- values[[3]]

  vapply(
    seq_along(x),
    function(j) {
      upper <- .multilife_ext_upper(
        x = x[j],
        y = y[j],
        n = n[j],
        model = model,
        ...
      )

      .multilife_ext_integrate(
        function(t) {
          value <- .multilife_ext_joint_survival(
            t = t,
            x = x[j],
            y = y[j],
            model = model,
            ...
          ) *
            .multilife_ext_hazard(
              age = y[j] + t,
              model = model,
              ...
            )

          value[!is.finite(value)] <- 0
          value
        },
        upper = upper
      )
    },
    numeric(1)
  )
}

#' @rdname multilife_contingent_probabilities
#' @export
tqxy2 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  n <- .multilife_ext_check_nonnegative(n, "n")

  values <- .multilife_ext_recycle(
    x,
    y,
    n,
    .names = c("x", "y", "n")
  )

  x <- values[[1]]
  y <- values[[2]]
  n <- values[[3]]

  qx <- 1 - .multilife_ext_single_survival(
    t = n,
    x = x,
    model = model,
    ...
  )

  out <- qx - tqxy1(
    x = x,
    y = y,
    n = n,
    model = model,
    ...
  )

  pmax(out, 0)
}

#' @rdname multilife_contingent_probabilities
#' @export
tqyx2 <- function(x, y, n, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  n <- .multilife_ext_check_nonnegative(n, "n")

  values <- .multilife_ext_recycle(
    x,
    y,
    n,
    .names = c("x", "y", "n")
  )

  x <- values[[1]]
  y <- values[[2]]
  n <- values[[3]]

  qy <- 1 - .multilife_ext_single_survival(
    t = n,
    x = y,
    model = model,
    ...
  )

  out <- qy - tqyx1(
    x = x,
    y = y,
    n = n,
    model = model,
    ...
  )

  pmax(out, 0)
}


# -------------------------------------------------------------------------
# Continuous multi-life annuities
# -------------------------------------------------------------------------

#' Continuous multi-life annuities
#'
#' Computes continuous joint-life, last-survivor, and reversionary whole-life
#' annuities for two independent lives.
#'
#' \code{abarxy()} computes the joint-life annuity.
#'
#' \code{abarxybar()} computes the last-survivor annuity.
#'
#' \code{abarx_y()} computes the reversionary annuity payable to the second
#' life after the death of the first.
#'
#' \code{abary_x()} computes the reversionary annuity payable to the first
#' life after the death of the second.
#'
#' These functions require a parametric survival model.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object retained for backward compatibility.
#'   Continuous calculations currently require \code{model}.
#' @param model Parametric survival model.
#' @param ... Additional parameters passed to the survival functions.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @details
#' Under independence, \code{abarxy()} represents the continuous joint-life
#' annuity, payable while both lives survive.
#'
#' The last-survivor annuity satisfies
#' \deqn{\bar{a}_{\overline{xy}} =
#' \bar{a}_x + \bar{a}_y - \bar{a}_{xy}.}
#'
#' The reversionary annuities are obtained as the difference between the
#' corresponding single-life annuity and the joint-life annuity.
#'
#' @examples
#' abarxy(
#'   x = 40,
#'   y = 50,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' abarxybar(
#'   x = 40,
#'   y = 50,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' @name continuous_multilife_annuities
#' @rdname continuous_multilife_annuities
#' @export
abarxy <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  delta <- log1p(i)

  vapply(
    seq_along(x),
    function(j) {
      upper <- .multilife_ext_upper(
        x = x[j],
        y = y[j],
        model = model,
        ...
      )

      .multilife_ext_integrate(
        function(t) {
          value <- exp(-delta[j] * t) *
            .multilife_ext_joint_survival(
              t = t,
              x = x[j],
              y = y[j],
              model = model,
              ...
            )

          value[!is.finite(value)] <- 0
          value
        },
        upper = upper
      )
    },
    numeric(1)
  )
}


#' @rdname continuous_multilife_annuities
#' @export
abarxybar <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  abarx(
    x = x,
    i = i,
    model = model,
    ...
  ) +
    abarx(
      x = y,
      i = i,
      model = model,
      ...
    ) -
    abarxy(
      x = x,
      y = y,
      i = i,
      model = model,
      ...
    )
}


#' @rdname continuous_multilife_annuities
#' @export
abarx_y <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  abarx(
    x = y,
    i = i,
    model = model,
    ...
  ) -
    abarxy(
      x = x,
      y = y,
      i = i,
      model = model,
      ...
    )
}


#' @rdname continuous_multilife_annuities
#' @export
abary_x <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  abarx(
    x = x,
    i = i,
    model = model,
    ...
  ) -
    abarxy(
      x = x,
      y = y,
      i = i,
      model = model,
      ...
    )
}

# -------------------------------------------------------------------------
# Continuous multi-life insurance
# -------------------------------------------------------------------------

#' Continuous multi-life insurance
#'
#' Computes continuous joint-life, last-survivor, and contingent whole-life
#' insurance values for two independent lives.
#'
#' \code{Abarxy()} computes joint-life insurance payable at the first death.
#'
#' \code{Abarxybar()} computes last-survivor insurance payable at the second
#' death.
#'
#' \code{Abarxy1()} and \code{Abaryx1()} compute contingent insurance payable
#' when the specified life dies first.
#'
#' \code{Abarxy2()} and \code{Abaryx2()} compute contingent insurance payable
#' when the specified life dies second.
#'
#' These functions require a parametric survival model.
#'
#' @param x Age of the first life. May be scalar or vector.
#' @param y Age of the second life. May be scalar or vector.
#' @param i Effective annual interest rate. May be scalar or vector.
#' @param tbl Optional life table object retained for backward compatibility.
#'   Continuous calculations currently require \code{model}.
#' @param model Parametric survival model.
#' @param ... Additional parameters passed to the survival and hazard
#'   functions.
#'
#' @return A numeric vector of actuarial present values.
#'
#' @details
#' Under independence, \code{Abarxy()} represents insurance payable at the
#' first death and \code{Abarxybar()} represents insurance payable at the
#' second death.
#'
#' The contingent functions distinguish both the life whose death triggers
#' payment and whether that death occurs first or second.
#'
#' @examples
#' Abarxy(
#'   x = 40,
#'   y = 50,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' Abarxy1(
#'   x = 40,
#'   y = 50,
#'   i = 0.05,
#'   model = "uniform",
#'   omega = 100
#' )
#'
#' @name continuous_multilife_insurance
#' @rdname continuous_multilife_insurance
#' @export
Abarxy <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  1 - log1p(i) * abarxy(
    x = x,
    y = y,
    i = i,
    model = model,
    ...
  )
}


#' @rdname continuous_multilife_insurance
#' @export
Abarxybar <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  Abarx(
    x = x,
    i = i,
    model = model,
    ...
  ) +
    Abarx(
      x = y,
      i = i,
      model = model,
      ...
    ) -
    Abarxy(
      x = x,
      y = y,
      i = i,
      model = model,
      ...
    )
}


#' @rdname continuous_multilife_insurance
#' @export
Abarxy1 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  delta <- log1p(i)

  vapply(
    seq_along(x),
    function(j) {
      upper <- .multilife_ext_upper(
        x = x[j],
        y = y[j],
        model = model,
        ...
      )

      .multilife_ext_integrate(
        function(t) {
          value <- exp(-delta[j] * t) *
            .multilife_ext_joint_survival(
              t = t,
              x = x[j],
              y = y[j],
              model = model,
              ...
            ) *
            .multilife_ext_hazard(
              age = x[j] + t,
              model = model,
              ...
            )

          value[!is.finite(value)] <- 0
          value
        },
        upper = upper
      )
    },
    numeric(1)
  )
}


#' @rdname continuous_multilife_insurance
#' @export
Abaryx1 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  delta <- log1p(i)

  vapply(
    seq_along(x),
    function(j) {
      upper <- .multilife_ext_upper(
        x = x[j],
        y = y[j],
        model = model,
        ...
      )

      .multilife_ext_integrate(
        function(t) {
          value <- exp(-delta[j] * t) *
            .multilife_ext_joint_survival(
              t = t,
              x = x[j],
              y = y[j],
              model = model,
              ...
            ) *
            .multilife_ext_hazard(
              age = y[j] + t,
              model = model,
              ...
            )

          value[!is.finite(value)] <- 0
          value
        },
        upper = upper
      )
    },
    numeric(1)
  )
}


#' @rdname continuous_multilife_insurance
#' @export
Abarxy2 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  out <- Abarx(
    x = x,
    i = i,
    model = model,
    ...
  ) -
    Abarxy1(
      x = x,
      y = y,
      i = i,
      model = model,
      ...
    )

  pmax(out, 0)
}


#' @rdname continuous_multilife_insurance
#' @export
Abaryx2 <- function(x, y, i, tbl = NULL, model = NULL, ...) {
  .multilife_ext_check_basis(
    tbl = tbl,
    model = model,
    require_continuous = TRUE
  )

  x <- .multilife_ext_check_nonnegative(x, "x")
  y <- .multilife_ext_check_nonnegative(y, "y")
  i <- .multilife_ext_check_interest(i)

  values <- .multilife_ext_recycle(
    x,
    y,
    i,
    .names = c("x", "y", "i")
  )

  x <- values[[1]]
  y <- values[[2]]
  i <- values[[3]]

  out <- Abarx(
    x = y,
    i = i,
    model = model,
    ...
  ) -
    Abaryx1(
      x = x,
      y = y,
      i = i,
      model = model,
      ...
    )

  pmax(out, 0)
}
