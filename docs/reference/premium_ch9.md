# Chapter 9 premium, loss, and expense functions

Functions in this file implement Chapter 9 funding-plan formulas,
including:

- net annual premiums under the equivalence principle,

- limited-payment premiums,

- continuous-payment premium rates,

- fully continuous premium rates,

- true fractional premiums,

- present-value-of-loss means and variances,

- a basic gross premium formula for whole life insurance.

## Usage

``` r
Px(x, i, tbl = NULL, model = NULL, ...)

Pxn1(x, n, i, tbl = NULL, model = NULL, ...)

PnEx(x, n, i, tbl = NULL, model = NULL, ...)

Pxn(x, n, i, tbl = NULL, model = NULL, ...)

tPx(x, t, i, tbl = NULL, model = NULL, ...)

tPxn1(x, n, t, i, tbl = NULL, model = NULL, ...)

tPnEx(x, n, t, i, tbl = NULL, model = NULL, ...)

tPxn(x, n, t, i, tbl = NULL, model = NULL, ...)

PnAx(x, n, i, tbl = NULL, model = NULL, ...)

tPnAx(x, n, t, i, tbl = NULL, model = NULL, ...)

Pbarx(x, i, model, ..., tol = 1e-10)

Pbarxn1(x, n, i, model, ...)

Pbarxn(x, n, i, model, ...)

PbarAbarx(x, i, model, ..., tol = 1e-10)

PbarAbarxn1(x, n, i, model, ...)

PbarAbarxn(x, n, i, model, ...)

Px_m(x, m, i, tbl = NULL, model = NULL, ...)

Pxn1_m(x, n, m, i, tbl = NULL, model = NULL, ...)

Pxn_m(x, n, m, i, tbl = NULL, model = NULL, ...)

PnAx_m(x, n, m, i, tbl = NULL, model = NULL, ...)

EL0x(x, P, i, tbl = NULL, model = NULL, ...)

varL0x(x, P, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000)

EL0xn1(x, n, P, i, tbl = NULL, model = NULL, ...)

varL0xn1(x, n, P, i, tbl = NULL, model = NULL, ...)

EL0xn(x, n, P, i, tbl = NULL, model = NULL, ...)

varL0xn(x, n, P, i, tbl = NULL, model = NULL, ...)

EL0barAbarx(x, P, i, model, ..., tol = 1e-10)

varL0barAbarx(x, P, i, model, ...)

Gx(
  x,
  i,
  benefit = 1,
  first_premium_pct = 0,
  renewal_premium_pct = 0,
  first_policy_exp = 0,
  renewal_policy_exp = 0,
  settlement_exp = 0,
  tbl = NULL,
  model = NULL,
  ...
)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to survival-model functions.

- n:

  Term.

- t:

  Premium-paying period.

- tol:

  Numerical tolerance for functions that truncate infinite sums.

- m:

  Number of payments per year.

- P:

  Premium amount or premium rate.

- k_max:

  Maximum summation horizon for functions that truncate infinite sums.

- benefit:

  Benefit amount.

- first_premium_pct:

  First-year premium expense proportion.

- renewal_premium_pct:

  Renewal premium expense proportion.

- first_policy_exp:

  First-year fixed expense.

- renewal_policy_exp:

  Renewal fixed expense each year after the first.

- settlement_exp:

  Settlement expense incurred at benefit payment.

## Value

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

Numeric vector.

## Details

Naming follows Chapter 9 notation as closely as possible:

- `Px()` = whole life annual premium

- `Pxn1()` = term insurance annual premium

- `PnEx()` = pure endowment annual premium

- `Pxn()` = endowment insurance annual premium

- `tPx()` = limited-payment whole life annual premium

- `tPxn1()` = limited-payment term insurance annual premium

- `tPnEx()` = limited-payment pure endowment annual premium

- `tPxn()` = limited-payment endowment insurance annual premium

- `PnAx()` = deferred insurance annual premium

- `tPnAx()` = limited-payment deferred insurance annual premium

- `Pbarx()` = continuous-payment premium for discrete whole life
  insurance

- `Pbarxn1()` = continuous-payment premium for discrete term insurance

- `Pbarxn()` = continuous-payment premium for discrete endowment
  insurance

- `PbarAbarx()` = fully continuous premium for continuous whole life
  insurance

- `PbarAbarxn1()` = fully continuous premium for continuous term
  insurance

- `PbarAbarxn()` = fully continuous premium for continuous endowment
  insurance

- `Px_m()` = true fractional whole life annual premium

- `Pxn1_m()` = true fractional term insurance annual premium

- `Pxn_m()` = true fractional endowment insurance annual premium

- `PnAx_m()` = true fractional deferred insurance annual premium

The discrete premium functions can be evaluated from either a life table
via `tbl = ...` or from a parametric model via `model = ...`.

The continuous-premium-rate functions use the continuous annuity
functions already in the package, so they are written for the parametric
survival model framework.
