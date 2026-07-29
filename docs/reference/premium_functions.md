# Premium, loss, and expense functions

Functions for net and gross premiums, present-value-of-loss moments,
continuous-payment premium rates, and premiums payable more frequently
than annually.

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

  Age. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to survival-model functions.

- n:

  Term. May be scalar or vector of nonnegative integers.

- t:

  Premium-paying period. May be scalar or vector of nonnegative
  integers.

- tol:

  Numerical tolerance for functions that truncate infinite sums.

- m:

  Number of payments per year. Must be a positive integer scalar.

- P:

  Premium amount or premium rate. May be scalar or vector.

- k_max:

  Maximum summation horizon for functions that truncate infinite sums.

- benefit:

  Benefit amount. May be scalar or vector.

- first_premium_pct:

  First-year premium expense proportion. May be scalar or vector.

- renewal_premium_pct:

  Renewal premium expense proportion. May be scalar or vector.

- first_policy_exp:

  First-year fixed expense. May be scalar or vector.

- renewal_policy_exp:

  Renewal fixed expense after the first year. May be scalar or vector.

- settlement_exp:

  Settlement expense incurred at benefit payment. May be scalar or
  vector.

## Value

Numeric vector.

## Details

The functions include:

- net annual premiums under the equivalence principle,

- limited-payment premiums,

- continuous-payment premium rates,

- fully continuous premium rates,

- true fractional premiums,

- present-value-of-loss means and variances,

- a basic gross premium formula for whole life insurance.

The discrete premium functions may be evaluated from either a life table
supplied through `tbl` or a parametric survival model supplied through
`model`. Continuous-payment premium functions are evaluated through the
corresponding continuous insurance and annuity functions.

Scalar inputs retain their existing behavior. Where mathematically
meaningful, numeric inputs may also be vectors. Inputs are evaluated
elementwise when they have a common length, and scalar inputs are
recycled.
