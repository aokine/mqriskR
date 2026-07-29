# m-thly contingent annuity functions

Life annuities payable m-thly.

## Usage

``` r
ax_m(x, m, i, model, ..., k_max = 2e+05, tol = 1e-12)

adotx_m(x, m, i, model, ..., k_max = 2e+05, tol = 1e-12)

axn_m(x, n, m, i, model, ...)

adotxn_m(x, n, m, i, model, ...)

nax_m(x, n, m, i, model, ..., k_max = 2e+05, tol = 1e-12)

nadotx_m(x, n, m, i, model, ..., k_max = 2e+05, tol = 1e-12)

sxn_m(x, n, m, i, model, ...)

sdotxn_m(x, n, m, i, model, ...)
```

## Arguments

- x:

  Age.

- m:

  Number of payments per year.

- i:

  Effective annual interest rate.

- model:

  Survival model name.

- ...:

  Additional parameters passed to the survival model.

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

- n:

  Term or deferral period in years.

## Value

Numeric vector containing the requested m-thly annuity value or
actuarial accumulated value.
