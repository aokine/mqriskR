# Varying-payment annuity functions

Increasing and decreasing life annuity functions.

## Usage

``` r
Iax(x, i, model = NULL, ..., tbl = NULL, k_max = 5000, tol = 1e-12)

Iaxn(x, n, i, model = NULL, ..., tbl = NULL)

Daxn(x, n, i, model = NULL, ..., tbl = NULL)

Iadotx(x, i, model = NULL, ..., tbl = NULL, k_max = 5000, tol = 1e-12)

Iadotxn(x, n, i, model = NULL, ..., tbl = NULL)

Dadotxn(x, n, i, model = NULL, ..., tbl = NULL)

Iabarx(x, i, model, ..., tol = 1e-10)

Iabarxn(x, n, i, model, ...)

Dabarxn(x, n, i, model, ...)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- model:

  Optional survival model name.

- ...:

  Additional model parameters.

- tbl:

  Optional life table object for discrete functions.

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

- n:

  Term in years.

## Value

Numeric vector containing the requested increasing or decreasing annuity
value.
