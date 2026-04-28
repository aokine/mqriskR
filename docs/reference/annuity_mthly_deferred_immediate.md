# Deferred whole life m-thly annuity-immediate

Computes the deferred whole life m-thly annuity-immediate using
\$\${}\_nE_x \\ a\_{x+n}^{(m)}\$\$

## Usage

``` r
nax_m(x, n, m, i, model, ..., k_max = 2e+05, tol = 1e-12)
```

## Arguments

- x:

  Age.

- n:

  Deferral period in years.

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

## Value

Numeric vector of annuity values.

## Examples

``` r
nax_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 6.540886
```
