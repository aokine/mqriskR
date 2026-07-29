# Actuarial present value of gross premiums

Computes the actuarial present value of premiums weighted by contract
persistency at the start of each policy year.

## Usage

``` r
APV_gross_premiums(G, r, p_tau)
```

## Arguments

- G:

  Nonnegative gross premium vector.

- r:

  Annual effective risk discount rate. May be scalar or vector; values
  must be greater than `-1`.

- p_tau:

  One-year in-force probabilities. For `n > 1`, this may have length
  `n - 1` or `n`; the final value is ignored when length `n`. For one
  premium, use `numeric(0)` or one value.

## Value

A numeric vector with one value for each rate in `r`.

## Examples

``` r
APV_gross_premiums(
  G = rep(95, 3),
  r = 0.10,
  p_tau = c(0.99858, 0.99847, 0.99834)
)
#> [1] 259.522
```
