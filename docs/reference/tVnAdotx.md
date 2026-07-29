# Reserve for a deferred annuity-due

Computes the reserve during the deferral period, where \\0 \le t \< n\\.

## Usage

``` r
tVnAdotx(x, n, t, i, model = NULL, ..., tbl = NULL)
```

## Arguments

- x:

  Issue age. May be scalar or vector.

- n:

  Positive integer deferral period. May be scalar or vector.

- t:

  Nonnegative integer duration. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the actuarial functions.

- tbl:

  Optional life table object. Supply by name.

## Value

Numeric vector of reserve values.

## Examples

``` r
tVnAdotx(
  40,
  n = 20,
  t = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 3.915575
```
