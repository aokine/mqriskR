# Retrospective term insurance reserve

Computes the retrospective reserve for a term insurance at a duration
satisfying \\0 \le t \le n\\. At expiry, the reserve is 0.

## Usage

``` r
tVxn1_ret(x, n, t, i, model = NULL, ..., tbl = NULL)
```

## Arguments

- x:

  Issue age. May be scalar or vector.

- n:

  Nonnegative integer contract term. May be scalar or vector.

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

Numeric vector of retrospective reserve values.

## Examples

``` r
tVxn1_ret(
  40,
  n = 20,
  t = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.01836757
```
