# Full preliminary term modified premiums and reserves

Computes modified premiums and reserves under the full preliminary term
method for whole-life insurance.

## Usage

``` r
alphaF(x, i, tbl = NULL, model = NULL, ...)

betaF(x, i, tbl = NULL, model = NULL, ...)

tVFx(x, t, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the actuarial functions.

- t:

  Nonnegative integer duration. May be scalar or vector.

## Value

A numeric vector of modified premiums or reserves.

## Details

`alphaF()` computes the first-year modified premium \\\alpha^F = vq_x\\.

`betaF()` computes the renewal modified premium \\\beta^F = P\_{x+1}\\.

`tVFx()` computes the full preliminary term reserve. The reserve is zero
at durations 0 and 1. For \\t \> 1\\, it equals the net level premium
reserve at duration \\t - 1\\ for a policy issued at age \\x + 1\\.

Under the full preliminary term method, the first policy year is treated
as one-year term insurance. The first-year modified premium is therefore
the actuarial present value of one-year term insurance at age \\x\\.

Renewal premiums are based on a whole-life policy issued one year later,
at age \\x + 1\\. Accordingly, reserves after the first policy year are
obtained from the corresponding net level premium reserve for that
deferred issue age.

## Examples

``` r
alphaF(
  x = 40,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.01587302

betaF(
  x = 40,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.02240155

tVFx(
  x = 40,
  t = 5,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.02773584
```
