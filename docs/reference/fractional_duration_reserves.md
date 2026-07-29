# Fractional-duration whole life reserves

Computes fractional-duration whole life reserves using linear
interpolation between the reserve immediately after the premium at
duration \\t\\ and the reserve at duration \\t+1\\.

## Usage

``` r
tsVx(x, t, s, i, tbl = NULL, model = NULL, ...)

meanVx(x, t, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age. May be scalar or vector.

- t:

  Nonnegative integer duration. May be scalar or vector.

- s:

  Fractional duration in \\\[0,1\]\\. May be scalar or vector.

- i:

  Effective annual interest rate. May be scalar or vector.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model.

- ...:

  Additional parameters passed to the actuarial functions.

## Value

A numeric vector of reserve values.

## Details

`tsVx()` computes the reserve at fractional duration \\t+s\\, where \\0
\le s \le 1\\.

`meanVx()` computes the reserve at the midpoint of the policy year
(\\s=0.5\\).

The fractional reserve is computed using

\$\$ {}\_{t+s}V_x = ({}\_tV_x + P_x)(1-s) + {}\_{t+1}V_x s. \$\$

The function `meanVx()` is a convenience wrapper corresponding to
\\s=0.5\\.

## Examples

``` r
tsVx(
  40,
  t = 10,
  s = 0.5,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.08762133

meanVx(
  40,
  t = 10,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.08762133
```
