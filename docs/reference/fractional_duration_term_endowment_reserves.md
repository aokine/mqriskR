# Fractional-duration term and endowment reserves

Computes fractional-duration reserves for term and endowment insurance
using linear interpolation between the reserve immediately after the
premium at duration \\t\\ and the reserve at duration \\t+1\\.

## Usage

``` r
tsVxn(x, n, t, s, i, tbl = NULL, model = NULL, ...)

tsVxn1(x, n, t, s, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Issue age. May be scalar or vector.

- n:

  Positive integer term. May be scalar or vector.

- t:

  Nonnegative integer duration satisfying \\t \< n\\. May be scalar or
  vector.

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

`tsVxn()` computes the fractional-duration reserve for endowment
insurance.

`tsVxn1()` computes the fractional-duration reserve for term insurance.

The fractional reserve is computed using

\$\$ {}\_{t+s}V = ({}\_tV + P)(1-s) + {}\_{t+1}V s. \$\$

The premium \\P\\ is the net annual premium for the corresponding
insurance contract (term or endowment).

## Examples

``` r
tsVxn(
  40,
  n = 20,
  t = 10,
  s = 0.5,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.388849

tsVxn1(
  40,
  n = 20,
  t = 10,
  s = 0.5,
  i = 0.05,
  model = "uniform",
  omega = 100
)
#> [1] 0.02775327
```
