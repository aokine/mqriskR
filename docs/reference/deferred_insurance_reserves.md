# Deferred insurance reserves

Computes prospective reserves for deferred whole life insurance
contracts.

## Usage

``` r
tVnAx(x, n, t, i, model = NULL, ..., tbl = NULL)

htVnAx(x, n, h, t, i, model = NULL, ..., tbl = NULL)
```

## Arguments

- x:

  Issue age.

- n:

  Deferral period in years.

- t:

  Duration in years.

- i:

  Effective annual interest rate.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to the underlying actuarial functions.

- tbl:

  Optional life table object.

- h:

  Premium-paying period in years.

## Value

A numeric vector of prospective reserve values.

## Details

`tVnAx()` computes the reserve at duration `t` for an `n`-year deferred
whole life insurance funded by level annual premiums during the deferral
period.

`htVnAx()` computes the corresponding reserve when premiums are limited
to the first `h` years, where `h <= n`.

## Examples

``` r
tVnAx(
  x = 40, n = 20, t = 10, i = 0.05,
  model = "uniform", omega = 100
)
#> [1] 0.1548738

htVnAx(
  x = 40, n = 20, h = 10, t = 5, i = 0.05,
  model = "uniform", omega = 100
)
#> [1] 0.0874482
```
