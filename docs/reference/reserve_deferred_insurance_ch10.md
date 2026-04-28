# Deferred insurance reserve functions

Chapter 10 reserve functions for deferred insurance contracts.

## Usage

``` r
tVnAx(x, n, t, i, model, ...)

htVnAx(x, n, h, t, i, model, ...)
```

## Arguments

- x:

  Issue age.

- n:

  Deferral period.

- t:

  Duration.

- i:

  Effective annual interest rate.

- model:

  Survival model.

- ...:

  Additional model parameters.

- h:

  Premium-paying period.

## Value

Numeric vector.

Numeric vector.

## Details

`tVnAx()` computes the reserve for an \\n\\-year deferred insurance
funded over the deferred period.

`htVnAx()` computes the reserve when premiums are limited to the first
\\h\\ years, with \\h \le n\\.

## Examples

``` r
tVnAx(40, n = 20, t = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.1548738
htVnAx(40, n = 20, h = 10, t = 5, i = 0.05, model = "uniform", omega = 100)
#> [1] 0.0874482
```
