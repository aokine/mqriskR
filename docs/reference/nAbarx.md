# Continuous deferred insurance APV

Computes \\{}\_{n\mid}\bar{A}\_x = v^n\\{}\_np_x\\\bar{A}\_{x+n}\\.

## Usage

``` r
nAbarx(x, n, i, model, ...)
```

## Arguments

- x:

  Age.

- n:

  Deferral period.

- i:

  Effective annual interest rate.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

## Value

Numeric vector of APVs.
