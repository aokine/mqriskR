# Second moment of continuous deferred insurance PV

Computes \\{}^{2}{}\_{n\mid}\bar{A}\_x\\ by evaluating
\\{}\_{n\mid}\bar{A}\_x\\ at doubled force.

## Usage

``` r
A2nAbarx(x, n, i, model, ...)
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

Numeric vector of second moments.
