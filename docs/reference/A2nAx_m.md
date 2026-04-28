# Second moment of m-thly deferred insurance PV

Second moment of m-thly deferred insurance PV

## Usage

``` r
A2nAx_m(x, n, i, m, model, ..., tol = 1e-12, j_max = 100000L)
```

## Arguments

- x:

  Age.

- n:

  Deferral period.

- i:

  Effective annual interest rate.

- m:

  Positive integer payment frequency.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

- tol:

  Numerical tolerance for truncating the infinite sum.

- j_max:

  Maximum number of m-thly intervals in the sum.

## Value

Numeric vector of second moments.
