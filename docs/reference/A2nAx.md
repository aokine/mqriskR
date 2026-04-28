# Second moment of deferred insurance PV

Computes \\{}^{2}{}\_{n\mid}A_x\\ by evaluating \\{}\_{n\mid}A_x\\ at
doubled force.

## Usage

``` r
A2nAx(x, n, i, tbl = NULL, model = NULL, ..., tol = 1e-12, k_max = 5000)
```

## Arguments

- x:

  Age.

- n:

  Deferral period.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to survival-model functions.

- tol:

  Numerical tolerance for truncating infinite sums.

- k_max:

  Maximum number of terms in the sum.

## Value

Numeric vector of second moments.
