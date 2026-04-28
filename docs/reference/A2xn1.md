# Second moment of term insurance PV

Computes \\{}^{2}A\_{x:\overline{n}\|}^{1}\\ by evaluating
\\A\_{x:\overline{n}\|}^{1}\\ at doubled force.

## Usage

``` r
A2xn1(x, n, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age.

- n:

  Term.

- i:

  Effective annual interest rate.

- tbl:

  Optional life table object.

- model:

  Optional parametric survival model name.

- ...:

  Additional arguments passed to survival-model functions.

## Value

Numeric vector of second moments.
