# Second moment of continuous term insurance PV

Computes \\{}^{2}\bar{A}\_{x:\overline{n}\|}^{1}\\ by evaluating
\\\bar{A}\_{x:\overline{n}\|}^{1}\\ at doubled force.

## Usage

``` r
A2barxn1(x, n, i, model, ...)
```

## Arguments

- x:

  Age.

- n:

  Term.

- i:

  Effective annual interest rate.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

## Value

Numeric vector of second moments.
