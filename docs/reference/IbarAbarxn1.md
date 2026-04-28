# Fully continuous increasing n-year term insurance

Computes \$\$(\bar{I}\bar{A})\_{x:\overline{n}\|}^{1} = \int_0^n t\\
v^t\\ {}\_tp_x\\ \mu\_{x+t}\\ dt.\$\$

## Usage

``` r
IbarAbarxn1(x, n, i, model, ...)
```

## Arguments

- x:

  Age.

- n:

  Term.

- i:

  Effective annual interest rate.

- model:

  Parametric survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.
