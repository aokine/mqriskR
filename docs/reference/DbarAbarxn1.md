# Fully continuous decreasing n-year term insurance

Computes \$\$(\bar{D}\bar{A})\_{x:\overline{n}\|}^{1} = \int_0^n (n-t)\\
v^t\\ {}\_tp_x\\ \mu\_{x+t}\\ dt.\$\$

## Usage

``` r
DbarAbarxn1(x, n, i, model, ...)
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
