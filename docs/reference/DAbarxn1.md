# Piecewise-continuous decreasing n-year term insurance

Computes \$\$(D\bar{A})\_{x:\overline{n}\|}^{1} = \int_0^n \lfloor n+1-t
\rfloor\\ v^t\\ {}\_tp_x\\ \mu\_{x+t}\\ dt.\$\$

## Usage

``` r
DAbarxn1(x, n, i, model, ...)
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
