# Increasing n-year term insurance

Computes \$\$(IA)\_{x:\overline{n}\|}^{1} = \sum\_{k=0}^{n-1} (k+1)
v^{k+1} \Pr(K_x = k).\$\$

## Usage

``` r
IAxn1(x, n, i, tbl = NULL, model = NULL, ...)
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

  Optional parametric survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.
