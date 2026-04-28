# m-thly term insurance APV

Computes \\A\_{x:\overline{n}\|}^{1(m)} = \sum\_{j=0}^{mn-1}
v^{(j+1)/m}\Pr(j/m \< T_x \le (j+1)/m)\\.

## Usage

``` r
Axn1_m(x, n, i, m, model, ...)
```

## Arguments

- x:

  Age.

- n:

  Term.

- i:

  Effective annual interest rate.

- m:

  Positive integer payment frequency.

- model:

  Parametric survival model name.

- ...:

  Additional model parameters passed to survival-model functions.

## Value

Numeric vector of APVs.
