# Continuous endowment insurance APV

Computes \\\bar{A}\_{x:\overline{n}\|} =
\bar{A}\_{x:\overline{n}\|}^{1} + v^n\\{}\_np_x\\.

## Usage

``` r
Abarxn(x, n, i, model, ...)
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

Numeric vector of APVs.
