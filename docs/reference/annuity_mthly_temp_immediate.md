# Temporary m-thly annuity-immediate

Computes the exact temporary m-thly annuity-immediate.

## Usage

``` r
axn_m(x, n, m, i, model, ...)
```

## Arguments

- x:

  Age.

- n:

  Term in years.

- m:

  Number of payments per year.

- i:

  Effective annual interest rate.

- model:

  Survival model name.

- ...:

  Additional parameters passed to the survival model.

## Value

Numeric vector of annuity values.

## Details

\$\$a\_{x:\overline{n}\|}^{(m)} = \frac{1}{m}\sum\_{t=1}^{mn}
v^{t/m}\\{}\_{t/m}p_x\$\$

## Examples

``` r
axn_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 7.286853
```
