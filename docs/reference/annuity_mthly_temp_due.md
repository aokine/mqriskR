# Temporary m-thly annuity-due

Computes the exact temporary m-thly annuity-due.

## Usage

``` r
adotxn_m(x, n, m, i, model, ...)
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

\$\$\ddot{a}\_{x:\overline{n}\|}^{(m)} = \frac{1}{m}\sum\_{t=0}^{mn-1}
v^{t/m}\\{}\_{t/m}p_x\$\$

## Examples

``` r
adotxn_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 7.327554
```
