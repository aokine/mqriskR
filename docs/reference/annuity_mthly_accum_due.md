# Temporary m-thly annuity-due actuarial accumulated value

Computes \$\$\ddot{s}\_{x:\overline{n}\|}^{(m)} =
\ddot{a}\_{x:\overline{n}\|}^{(m)} / {}\_nE_x\$\$

## Usage

``` r
sdotxn_m(x, n, m, i, model, ...)
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

Numeric vector of actuarial accumulated values.

## Examples

``` r
sdotxn_m(40, n = 10, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 14.32298
```
