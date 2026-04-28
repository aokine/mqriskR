# Whole life m-thly annuity-due

Computes the exact whole life m-thly annuity-due.

## Usage

``` r
adotx_m(x, m, i, model, ..., k_max = 2e+05, tol = 1e-12)
```

## Arguments

- x:

  Age.

- m:

  Number of payments per year.

- i:

  Effective annual interest rate.

- model:

  Survival model name.

- ...:

  Additional parameters passed to the survival model.

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

## Value

Numeric vector of annuity values.

## Details

\$\$\ddot{a}\_x^{(m)} = \frac{1}{m}\sum\_{t=0}^{\infty}
v^{t/m}\\{}\_{t/m}p_x\$\$

## Examples

``` r
adotx_m(40, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 13.91107
```
