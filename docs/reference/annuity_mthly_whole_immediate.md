# Whole life m-thly annuity-immediate

Computes the exact whole life m-thly annuity-immediate, with annual
payment rate 1 split into m equal payments of size 1/m.

## Usage

``` r
ax_m(x, m, i, model, ..., k_max = 2e+05, tol = 1e-12)
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

\$\$a_x^{(m)} = \frac{1}{m}\sum\_{t=1}^{\infty}
v^{t/m}\\{}\_{t/m}p_x\$\$

## Examples

``` r
ax_m(40, m = 12, i = 0.05, model = "uniform", omega = 100)
#> [1] 13.82774
```
