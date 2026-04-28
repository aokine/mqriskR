# Accrued benefit for a career average earnings plan

Computes the accrued benefit using actual salary history only.

## Usage

``` r
AB_cae(salary_history, p)
```

## Arguments

- salary_history:

  Numeric vector of annual salaries to date.

- p:

  Accrual percentage.

## Value

Accrued benefit.

## Examples

``` r
AB_cae(salary_history = c(100000, 104000, 108160), p = 1)
#> [1] 3121.6
```
