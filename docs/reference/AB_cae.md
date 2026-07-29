# Accrued benefit under a career-average-earnings plan

Computes the accrued benefit using salary history through the valuation
date.

## Usage

``` r
AB_cae(salary_history, p)
```

## Arguments

- salary_history:

  Positive numeric vector of annual salaries.

- p:

  Nonnegative scalar accrual percentage.

## Value

A numeric scalar.

## Examples

``` r
AB_cae(
  salary_history = c(100000, 104000, 108160),
  p = 1
)
#> [1] 3121.6
```
