# Accrued benefit under a final-average-salary plan

Computes the accrued benefit using salary and service history through
the valuation date.

## Usage

``` r
AB_fas(salary_history, p, fas_years = 3)
```

## Arguments

- salary_history:

  Positive numeric vector of annual salaries.

- p:

  Nonnegative scalar accrual percentage.

- fas_years:

  Positive integer number of years in the salary average.

## Value

A numeric scalar.

## Examples

``` r
AB_fas(
  salary_history = c(150000, 156000),
  p = 1,
  fas_years = 2
)
#> [1] 3060
```
