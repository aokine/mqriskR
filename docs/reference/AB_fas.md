# Accrued benefit for a final average salary plan

Computes the accrued benefit at the current date using service and
salary history only.

## Usage

``` r
AB_fas(salary_history, p, fas_years = 3)
```

## Arguments

- salary_history:

  Numeric vector of annual salaries to date.

- p:

  Accrual percentage, e.g. `2` for 2 percent.

- fas_years:

  Number of years in the final average salary average.

## Value

Accrued benefit.

## Examples

``` r
AB_fas(salary_history = c(150000, 156000), p = 1, fas_years = 2)
#> [1] 3060
```
