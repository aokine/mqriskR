# Projected annual benefit under a final-average-salary plan

Projects the annual benefit using a final-average-salary formula.

## Usage

``` r
PAB_fas(x, z, CASx, p, fas_years = 3, past_service = 0, g = NULL, s = NULL)
```

## Arguments

- x:

  Scalar current or entry age.

- z:

  Scalar retirement age.

- CASx:

  Positive scalar current annual salary.

- p:

  Nonnegative scalar accrual percentage, such as `2` for 2 percent.

- fas_years:

  Positive integer number of years in the final salary average.

- past_service:

  Nonnegative scalar years of service already completed.

- g:

  Optional scalar annual salary growth rate greater than `-1`.

- s:

  Optional positive salary-scale vector of length `z - x`.

## Value

A numeric scalar.

## Examples

``` r
PAB_fas(
  x = 35,
  z = 65,
  CASx = 60000,
  p = 2,
  fas_years = 3,
  g = 0.04
)
#> [1] 108008.7
```
