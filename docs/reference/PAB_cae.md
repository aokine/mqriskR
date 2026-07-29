# Projected annual benefit under a career-average-earnings plan

Projects the annual benefit using a career-average-earnings formula.

## Usage

``` r
PAB_cae(x, z, CASx, p, past_salary_total = 0, g = NULL, s = NULL)
```

## Arguments

- x:

  Scalar current or entry age.

- z:

  Scalar retirement age.

- CASx:

  Positive scalar current annual salary.

- p:

  Nonnegative scalar accrual percentage.

- past_salary_total:

  Nonnegative total of actual prior salaries.

- g:

  Optional scalar annual salary growth rate greater than `-1`.

- s:

  Optional positive salary-scale vector of length `z - x`.

## Value

A numeric scalar.

## Examples

``` r
PAB_cae(
  x = 30,
  z = 65,
  CASx = 100000,
  p = 1,
  g = 0.04
)
#> [1] 73652.22
```
