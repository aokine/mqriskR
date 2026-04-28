# Projected annual benefit under a career average earnings DB plan

Computes the projected annual benefit for a CAE plan using Equation
(18.9).

## Usage

``` r
PAB_cae(x, z, CASx, p, past_salary_total = 0, g = NULL, s = NULL)
```

## Arguments

- x:

  Current or entry age.

- z:

  Retirement age.

- CASx:

  Current annual salary at age `x`.

- p:

  Accrual percentage, e.g. `1` for 1 percent.

- past_salary_total:

  Optional total of actual past salaries.

- g:

  Optional constant annual salary growth rate.

- s:

  Optional salary scale vector of length `z - x`.

## Value

Projected annual benefit.

## Examples

``` r
PAB_cae(x = 30, z = 65, CASx = 100000, p = 1, g = 0.04)
#> [1] 73652.22
```
