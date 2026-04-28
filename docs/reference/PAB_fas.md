# Projected annual benefit under a final average salary DB plan

Computes the projected annual benefit for a final average salary plan
using Equation (18.6).

## Usage

``` r
PAB_fas(x, z, CASx, p, fas_years = 3, past_service = 0, g = NULL, s = NULL)
```

## Arguments

- x:

  Current or entry age.

- z:

  Retirement age.

- CASx:

  Current annual salary at age `x`.

- p:

  Accrual percentage, e.g. `2` for 2 percent.

- fas_years:

  Number of years in the final average salary period.

- past_service:

  Past years of service already completed at age `x`.

- g:

  Optional constant annual salary growth rate.

- s:

  Optional salary scale vector of length `z - x`.

## Value

Projected annual benefit.

## Examples

``` r
PAB_fas(x = 35, z = 65, CASx = 60000, p = 2, fas_years = 3, g = 0.04)
#> [1] 108008.7
```
