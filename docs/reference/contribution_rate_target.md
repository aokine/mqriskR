# Target contribution rate for a defined contribution plan

Solves Equation (18.5) for the contribution rate required to achieve a
target replacement ratio.

## Usage

``` r
contribution_rate_target(x, z, Sx, RR_target, i, adue_z, g = NULL, s = NULL)
```

## Arguments

- x:

  Entry age.

- z:

  Retirement age.

- Sx:

  Salary at age `x`.

- RR_target:

  Target replacement ratio.

- i:

  Annual effective interest rate.

- adue_z:

  Whole life annuity-due factor at age `z`.

- g:

  Optional constant annual salary growth rate.

- s:

  Optional salary scale vector of length `z - x`.

## Value

Required contribution rate.

## Examples

``` r
contribution_rate_target(
  x = 30, z = 65, Sx = 60000, RR_target = 0.50,
  i = 0.06, adue_z = 11, g = 0.04
)
#> [1] 0.1052808
```
