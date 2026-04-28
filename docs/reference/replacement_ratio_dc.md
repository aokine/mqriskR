# Replacement ratio for a defined contribution plan

Computes the replacement ratio defined in Equation (18.4).

## Usage

``` r
replacement_ratio_dc(x, z, Sx, c, i, adue_z, g = NULL, s = NULL)
```

## Arguments

- x:

  Entry age.

- z:

  Retirement age.

- Sx:

  Salary at age `x`.

- c:

  Contribution rate.

- i:

  Annual effective interest rate.

- adue_z:

  Whole life annuity-due factor at age `z`.

- g:

  Optional constant annual salary growth rate.

- s:

  Optional salary scale vector of length `z - x`.

## Value

Replacement ratio.

## Examples

``` r
replacement_ratio_dc(
  x = 30, z = 65, Sx = 50000, c = 0.10, i = 0.05,
  adue_z = 12, g = 0.04
)
#> [1] 0.3620377
```
