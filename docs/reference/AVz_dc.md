# Accumulated value of defined contribution plan contributions

Computes the accumulated value at retirement age \\z\\ for the defined
contribution model in Equation (18.1).

## Usage

``` r
AVz_dc(x, z, Sx, c, i, g = NULL, s = NULL)
```

## Arguments

- x:

  Entry age.

- z:

  Retirement age.

- Sx:

  Salary at age `x`.

- c:

  Contribution rate as a proportion of salary.

- i:

  Annual effective interest rate.

- g:

  Optional constant annual salary growth rate.

- s:

  Optional salary scale vector of length `z - x`.

## Value

The accumulated value of contributions at age `z`.

## Examples

``` r
AVz_dc(x = 30, z = 65, Sx = 50000, c = 0.10, i = 0.05, g = 0.04)
#> [1] 824211.3
```
