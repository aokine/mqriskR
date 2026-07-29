# Accumulated value of defined contribution plan contributions

Computes the accumulated value at retirement from contributions paid at
the beginning of each year and accumulated to retirement.

## Usage

``` r
AVz_dc(x, z, Sx, c, i, g = NULL, s = NULL)
```

## Arguments

- x:

  Scalar entry age.

- z:

  Scalar retirement age.

- Sx:

  Positive scalar salary at age `x`.

- c:

  Scalar contribution rate in `[0, 1]`.

- i:

  Scalar annual effective investment return greater than `-1`.

- g:

  Optional scalar annual salary growth rate greater than `-1`.

- s:

  Optional positive salary-scale vector of length `z - x`.

## Value

A numeric scalar.

## Examples

``` r
AVz_dc(
  x = 30,
  z = 65,
  Sx = 50000,
  c = 0.10,
  i = 0.05,
  g = 0.04
)
#> [1] 824211.3
```
