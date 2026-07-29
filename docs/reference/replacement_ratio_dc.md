# Replacement ratio for a defined contribution plan

Computes annual retirement income divided by salary in the final
pre-retirement year.

## Usage

``` r
replacement_ratio_dc(x, z, Sx, c, i, adue_z, g = NULL, s = NULL)
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

- adue_z:

  Positive scalar whole-life annuity-due factor at retirement.

- g:

  Optional scalar annual salary growth rate greater than `-1`.

- s:

  Optional positive salary-scale vector of length `z - x`.

## Value

A numeric scalar.

## Examples

``` r
replacement_ratio_dc(
  x = 30,
  z = 65,
  Sx = 50000,
  c = 0.10,
  i = 0.05,
  adue_z = 12,
  g = 0.04
)
#> [1] 0.3620377
```
