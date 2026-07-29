# Target contribution rate for a defined contribution plan

Calculates the contribution rate required to achieve a target
replacement ratio.

## Usage

``` r
contribution_rate_target(x, z, Sx, RR_target, i, adue_z, g = NULL, s = NULL)
```

## Arguments

- x:

  Scalar entry age.

- z:

  Scalar retirement age.

- Sx:

  Positive scalar salary at age `x`.

- RR_target:

  Scalar target replacement ratio in `[0, 1]`.

- i:

  Scalar annual effective investment return greater than `-1`.

- adue_z:

  Positive scalar whole-life annuity-due factor at retirement.

- g:

  Optional scalar annual salary growth rate greater than `-1`.

- s:

  Optional positive salary-scale vector of length `z - x`.

## Value

A numeric scalar. The result may exceed one when the target cannot be
achieved with a contribution rate no greater than 100 percent.

## Examples

``` r
contribution_rate_target(
  x = 30,
  z = 65,
  Sx = 60000,
  RR_target = 0.50,
  i = 0.06,
  adue_z = 11,
  g = 0.04
)
#> [1] 0.1052808
```
