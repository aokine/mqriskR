# Salary scale under constant annual growth

Constructs salary-scale values under a constant annual growth rate.

## Usage

``` r
salary_scale(k, g, base_age = min(k), s_base = 1)
```

## Arguments

- k:

  Numeric vector of ages or durations.

- g:

  Annual salary growth rate greater than `-1`.

- base_age:

  Scalar age or duration at which the scale is normalized.

- s_base:

  Positive scalar salary-scale value at `base_age`.

## Value

A numeric vector with the same length as `k`.

## Examples

``` r
salary_scale(k = 30:34, g = 0.04, base_age = 30)
#> [1] 1.000000 1.040000 1.081600 1.124864 1.169859
```
