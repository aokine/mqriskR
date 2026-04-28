# Salary scale under constant annual growth

Constructs salary scale factors \\s_k\\ under a constant annual growth
rate.

## Usage

``` r
salary_scale(k, g, base_age = min(k), s_base = 1)
```

## Arguments

- k:

  Numeric vector of ages.

- g:

  Annual salary growth rate.

- base_age:

  Age at which the scale is normalized.

- s_base:

  Salary scale value at `base_age`.

## Value

Numeric vector of salary scale factors.

## Examples

``` r
salary_scale(k = 30:34, g = 0.04, base_age = 30)
#> [1] 1.000000 1.040000 1.081600 1.124864 1.169859
```
