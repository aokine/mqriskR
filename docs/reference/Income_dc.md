# Retirement income from a defined contribution accumulation

Converts an accumulated account value into annual annuity-due income.
Scalar arguments are recycled to a common length.

## Usage

``` r
Income_dc(AVz, adue_z)
```

## Arguments

- AVz:

  Nonnegative accumulated value at retirement.

- adue_z:

  Positive whole-life annuity-due factor at retirement.

## Value

A numeric vector.

## Examples

``` r
Income_dc(AVz = 824211.35, adue_z = 12)
#> [1] 68684.28
```
