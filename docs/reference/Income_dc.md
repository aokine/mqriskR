# Retirement income from a defined contribution accumulation

Converts a defined contribution accumulation to annual annuity-due
income.

## Usage

``` r
Income_dc(AVz, adue_z)
```

## Arguments

- AVz:

  Accumulated value at retirement.

- adue_z:

  Whole life annuity-due factor at retirement age.

## Value

Annual retirement income.

## Examples

``` r
Income_dc(AVz = 824211.35, adue_z = 12)
#> [1] 68684.28
```
