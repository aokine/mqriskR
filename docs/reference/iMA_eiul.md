# Monthly-average index growth rate

Computes the monthly-average index growth rate from an initial index
value and 12 monthly closing values, matching Equation (16.13).

## Usage

``` r
iMA_eiul(index)
```

## Arguments

- index:

  Numeric vector of length 13 containing the initial index value
  followed by the 12 monthly closing values.

## Value

Numeric scalar.

## Examples

``` r
idx <- c(1000, 1020, 1100, 1150, 1080, 1040, 960, 1030, 1000, 1070, 1150, 1200, 1150)
iMA_eiul(idx)
#> [1] 0.07916667
```
