# Monthly-average index growth rate

Computes the monthly-average growth rate from an initial index value and
twelve monthly closing values.

## Usage

``` r
iMA_eiul(index)
```

## Arguments

- index:

  Numeric vector of length 13 containing a strictly positive initial
  index value followed by twelve nonnegative monthly closing values.

## Value

A numeric scalar.

## Examples

``` r
index <- c(
  1000, 1020, 1100, 1150, 1080, 1040, 960,
  1030, 1000, 1070, 1150, 1200, 1150
)
iMA_eiul(index)
#> [1] 0.07916667
```
