# Point-to-point index growth rates

Computes consecutive point-to-point growth rates from index values.

## Usage

``` r
iP_eiul(index)
```

## Arguments

- index:

  Numeric vector of strictly positive index values.

## Value

A numeric vector with length one less than `index`.

## Examples

``` r
iP_eiul(c(1000, 1050, 1200, 1100))
#> [1]  0.05000000  0.14285714 -0.08333333
```
