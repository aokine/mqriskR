# Point-to-point index growth rates

Computes annual point-to-point index growth rates from successive index
values, matching Equation (16.12).

## Usage

``` r
iP_eiul(index)
```

## Arguments

- index:

  Numeric vector of index closing values.

## Value

Numeric vector of growth rates.

## Examples

``` r
iP_eiul(c(1000, 1050, 1200, 1100))
#> [1]  0.05000000  0.14285714 -0.08333333
```
