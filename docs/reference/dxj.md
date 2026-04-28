# Cause-specific decrements \\d_x^{(j)}\\

Cause-specific decrements \\d_x^{(j)}\\

## Usage

``` r
dxj(lxtau, qxj)
```

## Arguments

- lxtau:

  Number alive at age x in the multiple-decrement table.

- qxj:

  Numeric vector of cause-specific decrement probabilities.

## Value

Numeric vector.

## Examples

``` r
dxj(1000, c(0.011, 0.100))
#> [1]  11 100
```
