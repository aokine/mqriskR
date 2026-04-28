# Total decrements \\d_x^{(\tau)}\\

Total decrements \\d_x^{(\tau)}\\

## Usage

``` r
dxtau(lxtau, qxj)
```

## Arguments

- lxtau:

  Number alive at age x in the multiple-decrement table.

- qxj:

  Numeric vector of cause-specific decrement probabilities.

## Value

Numeric scalar.

## Examples

``` r
dxtau(1000, c(0.011, 0.100))
#> [1] 111
```
