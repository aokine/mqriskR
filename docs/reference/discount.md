# Discount factor for compound interest

Computes the discount factor \\v^t = (1+i)^{-t}\\.

## Usage

``` r
discount(i, t)
```

## Arguments

- i:

  Effective interest rate. May be scalar or vector.

- t:

  Time. May be scalar or vector.

## Value

Numeric vector of discount factors.

## Examples

``` r
discount(0.05, 0:5)
#> [1] 1.0000000 0.9523810 0.9070295 0.8638376 0.8227025 0.7835262
discount(c(0.03, 0.05), 1)
#> [1] 0.9708738 0.9523810
```
