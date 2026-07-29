# Account-value to guaranteed-fund ratio

Computes the ratio of account value to guaranteed maturity fund, capped
at one.

## Usage

``` r
rt_ul(AV, GMF)
```

## Arguments

- AV:

  Nonnegative account value.

- GMF:

  Positive guaranteed maturity fund.

## Value

A numeric vector.

## Examples

``` r
rt_ul(AV = 4918.20, GMF = 14678.57)
#> [1] 0.3350599
```
