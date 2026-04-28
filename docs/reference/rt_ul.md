# Ratio \\r_t = AV_t / GMF_t\\ capped at 1

Computes the ratio in Equation (16.16).

## Usage

``` r
rt_ul(AV, GMF)
```

## Arguments

- AV:

  Account value.

- GMF:

  Guaranteed maturity fund.

## Value

Numeric scalar.

## Examples

``` r
rt_ul(AV = 4918.20, GMF = 14678.57)
#> [1] 0.3350599
```
