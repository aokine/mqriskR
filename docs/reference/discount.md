# Discount factor for compound interest

Discount factor for compound interest

## Usage

``` r
discount(i, t)
```

## Arguments

- i:

  Effective interest rate.

- t:

  Time (can be vector).

## Value

Discount factor at time t.

## Examples

``` r
discount(0.05, 0:5)
#> [1] 1.0000000 0.9523810 0.9070295 0.8638376 0.8227025 0.7835262
```
