# Net present value of a profit signature

Computes the Chapter 17 net present value: \$\$ NPV =
\sum\_{t=0}^{n}\frac{\Pi_t}{(1+r)^t}. \$\$

## Usage

``` r
NPV_profit(Pi, r)
```

## Arguments

- Pi:

  Profit signature vector.

- r:

  Risk discount rate.

## Value

Numeric scalar.

## Examples

``` r
Pi <- c(-15.00, 8.42, 8.39, 8.58)
NPV_profit(Pi, r = 0.10)
#> [1] 6.034711
```
