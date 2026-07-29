# Profit margin

Computes net present value divided by the actuarial present value of
gross premiums. Scalar arguments are recycled to the common length.

## Usage

``` r
profit_margin(NPV, APV_GP)
```

## Arguments

- NPV:

  Net present value of profits.

- APV_GP:

  Positive actuarial present value of gross premiums.

## Value

A numeric vector.

## Examples

``` r
profit_margin(6.03, 259.52)
#> [1] 0.0232352
```
