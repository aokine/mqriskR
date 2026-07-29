# AG 38 prefunding ratio

Computes the excess-payment-to-required-net-single-premium ratio, capped
at one.

## Usage

``` r
ag38_prefunding_ratio(excess_payment, nsp_required)
```

## Arguments

- excess_payment:

  Nonnegative excess payment or shadow-fund amount.

- nsp_required:

  Positive net single premium required to fully fund the guarantee.

## Value

A numeric vector.

## Examples

``` r
ag38_prefunding_ratio(60000, 100000)
#> [1] 0.6
```
