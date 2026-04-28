# AG 38 prefunding ratio

Computes the prefunding ratio in Equation (16.18), capped at 1.

## Usage

``` r
ag38_prefunding_ratio(excess_payment, nsp_required)
```

## Arguments

- excess_payment:

  Excess payment or shadow-fund amount.

- nsp_required:

  Net single premium required to fully fund the guarantee.

## Value

Numeric scalar.

## Examples

``` r
ag38_prefunding_ratio(60000, 100000)
#> [1] 0.6
```
