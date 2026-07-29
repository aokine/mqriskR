# AG 38 reserve calculation

Computes the prefunding ratio, net additional amount, reduced deficiency
reserve, intermediate reserve, and increased basic reserve.

## Usage

``` r
ag38_reserve_ul(
  basic_reserve,
  deficiency_reserve = 0,
  excess_payment,
  nsp_required,
  valuation_nsp,
  surrender_charge = 0
)
```

## Arguments

- basic_reserve:

  Nonnegative basic reserve.

- deficiency_reserve:

  Nonnegative deficiency reserve.

- excess_payment:

  Nonnegative excess payment or shadow-fund amount.

- nsp_required:

  Positive net single premium required to fully fund the guarantee.

- valuation_nsp:

  Nonnegative valuation net single premium.

- surrender_charge:

  Nonnegative surrender charge.

## Value

For scalar inputs, a named list. For vectorized inputs, a data frame
containing the same calculated quantities.

## Details

Scalar inputs preserve the original named-list output. Vectorized inputs
return a data frame with one row per calculation.

## Examples

``` r
ag38_reserve_ul(
  basic_reserve = 10000,
  deficiency_reserve = 0,
  excess_payment = 60000,
  nsp_required = 100000,
  valuation_nsp = 150000,
  surrender_charge = 5000
)
#> $prefunding_ratio
#> [1] 0.6
#> 
#> $net_amount_additional
#> [1] 84000
#> 
#> $reduced_deficiency_reserve
#> [1] 0
#> 
#> $step8_reserve
#> [1] 89000
#> 
#> $increased_basic_reserve
#> [1] 89000
#> 
#> $final_reserve
#> [1] 89000
#> 
```
