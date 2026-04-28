# AG 38 reserve calculation

Computes the main quantities in the Chapter 16 AG 38 reserve
calculation, including the prefunding ratio, reduced deficiency reserve,
Step (8) reserve, and final increased basic reserve.

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

  Basic reserve.

- deficiency_reserve:

  Deficiency reserve.

- excess_payment:

  Excess payment or shadow-fund amount.

- nsp_required:

  Net single premium required to fully fund the guarantee.

- valuation_nsp:

  Valuation net single premium.

- surrender_charge:

  Applicable surrender charge.

## Value

A named list.

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
