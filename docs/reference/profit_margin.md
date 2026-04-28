# Profit margin

Computes the Chapter 17 profit margin: \$\$ \text{Profit Margin} =
\frac{NPV}{APV\_{GP}}. \$\$

## Usage

``` r
profit_margin(NPV, APV_GP)
```

## Arguments

- NPV:

  Net present value of profits.

- APV_GP:

  Actuarial present value of gross premiums.

## Value

Numeric scalar.

## Examples

``` r
profit_margin(6.03, 259.52)
#> [1] 0.0232352
```
