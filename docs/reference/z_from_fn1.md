# Spot rates from forward one-year rates

Converts annual effective forward one-year rates \\f\_{0,1}, f\_{1,1},
\dots, f\_{n-1,1}\\ into annual effective spot rates \\z_1, z_2, \dots,
z_n\\.

## Usage

``` r
z_from_fn1(fn1)
```

## Arguments

- fn1:

  Numeric vector of annual effective forward one-year rates.

## Value

Numeric vector of annual effective spot rates.

## Details

Since \$\$ (1+z_n)^n = \prod\_{j=0}^{n-1}(1+f\_{j,1}), \$\$ the spot
rates are recovered directly.

## Examples

``` r
z_from_fn1(c(0.04, 0.05, 0.06, 0.07, 0.08))
#> [1] 0.04000000 0.04498804 0.04996825 0.05494075 0.05990565
```
