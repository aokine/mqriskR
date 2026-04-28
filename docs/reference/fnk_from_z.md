# Forward rate \\f\_{n,k}\\ from spot rates

Computes the \\n\\-year forward \\k\\-year annual effective rate implied
by annual effective spot rates: \$\$ (1+z\_{n+k})^{n+k} = (1+z_n)^n
(1+f\_{n,k})^k. \$\$

## Usage

``` r
fnk_from_z(z, n, k)
```

## Arguments

- z:

  Numeric vector of annual effective spot rates.

- n:

  Forward start in years.

- k:

  Forward maturity in years.

## Value

A numeric scalar.

## Details

The input vector `z` should contain annual effective spot rates for
maturities 1, 2, ..., length(z).

## Examples

``` r
z <- c(0.03, 0.04, 0.05, 0.06, 0.07)
fnk_from_z(z, n = 1, k = 4)
#> [1] 0.0802404
fnk_from_z(z, n = 2, k = 2)
#> [1] 0.08038462
```
