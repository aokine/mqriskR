# Forward rate implied by spot rates

Computes the annual effective forward rate for the interval from time
`n` to time `n + k`: \$\$ (1+z\_{n+k})^{n+k} = (1+z_n)^n(1+f\_{n,k})^k.
\$\$

## Usage

``` r
fnk_from_z(z, n, k)
```

## Arguments

- z:

  Numeric vector of annual effective spot rates for maturities
  `1, ..., length(z)`. Each value must be greater than `-1`.

- n:

  Nonnegative integer forward-start time.

- k:

  Positive integer forward period.

## Value

A numeric scalar.

## Examples

``` r
z <- c(0.03, 0.04, 0.05, 0.06, 0.07)
fnk_from_z(z, n = 1, k = 4)
#> [1] 0.0802404
fnk_from_z(z, n = 2, k = 2)
#> [1] 0.08038462
```
