# Single-decrement survival under a constant force

Computes \$\${}\_tp_x^{\prime(j)}=\exp(-\mu_jt).\$\$

## Usage

``` r
tpxprimej_cf(mu, t)
```

## Arguments

- mu:

  Nonnegative constant force of decrement. May be scalar or vector.

- t:

  Nonnegative time. May be scalar or vector.

## Value

A numeric vector.

## Examples

``` r
tpxprimej_cf(0.10, 5)
#> [1] 0.6065307
tpxprimej_cf(0.10, c(1, 5, 10))
#> [1] 0.9048374 0.6065307 0.3678794
```
