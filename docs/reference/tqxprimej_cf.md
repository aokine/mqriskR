# Single-decrement failure under a constant force

Computes \$\${}\_tq_x^{\prime(j)} =1-\exp(-\mu_jt).\$\$

## Usage

``` r
tqxprimej_cf(mu, t)
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
tqxprimej_cf(0.10, 5)
#> [1] 0.3934693
```
