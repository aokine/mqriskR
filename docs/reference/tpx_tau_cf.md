# Total survival under constant cause-specific forces

Computes \$\$ {}\_tp_x^{(\tau)} = \exp\left(-t\sum_j\mu_j\right). \$\$

## Usage

``` r
tpx_tau_cf(mu, t)
```

## Arguments

- mu:

  Numeric vector of nonnegative cause-specific forces.

- t:

  Nonnegative time. May be scalar or vector.

## Value

A numeric vector.

## Examples

``` r
tpx_tau_cf(c(0.10, 0.20), 5)
#> [1] 0.2231302
tpx_tau_cf(c(0.10, 0.20), c(1, 5, 10))
#> [1] 0.74081822 0.22313016 0.04978707
```
