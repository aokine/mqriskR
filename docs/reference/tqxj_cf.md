# Cause-specific decrement probability under constant forces

Computes \$\$ {}\_tq_x^{(j)} = \frac{\mu_j}{\sum_k\mu_k} \left\[
1-\exp\left(-t\sum_k\mu_k\right) \right\]. \$\$

## Usage

``` r
tqxj_cf(mu, j, t)
```

## Arguments

- mu:

  Numeric vector of nonnegative cause-specific forces.

- j:

  Positive integer cause index.

- t:

  Nonnegative time. May be scalar or vector.

## Value

A numeric vector.

## Examples

``` r
tqxj_cf(c(0.10, 0.20), j = 1, t = 5)
#> [1] 0.2589566
```
