# Reserve derivative from Thiele's equation

Computes \\dV/dt = P + \delta V - \mu(B - V)\\.

## Usage

``` r
thiele_dVdt(V, P, delta, mu, benefit = 1)
```

## Arguments

- V:

  Reserve at time t.

- P:

  Premium rate.

- delta:

  Force of interest.

- mu:

  Force of mortality.

- benefit:

  Benefit amount. Defaults to 1.

## Value

Numeric vector.

## Examples

``` r
thiele_dVdt(V = 900, P = 25, delta = 0.05, mu = 0.002, benefit = 1000)
#> [1] 69.8
```
