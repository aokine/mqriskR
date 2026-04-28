# Backward Euler reserve path from maturity

Starting from a terminal reserve value at time T, computes reserves
backward on a grid using the backward Euler-style Thiele step.

## Usage

``` r
thiele_backward_path(times, V_terminal, P, delta, mu, benefit = 1)
```

## Arguments

- times:

  Vector of times in increasing order.

- V_terminal:

  Reserve at the final time.

- P:

  Premium rate, scalar or vector of length length(times)-1.

- delta:

  Force of interest, scalar or vector of length length(times)-1.

- mu:

  Force of mortality, scalar or vector of length length(times)-1.

- benefit:

  Benefit amount, scalar or vector of length length(times)-1.

## Value

Numeric vector of reserve values on the grid.

## Examples

``` r
times <- seq(19, 20, by = 0.25)
thiele_backward_path(times, V_terminal = 1000, P = 26.96, delta = 0.058, mu = 0.002, benefit = 1000)
#> [1]  918.1329  938.1449  958.4570  979.0739 1000.0000
```
