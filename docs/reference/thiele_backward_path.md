# Backward reserve path from a terminal value

Starting from a terminal reserve at the final time, computes reserves
backward over a strictly increasing time grid using
\`thiele_backward_step()\`.

## Usage

``` r
thiele_backward_path(times, V_terminal, P, delta, mu, benefit = 1)
```

## Arguments

- times:

  Finite numeric vector of strictly increasing times.

- V_terminal:

  Finite scalar reserve at the final time.

- P:

  Premium rate, scalar or vector with one value per time step.

- delta:

  Force of interest, scalar or vector with one value per step.

- mu:

  Nonnegative force of mortality, scalar or vector with one value per
  step.

- benefit:

  Benefit amount, scalar or vector with one value per step.

## Value

Numeric vector of reserve values corresponding to \`times\`.

## Examples

``` r
times <- seq(19, 20, by = 0.25)

thiele_backward_path(
  times,
  V_terminal = 1000,
  P = 26.96,
  delta = 0.058,
  mu = 0.002,
  benefit = 1000
)
#> [1]  918.1329  938.1449  958.4570  979.0739 1000.0000
```
