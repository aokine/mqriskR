# One backward numerical step for Thiele's equation

Approximates the reserve at time \`t\` from a known reserve at time
\`t + h\`.

## Usage

``` r
thiele_backward_step(V_next, P, delta, mu, benefit = 1, h = 1)
```

## Arguments

- V_next:

  Reserve at time \`t + h\`.

- P:

  Premium rate.

- delta:

  Force of interest.

- mu:

  Nonnegative force of mortality at time \`t\`.

- benefit:

  Benefit amount.

- h:

  Positive step size.

## Value

Numeric vector of reserve approximations.

## Details

The implemented step is \$\$ V_t = \frac{ V\_{t+h} - hP + h\mu B }{ 1 +
h(\delta + \mu) }. \$\$

## Examples

``` r
thiele_backward_step(
  V_next = 1000,
  P = 26.96,
  delta = 0.058,
  mu = 0.002,
  benefit = 1000,
  h = 1
)
#> [1] 919.8491
```
