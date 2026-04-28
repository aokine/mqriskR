# Conditional density for Tx

Computes \\f_x(t)=f_0(x+t)/S_0(x)\\.

## Usage

``` r
fx(t, x, model, ...)
```

## Arguments

- t:

  Numeric vector of durations (\\t \ge 0\\).

- x:

  Numeric vector of ages (\\x \ge 0\\).

- model:

  One of `"uniform"`, `"exponential"`, `"gompertz"`, `"makeham"`,
  `"weibull"`.

- ...:

  Model parameters:

  - uniform: `omega`

  - exponential: `lambda`

  - gompertz: `B`, `c`

  - makeham: `A`, `B`, `c`

  - weibull: `shape`, `scale`

## Value

Numeric vector of densities (\\\ge 0\\).

## Details

Vectorization rule: - If `t` and `x` are the same length, values are
computed elementwise. - If one of `t` or `x` has length 1, it is
recycled to match the other.
