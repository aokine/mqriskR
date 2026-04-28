# Conditional failure probability for Tx

Computes \\{}\_t q_x = 1 - {}\_t p_x\\.

## Usage

``` r
tqx(t, x, model, ...)
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

Numeric vector in \\\[0,1\]\\.
