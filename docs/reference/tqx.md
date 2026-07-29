# Conditional failure probability

Computes \\{}\_t q_x = 1 - {}\_t p_x\\.

## Usage

``` r
tqx(t, x, model, ...)
```

## Arguments

- t:

  Numeric vector of durations.

- x:

  Numeric vector of ages.

- model:

  One of `"uniform"`, `"exponential"`, `"gompertz"`, `"makeham"`, or
  `"weibull"`.

- ...:

  Model parameters.

## Value

Numeric vector of failure probabilities.
