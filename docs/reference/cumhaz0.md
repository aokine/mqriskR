# Cumulative hazard for age-at-failure T0

Computes \\\Lambda_0(t)=\int_0^t \lambda_0(y)\\dy\\.

## Usage

``` r
cumhaz0(t, model, ...)
```

## Arguments

- t:

  Numeric vector of times (\\t \ge 0\\).

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

Numeric vector of cumulative hazard values.
