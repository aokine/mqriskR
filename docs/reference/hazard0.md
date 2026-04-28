# Hazard / force for age-at-failure T0

Computes \\\lambda_0(t)=f_0(t)/S_0(t)\\.

## Usage

``` r
hazard0(t, model, ...)
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

Numeric vector of hazard/force values.
