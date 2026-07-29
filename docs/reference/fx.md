# Conditional density

Computes \\f_x(t)=f_0(x+t)/S_0(x)\\.

## Usage

``` r
fx(t, x, model, ...)
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

Numeric vector of density values.
