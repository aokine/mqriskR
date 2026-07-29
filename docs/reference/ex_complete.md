# Complete expectation of life

Computes \\\overset{\circ}{e}\_x=\int_0^\infty {}\_t p_x\\dt\\.

## Usage

``` r
ex_complete(x, model, ..., tol = 1e-10)
```

## Arguments

- x:

  Numeric vector of ages.

- model:

  One of `"uniform"`, `"exponential"`, `"gompertz"`, `"makeham"`, or
  `"weibull"`.

- ...:

  Model parameters.

- tol:

  Tolerance used to choose a finite integration bound.

## Value

Numeric vector of complete expectations.
