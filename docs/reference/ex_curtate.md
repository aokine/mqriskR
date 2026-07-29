# Curtate expectation of life

Computes \\e_x = E\[K_x\] = \sum\_{k=1}^{\infty} {}\_k p_x\\.

## Usage

``` r
ex_curtate(x, model, ..., k_max = 5000, tol = 1e-12)
```

## Arguments

- x:

  Numeric vector of ages.

- model:

  One of `"uniform"`, `"exponential"`, `"gompertz"`, `"makeham"`, or
  `"weibull"`.

- ...:

  Model parameters.

- k_max:

  Maximum integer duration to sum to.

- tol:

  Stop early if the summand is smaller than this tolerance for several
  consecutive steps.

## Value

Numeric vector of curtate expectations.
