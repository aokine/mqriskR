# Curtate expectation of life at age x

Computes \\e_x = E\[K_x\] = \sum\_{k=1}^\infty {}\_k p_x\\ with
truncation.

## Usage

``` r
ex_curtate(x, model, ..., k_max = 5000, tol = 1e-12)
```

## Arguments

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

- k_max:

  Maximum integer duration to sum to.

- tol:

  Stop early if the summand is \< tol for several steps.

## Value

Numeric vector of curtate expectations.
