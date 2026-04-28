# Complete expectation of life at age x

Computes \\\overset{\circ}{e}\_x=\int_0^\infty {}\_t p_x \\ dt\\. Uses
closed form for uniform and exponential; numeric integration otherwise.

## Usage

``` r
ex_complete(x, model, ..., tol = 1e-10)
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

- tol:

  Tolerance used to choose a finite integration bound (numeric).

## Value

Numeric vector of complete expectations.
