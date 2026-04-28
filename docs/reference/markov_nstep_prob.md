# n-step transition probability for a discrete-time Markov chain

Computes the \\(i,j)\\ entry of \\P^n\\, useful for Chapter 14 examples
involving discrete-time multi-state models such as CCRC and risk-class
models.

## Usage

``` r
markov_nstep_prob(P, n, i, j)
```

## Arguments

- P:

  Transition probability matrix.

- n:

  Nonnegative integer number of steps.

- i:

  Starting state index.

- j:

  Ending state index.

## Value

A numeric scalar.

## Examples

``` r
P <- matrix(
  c(0.94, 0.03, 0.02, 0.01,
    0.50, 0.30, 0.18, 0.02,
    0.00, 0.00, 0.93, 0.07,
    0.00, 0.00, 0.00, 1.00),
  nrow = 4, byrow = TRUE
)

markov_nstep_prob(P, n = 3, i = 1, j = 1)
#> [1] 0.863284
markov_nstep_prob(P, n = 3, i = 1, j = 3)
#> [1] 0.064472
```
