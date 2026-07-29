# Multi-step transition probability

Computes an entry of the matrix power \\P^n\\.

## Usage

``` r
markov_nstep_prob(P, n, i, j)
```

## Arguments

- P:

  Square transition-probability matrix.

- n:

  Nonnegative integer number of steps.

- i:

  Starting-state index.

- j:

  Ending-state index.

## Value

A numeric scalar.

## Examples

``` r
P <- matrix(c(0.9, 0.1, 0, 1), nrow = 2, byrow = TRUE)
markov_nstep_prob(P, n = 3, i = 1, j = 2)
#> [1] 0.271
```
