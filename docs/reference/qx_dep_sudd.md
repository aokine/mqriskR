# Dependent probabilities \\q_x^{(j)}\\ from independent probabilities \\q_x^{\prime(j)}\\ under SUDD

Two-decrement case only.

## Usage

``` r
qx_dep_sudd(q1prime, q2prime)
```

## Arguments

- q1prime:

  Independent probability for decrement 1.

- q2prime:

  Independent probability for decrement 2.

## Value

Numeric vector `c(q1, q2)`.

## Examples

``` r
qx_dep_sudd(0.20, 0.10)
#>   q1   q2 
#> 0.19 0.09 
```
