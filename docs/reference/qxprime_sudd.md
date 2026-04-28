# Independent probabilities \\q_x^{\prime(j)}\\ from dependent probabilities \\q_x^{(j)}\\ under SUDD

Two-decrement case only.

## Usage

``` r
qxprime_sudd(q1, q2)
```

## Arguments

- q1:

  Dependent probability for decrement 1.

- q2:

  Dependent probability for decrement 2.

## Value

Numeric vector `c(q1prime, q2prime)`.

## Examples

``` r
qxprime_sudd(0.20, 0.10)
#>   q1prime   q2prime 
#> 0.2118473 0.1118473 
```
