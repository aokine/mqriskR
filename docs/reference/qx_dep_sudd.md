# Multiple-decrement probabilities under SUDD

Converts two associated single-decrement probabilities to dependent
multiple-decrement probabilities under the single-decrement uniform
distribution assumption.

## Usage

``` r
qx_dep_sudd(q1prime, q2prime)
```

## Arguments

- q1prime:

  Associated single-decrement probability for cause 1. May be scalar or
  vector.

- q2prime:

  Associated single-decrement probability for cause 2. May be scalar or
  vector.

## Value

For scalar input, a named numeric vector containing `q1` and `q2`. For
vectorized input, a numeric matrix with columns `q1` and `q2`.

## Examples

``` r
qx_dep_sudd(0.20, 0.10)
#>   q1   q2 
#> 0.19 0.09 
```
