# Associated single-decrement probabilities under SUDD

Converts two dependent multiple-decrement probabilities to associated
single-decrement probabilities under the single-decrement uniform
distribution assumption.

## Usage

``` r
qxprime_sudd(q1, q2)
```

## Arguments

- q1:

  Dependent decrement probability for cause 1. May be scalar or vector.

- q2:

  Dependent decrement probability for cause 2. May be scalar or vector.

## Value

For scalar input, a named numeric vector containing `q1prime` and
`q2prime`. For vectorized input, a numeric matrix with columns `q1prime`
and `q2prime`.

## Examples

``` r
qxprime_sudd(0.20, 0.10)
#>   q1prime   q2prime 
#> 0.2118473 0.1118473 
```
