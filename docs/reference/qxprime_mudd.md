# Associated single-decrement probabilities under MUDD

Converts dependent multiple-decrement probabilities to associated
single-decrement probabilities under the multiple-decrement uniform
distribution assumption.

## Usage

``` r
qxprime_mudd(qxj)
```

## Arguments

- qxj:

  Numeric vector of dependent cause-specific decrement probabilities.

## Value

A numeric vector of associated single-decrement probabilities.

## Examples

``` r
qxprime_mudd(c(0.20, 0.10))
#> [1] 0.2116265 0.1120960
```
