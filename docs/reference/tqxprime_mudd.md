# Independent probabilities \\{}\_t q_x^{\prime(j)}\\ under MUDD

Independent probabilities \\{}\_t q_x^{\prime(j)}\\ under MUDD

## Usage

``` r
tqxprime_mudd(qxj, t)
```

## Arguments

- qxj:

  Numeric vector of dependent probabilities \\q_x^{(j)}\\.

- t:

  Time in \[0,1\].

## Value

Numeric vector.

## Examples

``` r
tqxprime_mudd(c(.20, .10), t = 0.5)
#> [1] 0.10268289 0.05273176
```
