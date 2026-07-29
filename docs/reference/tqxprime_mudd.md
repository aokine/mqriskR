# Fractional-year associated single-decrement probabilities under MUDD

Computes the associated single-decrement probabilities through time `t`
under the multiple-decrement uniform distribution assumption.

## Usage

``` r
tqxprime_mudd(qxj, t)
```

## Arguments

- qxj:

  Numeric vector of dependent cause-specific decrement probabilities.

- t:

  A single time in \\\[0,1\]\\.

## Value

A numeric vector of associated single-decrement probabilities.

## Examples

``` r
tqxprime_mudd(c(0.20, 0.10), t = 0.5)
#> [1] 0.10268289 0.05273176
```
