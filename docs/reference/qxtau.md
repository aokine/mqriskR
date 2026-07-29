# Total one-year decrement probability

Computes the total one-year multiple-decrement probability
\$\$q_x^{(\tau)}=\sum_j q_x^{(j)}.\$\$

## Usage

``` r
qxtau(qxj)
```

## Arguments

- qxj:

  Numeric vector of cause-specific decrement probabilities.

## Value

A numeric scalar.

## Examples

``` r
qxtau(c(0.011, 0.100))
#> [1] 0.111
```
