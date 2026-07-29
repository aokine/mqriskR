# Total one-year survival probability

Computes \$\$p_x^{(\tau)}=1-q_x^{(\tau)}.\$\$

## Usage

``` r
pxtau(qxj)
```

## Arguments

- qxj:

  Numeric vector of cause-specific decrement probabilities.

## Value

A numeric scalar.

## Examples

``` r
pxtau(c(0.011, 0.100))
#> [1] 0.889
```
