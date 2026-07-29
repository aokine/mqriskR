# Total number of decrements

Computes \$\$d_x^{(\tau)}=\sum_j d_x^{(j)}.\$\$

## Usage

``` r
dxtau(lxtau, qxj)
```

## Arguments

- lxtau:

  Number alive at age or duration \\x\\ in the multiple-decrement table.

- qxj:

  Numeric vector of cause-specific decrement probabilities.

## Value

A numeric scalar.

## Examples

``` r
dxtau(1000, c(0.011, 0.100))
#> [1] 111
```
