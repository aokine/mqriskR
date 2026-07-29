# Cause-specific numbers of decrements

Computes \$\$d_x^{(j)}=l_x^{(\tau)}q_x^{(j)}.\$\$

## Usage

``` r
dxj(lxtau, qxj)
```

## Arguments

- lxtau:

  Number alive at age or duration \\x\\ in the multiple-decrement table.

- qxj:

  Numeric vector of cause-specific decrement probabilities.

## Value

A numeric vector containing the number of decrements from each cause.

## Examples

``` r
dxj(1000, c(0.011, 0.100))
#> [1]  11 100
```
