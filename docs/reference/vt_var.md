# Discount factors under variable annual interest rates

Computes cumulative discount factors for a sequence of annual effective
interest rates: \$\$ v_t = \prod\_{k=1}^{t}(1+i_k)^{-1}, \qquad
t=1,\ldots,n. \$\$

## Usage

``` r
vt_var(i)
```

## Arguments

- i:

  Numeric vector of annual effective interest rates. Each value must be
  greater than `-1`.

## Value

A numeric vector of cumulative discount factors with the same length as
`i`.

## Examples

``` r
vt_var(c(0.06, 0.07, 0.08))
#> [1] 0.9433962 0.8816787 0.8163692
vt_var(c(-0.01, 0.02, 0.03))
#> [1] 1.0101010 0.9902951 0.9614516
```
