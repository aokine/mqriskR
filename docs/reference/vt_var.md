# Discount factors under a variable annual interest scenario

Computes the sequence of discount factors \$\$v_1,\\ v_2,\\ \dots,\\
v_n\$\$ where \$\$v_t = \prod\_{k=1}^{t}(1+i_k)^{-1}.\$\$

## Usage

``` r
vt_var(i)
```

## Arguments

- i:

  Numeric vector of annual effective interest rates \\i_1, i_2, \dots,
  i_n\\.

## Value

Numeric vector of discount factors of the same length as `i`.

## Details

This corresponds to the Chapter 15 notation \\{}\_j v^t\\ for a fixed
scenario \\j\\.

## Examples

``` r
vt_var(c(0.06, 0.07, 0.08))
#> [1] 0.9433962 0.8816787 0.8163692
```
