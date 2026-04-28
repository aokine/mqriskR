# Matrix of all determinable forward rates from spot rates

Constructs an upper-left triangular matrix of annual effective forward
rates \\f\_{n,k}\\ implied by annual effective spot rates
\\z_1,\dots,z_m\\.

## Usage

``` r
forward_matrix_from_z(z)
```

## Arguments

- z:

  Numeric vector of annual effective spot rates.

## Value

A numeric matrix.

## Details

Rows correspond to \\n = 1,\dots,m-1\\ and columns correspond to \\k =
1,\dots,m-1\\. Entries that are not determinable are returned as `NA`.

## Examples

``` r
z <- c(0.03, 0.04, 0.05, 0.06, 0.07)
forward_matrix_from_z(z)
#>            k=1        k=2        k=3       k=4
#> n=1 0.05009709 0.06014516 0.07019293 0.0802404
#> n=2 0.07028939 0.08038462 0.09047924        NA
#> n=3 0.09057507 0.10071655         NA        NA
#> n=4 0.11095234         NA         NA        NA
```
