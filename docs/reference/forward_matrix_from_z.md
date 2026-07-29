# Matrix of forward rates implied by spot rates

Constructs a matrix of annual effective forward rates. Rows correspond
to forward-start times \\n=1,\ldots,m-1\\; columns correspond to forward
periods \\k=1,\ldots,m-1\\. Entries requiring maturities beyond \\m\\
are returned as `NA`.

## Usage

``` r
forward_matrix_from_z(z)
```

## Arguments

- z:

  Numeric vector of at least two annual effective spot rates. Each value
  must be greater than `-1`.

## Value

A numeric matrix.

## Examples

``` r
forward_matrix_from_z(c(0.03, 0.04, 0.05, 0.06, 0.07))
#>            k=1        k=2        k=3       k=4
#> n=1 0.05009709 0.06014516 0.07019293 0.0802404
#> n=2 0.07028939 0.08038462 0.09047924        NA
#> n=3 0.09057507 0.10071655         NA        NA
#> n=4 0.11095234         NA         NA        NA
```
