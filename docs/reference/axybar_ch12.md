# Last-survivor whole life annuity-immediate

Computes \\a\_{\overline{xy}}\\.

## Usage

``` r
axybar(x, y, i, tbl = NULL, model = NULL, ..., k_max = 5000, tol = 1e-12)
```

## Arguments

- x:

  Age of first life.

- y:

  Age of second life.

- i:

  Effective annual interest rate.

- tbl:

  Life table.

- model:

  Survival model.

- ...:

  Additional model parameters.

- k_max:

  Maximum number of terms.

- tol:

  Convergence tolerance.

## Value

Numeric vector.

## Details

Shared documentation topic used to avoid filename collisions with
case-distinct function names on case-insensitive file systems.

## Examples

``` r
axybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 15.74572
```
