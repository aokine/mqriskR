# Last-survivor temporary annuity-immediate

Computes \\a\_{\overline{xy}:\overline{n}\|}\\.

## Usage

``` r
axybarn(x, y, n, i, tbl = NULL, model = NULL, ...)
```

## Arguments

- x:

  Age of first life.

- y:

  Age of second life.

- n:

  Term.

- i:

  Effective annual interest rate.

- tbl:

  Life table.

- model:

  Survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

## Details

Shared documentation topic used to avoid filename collisions with
case-distinct function names on case-insensitive file systems.

## Examples

``` r
axybarn(40, 50, n = 10, i = 0.05, model = "uniform", omega = 100)
#> [1] 7.452381
```
