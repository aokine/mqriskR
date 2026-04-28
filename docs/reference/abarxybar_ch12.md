# Continuous last-survivor whole life annuity

Computes \\\overline{a}\_{\overline{xy}} = \overline{a}\_x +
\overline{a}\_y - \overline{a}\_{xy}\\.

## Usage

``` r
abarxybar(x, y, i, tbl = NULL, model = NULL, ...)
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

## Value

Numeric vector.

## Details

Shared documentation topic used to avoid filename collisions with
case-distinct function names on case-insensitive file systems.

## Examples

``` r
abarxybar(40, 50, i = 0.05, model = "uniform", omega = 100)
#> [1] 16.24185
```
