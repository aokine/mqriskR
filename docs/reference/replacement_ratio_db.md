# Replacement ratio for a defined benefit plan

Computes a DB replacement ratio as benefit divided by a chosen salary
measure.

## Usage

``` r
replacement_ratio_db(benefit, salary)
```

## Arguments

- benefit:

  Annual retirement benefit.

- salary:

  Salary measure used in the denominator.

## Value

Replacement ratio.

## Examples

``` r
replacement_ratio_db(benefit = 108008.66, salary = 187119.09)
#> [1] 0.5772188
```
