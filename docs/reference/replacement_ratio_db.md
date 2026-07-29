# Replacement ratio for a defined benefit plan

Computes annual benefit divided by a selected salary measure. Scalar
arguments are recycled to a common length.

## Usage

``` r
replacement_ratio_db(benefit, salary)
```

## Arguments

- benefit:

  Nonnegative annual retirement benefit.

- salary:

  Positive salary measure used in the denominator.

## Value

A numeric vector.

## Examples

``` r
replacement_ratio_db(
  benefit = 108008.66,
  salary = 187119.09
)
#> [1] 0.5772188
```
