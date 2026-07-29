# Multiple-decrement probabilities under constant forces

Converts associated single-decrement probabilities to dependent
cause-specific decrement probabilities under constant forces.

## Usage

``` r
qx_dep_cf(qxprime)
```

## Arguments

- qxprime:

  Numeric vector of associated single-decrement probabilities. Each
  value must be less than one.

## Value

A numeric vector of dependent cause-specific decrement probabilities.

## Examples

``` r
qx_dep_cf(c(0.20, 0.10))
#> [1] 0.19019611 0.08980389
```
