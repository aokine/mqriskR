# Profit signature

Converts policy-year expected profits into a profit signature by
weighting each future expected profit by the probability that the
contract is in force at the start of that policy year.

## Usage

``` r
Pi_signature(Pr, p_tau)
```

## Arguments

- Pr:

  Profit vector of length `n + 1`.

- p_tau:

  One-year in-force probabilities. For `n > 1`, this may have length
  `n - 1` or `n`; the final value is ignored when length `n`. For a
  one-year contract, use `numeric(0)` or a single probability.

## Value

A named numeric vector with the same length as `Pr`.

## Examples

``` r
Pr <- c(-15.00, 8.42, 8.40, 8.61)
Pi_signature(Pr, p_tau = c(0.99858, 0.99847, 0.99834))
#>        Pi0        Pi1        Pi2        Pi3 
#> -15.000000   8.420000   8.388072   8.584619 
```
