# Profit signature from a profit vector

Converts a Chapter 17 profit vector \\\mathbf{Pr}=(Pr_0,\dots,Pr_n)\\
into the corresponding profit signature
\\\mathbf{\Pi}=(\Pi_0,\dots,\Pi_n)\\ using Equation (17.3).

## Usage

``` r
Pi_signature(Pr, p_tau)
```

## Arguments

- Pr:

  Profit vector of length \\n+1\\.

- p_tau:

  One-year in-force probabilities. This may have length \\n-1\\ or
  \\n\\. If length \\n\\, the final entry is ignored for the
  profit-signature calculation.

## Value

Numeric vector of length \\n+1\\.

## Examples

``` r
Pr <- c(-15.00, 8.42, 8.40, 8.61)
Pi_signature(Pr, p_tau = c(0.99858, 0.99847, 0.99834))
#>        Pi0        Pi1        Pi2        Pi3 
#> -15.000000   8.420000   8.388072   8.584619 
```
