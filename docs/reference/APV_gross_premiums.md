# APV of gross premiums under a risk discount rate

Computes the actuarial present value of gross premiums: \$\$ APV\_{GP} =
\sum\_{t=0}^{n-1} \frac{G\_{t+1} \cdot {}\_tp_x^{(\tau)}}{(1+r)^t}. \$\$

## Usage

``` r
APV_gross_premiums(G, r, p_tau)
```

## Arguments

- G:

  Gross premium vector for policy years 1 through \\n\\.

- r:

  Risk discount rate.

- p_tau:

  One-year in-force probabilities. This may have length \\n-1\\ or
  \\n\\. If length \\n\\, the final entry is ignored.

## Value

Numeric scalar.

## Examples

``` r
APV_gross_premiums(G = rep(95, 3),r = 0.10,p_tau = c(0.99858, 0.99847, 0.99834))
#> [1] 259.522
```
