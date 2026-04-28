# Temporary annuity-immediate under mortality improvement

Computes \$\$a\_{x:\overline{n}\|} = \sum\_{t=1}^n v^t \\
{}\_tp_x^{(\mathrm{improved})}\$\$ using projected mortality rates.

## Usage

``` r
axn_improved(x0, n, i, qx_base_vec, AAx_vec, base_year, issue_year)
```

## Arguments

- x0:

  Issue age.

- n:

  Term in years.

- i:

  Effective annual interest rate.

- qx_base_vec:

  Base-year one-year death probabilities for ages
  `x0, x0+1, ..., x0+n-1`.

- AAx_vec:

  Mortality improvement factors for ages `x0, x0+1, ..., x0+n-1`.

- base_year:

  Base year.

- issue_year:

  Issue year.

## Value

Numeric scalar.
