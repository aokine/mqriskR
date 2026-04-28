# Deferred temporary annuity-immediate under mortality improvement

Computes \$\${}\_{u\mid}a\_{x:\overline{n}\|} = \sum\_{t=u+1}^{u+n} v^t
\\ {}\_tp_x^{(\mathrm{improved})}\$\$

## Usage

``` r
naxn_improved(x0, u, n, i, qx_base_vec, AAx_vec, base_year, issue_year)
```

## Arguments

- x0:

  Issue age.

- u:

  Deferral period in years.

- n:

  Temporary payment period in years.

- i:

  Effective annual interest rate.

- qx_base_vec:

  Base-year one-year death probabilities for ages
  `x0, x0+1, ..., x0+u+n-1`.

- AAx_vec:

  Mortality improvement factors for ages `x0, x0+1, ..., x0+u+n-1`.

- base_year:

  Base year.

- issue_year:

  Issue year.

## Value

Numeric scalar.
