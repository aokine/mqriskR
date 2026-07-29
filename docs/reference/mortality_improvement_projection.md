# Mortality improvement projection functions

Functions for projecting one-year death and survival probabilities under
mortality improvement and evaluating annuity values under projected
mortality.

## Usage

``` r
qx_proj(qx_base, AAx, base_year, proj_year)

px_proj(qx_base, AAx, base_year, proj_year)

tpx_improved(x0, n, qx_base_vec, AAx_vec, base_year, issue_year)

axn_improved(x0, n, i, qx_base_vec, AAx_vec, base_year, issue_year)

naxn_improved(x0, u, n, i, qx_base_vec, AAx_vec, base_year, issue_year)

ax_improved(x0, i, qx_base_vec, AAx_vec, base_year, issue_year)
```

## Arguments

- qx_base:

  Base-year one-year death probability.

- AAx:

  Mortality improvement factor.

- base_year:

  Base year.

- proj_year:

  Projection year. May be scalar or vector.

- x0:

  Issue age.

- n:

  Number of years.

- qx_base_vec:

  Base-year death probabilities for successive ages.

- AAx_vec:

  Mortality improvement factors for successive ages.

- issue_year:

  Issue year.

- i:

  Effective annual interest rate.

- u:

  Deferral period in years.

## Value

Numeric vector of projected one-year death probabilities.

## Details

The standard projection used is \$\$q_x^{\[Y\]} = q_x^{\[B\]} (1 -
AA_x)^{Y-B},\$\$ where \\B\\ is the base year, \\Y\\ is the projection
year, and \\AA_x\\ is the mortality improvement factor at age \\x\\.
