# Project one-year survival probability under mortality improvement

Computes projected one-year survival probability \$\$p_x^{\[Y\]} = 1 -
q_x^{\[Y\]}\$\$.

## Usage

``` r
px_proj(qx_base, AAx, base_year, proj_year)
```

## Arguments

- qx_base:

  Base-year one-year death probability \\q_x^{\[B\]}\\.

- AAx:

  Mortality improvement factor \\AA_x\\.

- base_year:

  Base year \\B\\.

- proj_year:

  Projection year \\Y\\. May be scalar or vector.

## Value

Numeric vector of projected one-year survival probabilities.
