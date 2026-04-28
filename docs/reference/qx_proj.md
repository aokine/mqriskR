# Project one-year death probability under mortality improvement

Computes projected one-year death probability \$\$q_x^{\[Y\]} =
q_x^{\[B\]} (1-AA_x)^{Y-B}\$\$.

## Usage

``` r
qx_proj(qx_base, AAx, base_year, proj_year)
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

Numeric vector of projected one-year death probabilities.
