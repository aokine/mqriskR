# Multi-year survival probability under mortality improvement

Computes survival over `n` years starting at age `x0` in year
`issue_year`, using base-year mortality `qx_base_vec` and improvement
factors `AAx_vec`.

## Usage

``` r
tpx_improved(x0, n, qx_base_vec, AAx_vec, base_year, issue_year)
```

## Arguments

- x0:

  Issue age.

- n:

  Number of years.

- qx_base_vec:

  Base-year one-year death probabilities for successive ages.

- AAx_vec:

  Mortality improvement factors for successive ages.

- base_year:

  Base year.

- issue_year:

  Issue year.

## Value

Numeric scalar.

## Details

The vectors should correspond to ages `x0, x0+1, ..., x0+n-1`.

The survival probability is \$\$\prod\_{j=0}^{n-1} \left(1 -
q\_{x+j}^{\[issue\\year+j\]}\right)\$\$.
