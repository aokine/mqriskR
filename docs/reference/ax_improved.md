# Whole life annuity-immediate under mortality improvement

Computes a truncated whole life annuity-immediate under projected
mortality.

## Usage

``` r
ax_improved(x0, i, qx_base_vec, AAx_vec, base_year, issue_year)
```

## Arguments

- x0:

  Issue age.

- i:

  Effective annual interest rate.

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
