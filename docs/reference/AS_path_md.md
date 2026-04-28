# General projected asset share path (multiple decrements)

General projected asset share path (multiple decrements)

## Usage

``` r
AS_path_md(AS0, G, r, e, b_mat, q_mat, p_tau, i, b_surv = NULL)
```

## Arguments

- AS0:

  Initial asset share.

- G:

  Premium.

- r:

  Expense percentages.

- e:

  Fixed expenses.

- b_mat:

  Matrix of benefits (rows = years, cols = causes).

- q_mat:

  Matrix of decrement probabilities (same shape as b_mat).

- p_tau:

  In-force probabilities.

- i:

  Interest rate.

- b_surv:

  Optional survival benefits.

## Value

A data frame with columns `k` and `AS`. Column `k` gives the policy year
from 0 to `n`, and column `AS` gives the corresponding projected asset
share at each year.
