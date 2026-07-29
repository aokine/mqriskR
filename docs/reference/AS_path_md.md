# General projected asset-share path

Computes projected asset shares for a multiple-decrement contract. Rows
of `b_mat` and `q_mat` represent policy years and columns represent
decrement causes.

## Usage

``` r
AS_path_md(AS0, G, r, e, b_mat, q_mat, p_tau, i, b_surv = NULL)
```

## Arguments

- AS0:

  Initial asset share.

- G:

  Premium amount by policy year.

- r:

  Percent-of-premium expense rate by policy year.

- e:

  Fixed expense by policy year.

- b_mat:

  Matrix of decrement benefits.

- q_mat:

  Matrix of decrement probabilities.

- p_tau:

  In-force probability by policy year.

- i:

  Effective annual interest rate by policy year.

- b_surv:

  Survival benefit by policy year.

## Value

A data frame with policy year `k` and asset share `AS`.
