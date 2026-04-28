# Present value of a level annuity-certain

Computes the present value of an n-period annuity certain with level
payments of 1 per period.

## Usage

``` r
annuity_certain(n, i, due = FALSE, m = 1, cont = FALSE)
```

## Arguments

- n:

  Number of payments or periods.

- i:

  Effective interest rate per period.

- due:

  If TRUE, annuity-due; otherwise annuity-immediate.

- m:

  Payment frequency per period (m = 1 means annual).

- cont:

  If TRUE, continuous payment model.

## Value

Present value.

## Examples

``` r
annuity_certain(n = 10, i = 0.05)
#> [1] 7.721735
annuity_certain(n = 10, i = 0.05, due = TRUE)
#> [1] 8.107822
annuity_certain(n = 10, i = 0.05, cont = TRUE)
#> [1] 7.913209
```
