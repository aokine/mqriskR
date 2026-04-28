# Interest gain for a discrete insurance contract

Interest gain for a discrete insurance contract

## Usage

``` r
GI_disc(Vt, Vt1, P, i_actual, q_assumed, B = 1)
```

## Arguments

- Vt:

  Reserve at duration t.

- Vt1:

  Reserve at duration t+1.

- P:

  Net premium for the year.

- i_actual:

  Actual annual effective interest rate.

- q_assumed:

  Assumed mortality rate for the year.

- B:

  Benefit amount. Defaults to 1.

## Value

Numeric vector.

## Examples

``` r
GI_disc(Vt = 0.1, Vt1 = 0.11, P = 0.02, i_actual = 0.05, q_assumed = 0.01)
#> [1] 0.0071
```
