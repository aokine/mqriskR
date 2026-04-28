# Mortality gain for a discrete insurance contract

Mortality gain for a discrete insurance contract

## Usage

``` r
GM_disc(Vt, Vt1, P, i_assumed, q_actual, B = 1)
```

## Arguments

- Vt:

  Reserve at duration t.

- Vt1:

  Reserve at duration t+1.

- P:

  Net premium for the year.

- i_assumed:

  Assumed annual effective interest rate.

- q_actual:

  Actual mortality rate for the year.

- B:

  Benefit amount. Defaults to 1.

## Value

Numeric vector.

## Examples

``` r
GM_disc(Vt = 0.1, Vt1 = 0.11, P = 0.02, i_assumed = 0.04, q_actual = 0.01)
#> [1] 0.0059
```
