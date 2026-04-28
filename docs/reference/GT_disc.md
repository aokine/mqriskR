# Total gain for a discrete insurance contract

Computes the Chapter 10 total gain: amount on hand at year-end minus
amount required.

## Usage

``` r
GT_disc(Vt, Vt1, P, i_actual, q_actual, B = 1)
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

- q_actual:

  Actual mortality rate for the year.

- B:

  Benefit amount. Defaults to 1.

## Value

Numeric vector.

## Examples

``` r
GT_disc(Vt = 0.1, Vt1 = 0.11, P = 0.02, i_actual = 0.05, q_actual = 0.01)
#> [1] 0.0071
```
