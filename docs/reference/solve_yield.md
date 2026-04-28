# Solve the yield rate by the equation of value

Finds the interest rate i such that the present value of the cash flows
is 0.

## Usage

``` r
solve_yield(cf, t, interval = c(-0.99, 1), tol = 1e-10)
```

## Arguments

- cf:

  Cash flows.

- t:

  Times.

- interval:

  Two-length numeric vector bracketing the root.

- tol:

  Tolerance passed to uniroot.

## Value

Yield rate i.

## Examples

``` r
solve_yield(c(-100, 60, 60), c(0, 1, 2), interval = c(-0.5, 1))
#> [1] 0.1306624
```
