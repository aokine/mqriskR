# Present value of cash flows at time 0

Present value of cash flows at time 0

## Usage

``` r
pv_cashflows(cf, t, i)
```

## Arguments

- cf:

  Cash flow amounts (positive = inflow, negative = outflow).

- t:

  Times of cash flows (same length as cf).

- i:

  Effective interest rate.

## Value

Present value at time 0.

## Examples

``` r
pv_cashflows(c(-100, 60, 60), c(0, 1, 2), i = 0.10)
#> [1] 4.132231
```
