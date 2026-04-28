# Internal: choose a practical upper bound for integrals to infinity

Strategy: find t such that S0(t) is very small, or fall back to a cap.

## Usage

``` r
.find_upper_t(x, model, ..., eps = 1e-10, t_max = 500)
```
