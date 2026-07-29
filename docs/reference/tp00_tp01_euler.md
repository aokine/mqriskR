# Euler approximation of disability-state probabilities

Approximates probabilities of being healthy, disabled, or deceased in a
three-state model that allows recovery from disability. The final time
is always included, with a shorter final step when needed.

## Usage

``` r
tp00_tp01_euler(h, n, mu01, mu02, mu10, mu12, p00_0 = 1, p01_0 = 0)
```

## Arguments

- h:

  Positive step size.

- n:

  Nonnegative final time.

- mu01:

  Healthy-to-disabled intensity function.

- mu02:

  Healthy-to-deceased intensity function.

- mu10:

  Disabled-to-healthy intensity function.

- mu12:

  Disabled-to-deceased intensity function.

- p00_0:

  Initial healthy-state probability.

- p01_0:

  Initial disabled-state probability.

## Value

A data frame with columns `t`, `tp00`, `tp01`, and `tp02`.

## Examples

``` r
mu01 <- function(t) 0.10 * t + 0.20
mu02 <- function(t) 0.20
mu10 <- function(t) 0.50
mu12 <- function(t) 0.125 * t + 0.20
tp00_tp01_euler(0.10, 2, mu01, mu02, mu10, mu12)
#>      t      tp00       tp01       tp02
#> 1  0.0 1.0000000 0.00000000 0.00000000
#> 2  0.1 0.9600000 0.02000000 0.02000000
#> 3  0.2 0.9216400 0.03873500 0.03962500
#> 4  0.3 0.8848679 0.05620279 0.05892934
#> 5  0.4 0.8496287 0.07240980 0.07796151
#> 6  0.5 0.8158655 0.08737015 0.09676433
#> 7  0.6 0.7835201 0.10110482 0.11537511
#> 8  0.7 0.7525334 0.11364071 0.13382589
#> 9  0.8 0.7228464 0.12500991 0.15214373
#> 10 0.9 0.6944002 0.13524881 0.17035095
#> 11 1.0 0.6671371 0.14439746 0.18846548
#> 12 1.1 0.6410001 0.15249878 0.20650114
#> 13 1.2 0.6159340 0.15959801 0.22446798
#> 14 1.3 0.5918853 0.16574206 0.24237259
#> 15 1.4 0.5688025 0.17097903 0.26021845
#> 16 1.5 0.5466361 0.17535765 0.27800621
#> 17 1.6 0.5253390 0.17892692 0.29573404
#> 18 1.7 0.5048664 0.18173571 0.31339790
#> 19 1.8 0.4851758 0.18383238 0.33099182
#> 20 1.9 0.4662272 0.18526456 0.34850822
#> 21 2.0 0.4479830 0.18607887 0.36593809
```
