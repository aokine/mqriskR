# Mean present value of loss at duration t for whole life insurance

Computes the Chapter 10 conditional mean \\E\[{}\_tL_x \mid K_x \ge
t\]\\ for a fully discrete whole life insurance.

## Usage

``` r
ELtx(x, t, i, P, model, ...)
```

## Arguments

- x:

  Issue age.

- t:

  Duration.

- i:

  Effective annual interest rate.

- P:

  Annual premium.

- model:

  Survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.

## Details

Under the equivalence-principle premium, this equals the prospective
reserve \\{}\_tV_x\\.

## Examples

``` r
prem <- Px(40, i = 0.05, model = "uniform", omega = 100)
ELtx(40, t = 10, i = 0.05, P = prem, model = "uniform", omega = 100)
#> [1] 0.07250474
```
