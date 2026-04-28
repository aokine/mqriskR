# Piecewise-continuous increasing whole life insurance

Computes \$\$(I\bar{A})\_x = \int_0^\infty \lfloor t+1 \rfloor\\ v^t\\
{}\_tp_x\\ \mu\_{x+t}\\ dt.\$\$

## Usage

``` r
IAbarx(x, i, model, ...)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- model:

  Parametric survival model.

- ...:

  Additional model parameters.

## Value

Numeric vector.
