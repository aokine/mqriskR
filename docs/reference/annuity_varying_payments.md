# Varying-payment annuity functions (Chapter 8)

Chapter 8 non-level annuity functions for increasing and decreasing life
annuities.

Computes \$\$(Ia)\_x = \sum\_{t=1}^{\infty} t \\ v^t \\ {}\_tp_x.\$\$

Computes \$\$(Ia)\_{x:\overline{n}\|} = \sum\_{t=1}^{n} t \\ v^t \\
{}\_tp_x.\$\$

Computes \$\$(Da)\_{x:\overline{n}\|} = \sum\_{t=1}^{n} (n+1-t)\\ v^t \\
{}\_tp_x.\$\$

Computes \$\$(I\ddot{a})\_x = \sum\_{t=0}^{\infty} (t+1)\\ v^t \\
{}\_tp_x.\$\$

Computes \$\$(I\ddot{a})\_{x:\overline{n}\|} = \sum\_{t=0}^{n-1} (t+1)\\
v^t \\ {}\_tp_x.\$\$

Computes \$\$(D\ddot{a})\_{x:\overline{n}\|} = \sum\_{t=0}^{n-1} (n-t)\\
v^t \\ {}\_tp_x.\$\$

Computes \$\$(\bar{I}\bar{a})\_x = \int_0^\infty
t\\v^t\\{}\_tp_x\\dt.\$\$

Computes \$\$(\bar{I}\bar{a})\_{x:\overline{n}\|} = \int_0^n
t\\v^t\\{}\_tp_x\\dt.\$\$

Computes \$\$(\bar{D}\bar{a})\_{x:\overline{n}\|} = \int_0^n
(n-t)\\v^t\\{}\_tp_x\\dt.\$\$

## Usage

``` r
Iax(x, i, model, ..., k_max = 5000, tol = 1e-12)

Iaxn(x, n, i, model, ...)

Daxn(x, n, i, model, ...)

Iadotx(x, i, model, ..., k_max = 5000, tol = 1e-12)

Iadotxn(x, n, i, model, ...)

Dadotxn(x, n, i, model, ...)

Iabarx(x, i, model, ..., tol = 1e-10)

Iabarxn(x, n, i, model, ...)

Dabarxn(x, n, i, model, ...)
```

## Arguments

- x:

  Age.

- i:

  Effective annual interest rate.

- model:

  Survival model.

- ...:

  Additional model parameters.

- k_max:

  Maximum summation horizon for non-terminating models.

- tol:

  Truncation tolerance for non-terminating models.

- n:

  Term in years.

## Value

Numeric vector.

Numeric vector.

Numeric vector.

## Details

The functions implemented here match the notation in Section 8.6:

- `Iax()` = \\(Ia)\_x\\

- `Iaxn()` = \\(Ia)\_{x:\overline{n}\|}\\

- `Daxn()` = \\(Da)\_{x:\overline{n}\|}\\

- `Iadotx()` = \\(I\ddot{a})\_x\\

- `Iadotxn()` = \\(I\ddot{a})\_{x:\overline{n}\|}\\

- `Dadotxn()` = \\(D\ddot{a})\_{x:\overline{n}\|}\\

- `Iabarx()` = \\(\bar{I}\bar{a})\_x\\

- `Iabarxn()` = \\(\bar{I}\bar{a})\_{x:\overline{n}\|}\\

- `Dabarxn()` = \\(\bar{D}\bar{a})\_{x:\overline{n}\|}\\
