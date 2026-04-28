# Survival models (Chapter 5): core survival functions + parametric models

This file provides actuarial survival-model utilities for Chapter 5.

## Details

\## Age-at-failure random variable (T0) - \`S0(t, ...)\` : survival
function S0(t) = Pr(T0 \> t) - \`F0(t, ...)\` : CDF F0(t) = Pr(T0 \<= t)
= 1 - S0(t) - \`f0(t, ...)\` : density f0(t) = d/dt F0(t) - \`hazard0(t,
...)\` : hazard lambda0(t) = f0(t) / S0(t) - \`cumhaz0(t, ...)\` : cum
hazard Lambda0(t) = \\\int_0^t \lambda_0(y) dy\\

\## Conditional (future lifetime) random variable (Tx), given alive at
age x - \`tpx(t, x, ...)\` : \\{}\_t p_x = S0(x+t)/S0(x)\\ - \`tqx(t, x,
...)\` : \\{}\_t q_x = 1 - {}\_t p_x\\ - \`fx(t, x, ...)\` : \\f_x(t) =
f0(x+t)/S0(x)\\

\## Expectations - \`ex_complete(x, ...)\` : \\\overset{\circ}{e}\_x =
\int_0^\infty {}\_t p_x\\dt\\ - \`ex_curtate(x, ...)\` : \\e_x =
E\[K_x\] = \sum\_{k\ge 1} {}\_k p_x\\ (truncated sum)

\## Models supported - \`"uniform"\` (de Moivre) parameter: \`omega\` -
\`"exponential"\` (constant force) parameter: \`lambda\` -
\`"gompertz"\` parameters: \`B\`, \`c\` - \`"makeham"\` parameters:
\`A\`, \`B\`, \`c\` - \`"weibull"\` parameters: \`shape\`, \`scale\`

All functions are vectorized over \`t\` (and over \`x\` where
applicable).
