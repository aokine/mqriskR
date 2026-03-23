## Listing 13.1: Multiple-decrement table construction
library(mqriskR)

x <- c(45:50)

qmat <- cbind(q1 = c(.011, .012, .013, .014, .015, .016),
              q2 = c(.100, .100, .100, .100, .100, .100))

tbl <- md_table(x = x, qxj = qmat, radix = 1000)
tbl

## Listing 13.2: Multi-year probabilities from a multiple-decrement table
# 3-year survival probability from age 46
npxtau_md(tbl, x = 46, n = 3)

# 2-year decrement due to Cause 2 from age 47
nqxj_md(tbl, x = 47, n = 2, j = 2)

# 2-year decrement due to Cause 1 from age 46
nqxj_md(tbl, x = 46, n = 2, j = 1)



## Listing 13.3: Constant-force multiple-decrement calculations
mu <- c(0.10, 0.20)
t <- 5

# total survival probability
tpx_tau_cf(mu, t)

# single-decrement probability
tqxprimej_cf(mu = 0.10, t = t)

# dependent probability for Cause 1
tqxj_cf(mu, j = 1, t = t)

# ultimate probability as t -> infinity
tqxj_cf(mu, j = 1, t = 100)


## Listing 13.4: MUDD conversion of probabilities
qxj <- c(0.20, 0.10)

qxprime_mudd(qxj)


## Listing 13.5: SUDD conversion of probabilities
qxprime_sudd(q1 = 0.20, q2 = 0.10)



## Listing 13.5: Ratio property under constant forces
mu <- c(0.10, 0.20)
t <- 2

lhs <- tqxj_cf(mu, j = 1, t = t) / (1 - tpx_tau_cf(mu, t))
rhs <- mu[1] / sum(mu)

c(lhs = lhs, rhs = rhs)
