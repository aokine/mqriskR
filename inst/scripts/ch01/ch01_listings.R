## Listing 1.1: Converting between $i$, $d$, and $\delta$ (compound interest)

library(mqriskR)

# Convert from i to d and delta
interest_convert(i = 0.05)

# Convert from d to i and delta
interest_convert(d = 0.05/1.05)

# Convert from delta to i and d
interest_convert(delta = log(1.05))

# Nominal annual rate convertible m-thly (e.g., m = 12)
interest_convert(i = 0.06, m = 12)



## Listing 1.2: Verifying annuity-certain identities by computation
library(mqriskR)

i <- 0.05
n <- 10

# Level annuity-immediate and annuity-due (annual)
a_immediate <- annuity_certain(n = n, i = i, due = FALSE)
a_due<- annuity_certain(n = n, i = i, due = TRUE)

# Identity: ä_n = (1+i) a_n
c(a_immediate = a_immediate, a_due = a_due, ratio = a_due / a_immediate)

# Continuous annuity-certain using delta = log(1+i)
a_cont <- annuity_certain(n = n, i = i, cont = TRUE)
a_cont





## Listing 1.3: Non-level annuities via a payment vector
library(mqriskR)

i <- 0.05
n <- 10

# Payments p_k at times t_k: PV = sum(p_k * (1+i)^(-t_k))
pmt <- 1:n      # arithmetic increasing payments 1,2,...,n
t   <- 1:n      # paid at t=1,...,n (annuity-immediate)
pv_cashflows(cf = pmt, t = t, i = i)




## Listing 1.4: Equation of value and yield rate via root finding
library(mqriskR)

# Net cash flows (investor perspective):
# Four payments of X at t=0,1,2,3 and five payments of Y at t=4,...,8
X <- 100
Y <- 120

cf <- c(-X, -X, -X, -X, rep(Y, 5))
t  <- c(0, 1, 2, 3, 4:8)

# Yield rate solves PV(i) = 0
i_hat <- solve_yield(cf, t, interval = c(-0.5, 1))
i_hat

# Plot NPV(i)
grid <- seq(-0.2, 0.6, by = 0.002)
vals <- vapply(grid, function(ii) pv_cashflows(cf, t, ii), numeric(1))

plot(grid, vals, xlab = "i", ylab = "NPV(i)",
     main = "Equation of Value: NPV(i) and Yield Rate")
abline(h = 0)
abline(v = i_hat, lty = 2)








