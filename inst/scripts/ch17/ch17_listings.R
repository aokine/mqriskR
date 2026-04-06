##Listing 17.1: Profit vector calculation
library(mqriskR)

V <- c(0, 5.66, 6.17, 0)
qx <- c(0.00142, 0.00153, 0.00166)

Pr <- Pr_vector_disc(V = V, G = 95, i = 0.06, r = 0.05, e = 10, q1 = qx, b1 = 50000,pre_contract_expense = 15)

round(Pr, 2)



##Listing 17.2: Profit signature
p_tau <- c(0.99858, 0.99847, 0.99834)

Pi <- Pi_signature(Pr, p_tau = p_tau)

round(Pi, 2)



##Listing 17.3: Net present value
NPV_profit(Pi, r = 0.10)


##Listing 17.4: Internal rate of return
IRR_profit(Pi, interval = c(0, 1))



##Listing 17.5: Profit margin
APV_GP <- APV_gross_premiums(G = rep(95, 3), r = 0.10, p_tau = p_tau)

NPV_val <- NPV_profit(Pi, r = 0.10)

profit_margin(NPV_val, APV_GP)



##Listing 17.6: Partial NPV and payback period
partial <- NPV_partial(Pi, r = 0.10)

payback <- discounted_payback_period(Pi, r = 0.10)

round(partial, 2)
payback


##Listing 17.7: Zeroized reserves
Vz <- V_zeroized(qx = c(0.015, 0.017, 0.019, 0.021, 0.024), i = 0.06, G = 19279, benefit = 1000000, e = 240)

round(Vz, 2)




