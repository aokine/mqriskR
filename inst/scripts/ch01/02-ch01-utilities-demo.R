library(mqriskR)

# ----------------------------
# 1) Interest measure conversions
# ----------------------------
out_i <- interest_convert(i = 0.05, m = 12)
out_i

out_d <- interest_convert(d = 0.04)
out_d

out_delta <- interest_convert(delta = log(1.05))
out_delta

# ----------------------------
# 2) Discount factors / accumulation
# ----------------------------
i <- 0.05
t <- 0:30
a_t <- (1 + i)^t
v_t <- discount(i, t)

head(data.frame(t = t, a_t = a_t, v_t = v_t))

# Plot accumulation function a(t) = (1+i)^t for Fig 1.1
plot(t, a_t,
     xlab = "t", ylab = "a(t)",
     main = expression(paste("Accumulation function  ", a(t) == (1+i)^t)))
grid()

# ----------------------------
# 3) Level annuities: annual, m-thly, continuous
# ----------------------------
n <- 10
a_immediate <- annuity_certain(n, i, due = FALSE)
a_due       <- annuity_certain(n, i, due = TRUE)
a_cont      <- annuity_certain(n, i, cont = TRUE)

c(a_immediate = a_immediate, a_due = a_due, a_cont = a_cont)

# Verify identity: a_due = (1+i)*a_immediate
(1 + i) * a_immediate
a_due

# m-thly (e.g., monthly payments of 1/12 for n years)
a_monthly <- annuity_certain(n = 10, i = 0.05, m = 12, due = FALSE)
a_monthly

# ----------------------------
# 4) Non-level payments (vector approach)
# ----------------------------
# Arithmetic increasing annuity-immediate: 1,2,...,n paid at t=1,...,n
n <- 10
pmt_arith_inc <- 1:n
t_arith_inc <- 1:n
pv_arith_inc <- pv_cashflows(cf = pmt_arith_inc, t = t_arith_inc, i = i)
pv_arith_inc

# Geometric annuity-immediate: 1, r, r^2, ..., r^(n-1) at t=1,...,n
r <- 1.03
pmt_geom <- r^(0:(n-1))
t_geom <- 1:n
pv_geom <- pv_cashflows(cf = pmt_geom, t = t_geom, i = i)
pv_geom

# Plot payment patterns
plot(1:n, pmt_arith_inc, type = "b",
     xlab = "Payment number", ylab = "Payment",
     main = "Non-level payment patterns")
lines(1:n, pmt_geom, type = "b")
legend("topleft", legend = c("Arithmetic increasing", "Geometric (r=1.03)"),
       lty = 1, pch = 1, bty = "n")

# ----------------------------
# 5) Equation of value: verify yield solves NPV=0
# ----------------------------
cf <- c(-100, -100, -100, -100, rep(120, 5))
t  <- c(0, 1, 2, 3, 4:8)
i_hat <- solve_yield(cf, t, interval = c(-0.5, 1))
i_hat
pv_cashflows(cf, t, i_hat)  # should be ~0
