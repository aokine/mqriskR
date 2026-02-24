library(mqriskR)


# -------------------------
# Figure 1.1: accumulation function a(t) = (1+i)^t
# -------------------------
i <- 0.05
t <- 0:30
a_t <- (1 + i)^t


plot(t, a_t,
     xlab = "t", ylab = "a(t)",
     main = expression(paste("Accumulation function  ", a(t) == (1+i)^t)))




# ----------------------------
# Level annuities:  m-thly,
# ----------------------------

# m-thly (e.g., monthly payments of 1/12 for n years)
a_monthly <- annuity_certain(n = 10, i = 0.05, m = 12, due = FALSE)
a_monthly


# ----------------------------
# Non-level payments (vector approach)
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
