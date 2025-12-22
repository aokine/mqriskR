library(mqriskR)

# Save figures here (adjust if your LaTeX expects a different path)
out_dir <- "Graphics"
if (!dir.exists(out_dir)) dir.create(out_dir)

# -------------------------
# Figure 1.1: accumulation function a(t) = (1+i)^t
# -------------------------
i <- 0.05
t <- 0:30
a_t <- (1 + i)^t

png(filename = file.path(out_dir, "fig1.1.png"),
    width = 1200, height = 800, res = 150)
plot(t, a_t,
     xlab = "t", ylab = "a(t)",
     main = expression(paste("Accumulation function  ", a(t) == (1+i)^t)))
grid()
dev.off()

# -------------------------
# Figure 1.npv: NPV curve for equation of value
# -------------------------
X <- 100
Y <- 120
cf <- c(-X, -X, -X, -X, rep(Y, 5))
tt <- c(0, 1, 2, 3, 4:8)

i_hat <- solve_yield(cf, tt, interval = c(-0.5, 1))

grid_i <- seq(-0.2, 0.6, by = 0.002)
vals <- vapply(grid_i, function(ii) pv_cashflows(cf, tt, ii), numeric(1))

png(filename = file.path(out_dir, "fig1.10.png"),
    width = 1200, height = 800, res = 150)
plot(grid_i, vals,
     xlab = "i", ylab = "NPV(i)",
     main = "Equation of Value: NPV(i) and Yield Rate")
abline(h = 0)
abline(v = i_hat, lty = 2)
grid()
dev.off()


