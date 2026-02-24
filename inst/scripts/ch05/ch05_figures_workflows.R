## Figure 5.1

library(mqriskR)

# Time grid (avoid the uniform model's endpoint to prevent division-by-zero in the hazard)
t_grid <- seq(0, 100, length.out = 501)
t_grid_u <- t_grid[t_grid < 100]   # for uniform hazard where omega = 100

# --- Choose parameters (adjust as desired for your narrative) ---
omega <- 100

lambda <- 0.03          # exponential
B_g <- 0.0005; c_g <- 1.08        # Gompertz
A_m <- 0.0002; B_m <- 0.0004; c_m <- 1.08   # Makeham

# --- Survival functions ---
S_u <- S0(t_grid, model = "uniform", omega = omega)
S_e <- S0(t_grid, model = "exponential", lambda = lambda)
S_g <- S0(t_grid, model = "gompertz", B = B_g, c = c_g)
S_m <- S0(t_grid, model = "makeham", A = A_m, B = B_m, c = c_m)

# --- Hazard functions ---
h_u <- hazard0(t_grid_u, model = "uniform", omega = omega)
h_e <- hazard0(t_grid,   model = "exponential", lambda = lambda)
h_g <- hazard0(t_grid,   model = "gompertz", B = B_g, c = c_g)
h_m <- hazard0(t_grid,   model = "makeham", A = A_m, B = B_m, c = c_m)

# Side-by-side plots
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# ---- Panel 1: survival functions ----
plot(t_grid, S_u, type = "l", lty = 1,
     xlab = "t (time / age)", ylab = expression(S[0](t) == Pr(T[0] > t)),
     main = expression(paste("Survival functions: ", S[0](t))), ylim = c(0, 1))
lines(t_grid, S_e, lty = 2)
lines(t_grid, S_g, lty = 3)
lines(t_grid, S_m, lty = 4)
#legend("topright",
#       legend = c("Uniform (de Moivre)", "Exponential", "Gompertz", "Makeham"),
#       lty = 1:4, bty = "n", cex = 0.85)

# ---- Panel 2: hazard functions ----
plot(t_grid_u, h_u, type = "l", lty = 1,
     xlab = "t (time / age)", ylab = expression(lambda[0](t)),
     main = expression("Hazard Functions: " * lambda[0](t)))
lines(t_grid, h_e, lty = 2)
lines(t_grid, h_g, lty = 3)
lines(t_grid, h_m, lty = 4)
legend("topleft",
       legend = c("Uniform (de Moivre)", "Exponential", "Gompertz", "Makeham"),
       lty = 1:4, bty = "n", cex = 0.85)

par(op)


