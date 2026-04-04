##Listing 16.1: Type B universal life account value path
library(mqriskR)

B <- 100000; G <- rep(5000, 5)
r <- c(0.75, rep(0.10, 4)); e <- c(100, rep(20, 4))
qx <- c(0.00076, 0.00081, 0.00085, 0.00090, 0.00095)
iq <- 0.03; ic <- 0.03

AV <- numeric(5); COI <- numeric(5)

for (t in 1:5) {prior <- if (t == 1) 0 else AV[t - 1]
COI[t] <- B * qx[t] / (1 + iq)
AV[t] <- (prior + G[t] * (1 - r[t]) - e[t] - COI[t]) * (1 + ic)}

direct_out <- data.frame(year = 1:5,COI = round(COI, 2),AV  = round(AV, 2))

pkg_out <- AV_path_ul_typeB(G = G, r = r, e = e, qx = qx,ic = ic, B = B, iq = iq)

direct_out
pkg_out[, c("t", "COI", "AV")]



##Listing 16.2: Type A universal life account value path
library(mqriskR)

B <- 100000; G <- rep(5000, 5)
r <- c(0.75, rep(0.10, 4)); e <- c(100, rep(20, 4))
qx <- c(0.00076, 0.00081, 0.00085, 0.00090, 0.00095)
i <- 0.03

AV <- numeric(5)

for (t in 1:5) {prior <- if (t == 1) 0 else AV[t - 1]
AV[t] <- ((prior + G[t] * (1 - r[t]) - e[t]) * (1 + i) - B * qx[t]) / (1 - qx[t])}

direct_out <- data.frame(year = 1:5,AV = round(AV, 2))

pkg_out <- AV_path_ul_typeA(G = G, r = r, e = e, qx = qx,ic = i, B = B)

direct_out
pkg_out[, c("t", "AV")]




##Listing 16.3: Point-to-point credited rates under EIUL
library(mqriskR)

index <- c(1000, 1050, 1200, 1100, 950, 1060, 1150)
cap <- 0.10
floor <- 0.01
part <- 1.10

raw <- index[-1] / index[-length(index)] - 1
after_participation <- part * raw
credited <- pmin(cap, pmax(floor, after_participation))

data.frame(year = 1:6,
  raw_growth = round(raw, 4),
  after_participation = round(after_participation, 4),
  credited_rate = round(credited, 4))



##Listing 16.4: Monthly-average credited rate under EIUL
library(mqriskR)

index <- c(1000, 1020, 1100, 1150, 1080, 1040, 960,
           1030, 1000, 1070, 1150, 1200, 1150)

raw <- iMA_eiul(index)

part <- 1.10
cap <- 0.10
floor <- 0.01

credited <- pmin(cap, pmax(floor, part * raw))

c(raw_growth = raw, credited_rate = credited)




##Listing 16.5: Persistency and cumulative in-force probabilities
library(mqriskR)

q_d <- c(.001, .002, .003, .004, .005)
q_w <- c(.02, .02, .03, .04, .05)

# Direct computation from formulas
p_tau_direct <- (1 - q_d) * (1 - q_w)
tp_tau_direct <- cumprod(p_tau_direct)

# Package functions
p_tau_pkg <- pxtau_ul(qd = q_d, qw = q_w)
tp_tau_pkg <- tpxtau_ul(qd = q_d, qw = q_w)

data.frame(year = 1:5,
           p_tau_direct = round(p_tau_direct, 5),
           p_tau_pkg    = round(p_tau_pkg, 5),
           tp_tau_direct = round(tp_tau_direct, 5),
           tp_tau_pkg    = round(tp_tau_pkg, 5))



