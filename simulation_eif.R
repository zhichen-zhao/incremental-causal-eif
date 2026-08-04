#!/usr/bin/env Rscript

# =============================================================================
# Simulation: EIF/AIPW estimator with (semi)parametric nuisances (NO cross-fitting)
#
# Paper: "Efficient Inference for Incremental Causal Effects of Time to
#         Treatment" (Zhao, Ying & Xu), Section 5.
#
# This script evaluates the AIPW estimator psi_hat(theta) built from the
# efficient influence function (EIF), where the two nuisance functions are fit
# with simple (semi)parametric models:
#   - Lambda (treatment-initiation cumulative hazard): Cox proportional hazards
#   - mu (outcome regression E[Y | U, L]):             linear model
# No cross-fitting is used here (that is the companion script simulation_ml.R).
#
# The EIF has three terms (term1 - term2 + term3), matching phi in Theorem 1 /
# Corollary 1 of the paper. psi_hat(theta) = mean over i of phi_i(theta).
#
# USAGE (run from the command line; N = sample size, seed = RNG seed):
#   Rscript simulation_eif.R <N> <seed>
# e.g.
#   Rscript simulation_eif.R 1000 1
# The per-seed result (point estimate, SE, coverage indicator) is written to
#   eif_<N>/seed_<seed>.rds
# Reproducing the paper's tables means running this across many seeds and
# sample sizes (R = 1000 replications, n in {200, 1000, 5000}) and averaging.
# =============================================================================

#------------------------------------------------
library(survival)   # coxph, basehaz for the Cox nuisance (Lambda)
library(flexsurv)
#------------------------------------------------

#------------------------------------------------
# ---- Data-generating parameters (see Section 5 of the paper) ----------------
# L ~ Unif(0,2); treatment time T has hazard exp(beta*L); outcome Y is linear in
# L and the (censored) treatment time U = T ^ tau, with Gaussian noise.
para_set <- list(beta = 0.2, 
                 mu_cY = 3, mu_LY = 0.6, mu_UY = 1, sigma_Y = 0.5)

# Draw T given L under an intervened hazard theta(t,l) = (a t + b) exp(beta l),
# by inverting the cumulative hazard (closed form for this linear-in-t shape).
rT_given_L <- function(L, a, b, beta) {
  n <- length(L)
  eta <- as.vector(L * (para_set$beta + beta))
  U <- runif(n)
  c <- (-log(U)) / exp(eta)
  
  time <- (-b + sqrt(b^2 + 2 * a * c)) / a
  
  return(time)
}

data_gen <- function(N, para_set, hz_ratio) {
  L <- runif(N, 0, 2)
  
  if (all(hz_ratio == 0)) {
    U <- rexp(N, exp(L * para_set$beta))
  }
  else {
    theta_a <- hz_ratio[1]; theta_b <- hz_ratio[2]; theta_beta <- hz_ratio[3]
    U <- rT_given_L(L, theta_a, theta_b, theta_beta)}
  
  Delta <- U < 2
  U <- pmin(U, 2)
  Y <- para_set$mu_cY - L * para_set$mu_LY - para_set$mu_UY * (2 - U) + rnorm(N, 0, para_set$sigma_Y)
  
  df <- data.frame(L, U, Y, Delta)
  
  return(df)
}

# ---- Ground-truth psi(theta) via a huge Monte Carlo sample ------------------
# Each element of cand_set is a triple (a, b, beta) defining an intervention
# theta(t,l) = (a t + b) exp(beta l). These are the eight scenarios theta_1..theta_8.
set.seed(12345)
N <- 10^7
cand_set <- 
  list(c(0.9, 0.3, -0.7), c(0.9, 0.5, -0.7), c(0.7, 0.3, -0.5), c(0.7, 0.5, -0.5), 
       c(0.5, 0.1, -0.1), c(0.5, 0.1, -0.2), c(0.3, 0.1,  0.4), c(0.3, 0.1,  0.6))
truth_set <- c()
p <- c()

for (hz_ratio in cand_set) {
  df <- data_gen(N, para_set, hz_ratio)
  truth_set <- c(truth_set, mean(df$Y))  # true psi(theta) for this scenario
  p <- c(p, mean(df$Delta))              # event (non-censoring) rate
}
#------------------------------------------------

#------------------------------------------------
# -- Lambda_hat
# return Lambda, dLambda, and dexp(Lambda) at all event time given L=l
Lambda_hat_cox <- function(l, fit, bh) {
  Lambda <- bh$hazard * exp(as.numeric(l * fit$coefficients))
  d_Lambda <- diff(c(0, Lambda))
  d_e_Lambda <- diff(c(1, exp(Lambda)))

  return(list(Lambda = Lambda,
              d_Lambda = d_Lambda,
              d_e_Lambda = d_e_Lambda))
}

# -- mu_hat(u, l)
mu_hat <- function(u, l, mu_model) {
  est <- predict(mu_model, newdata = data.frame(U = u, L = l))
  
  return(est)
}

# -- theta_fun(t, l)
theta_fun <- function(t, l, hz_ratio) {
  theta_a <- hz_ratio[1]
  theta_b <- hz_ratio[2]
  theta_beta <- hz_ratio[3]
  
  val <- (theta_a * t + theta_b) * exp((l * theta_beta))
  
  return(val)
}

# -- phi_i(theta) for a single individual i
phi_i <- function(i, hz_ratio, data, mu_model, tau, fit, bh, u_int) {
  L_i <- as.numeric(data[i, "L"])
  U_i <- data$U[i]
  Delta_i <- data$Delta[i]
  Y_i     <- data$Y[i]
  
  theta_i <- theta_fun(U_i, L_i, hz_ratio)
  
  Lam_dLam_deLam_vals <- Lambda_hat_cox(l = L_i, fit = fit, bh = bh)
  Lam_v <- Lam_dLam_deLam_vals$Lambda
  dLam_v <- Lam_dLam_deLam_vals$d_Lambda
  deLam_v <- Lam_dLam_deLam_vals$d_e_Lambda

  theta_v <- theta_fun(u_int, L_i, hz_ratio)
  mu_v <- mu_hat(u_int, L_i, mu_model)
  dfLam_v <- diff(c(1, exp(-cumsum(theta_v * dLam_v))))
  dfLam_v[length(dfLam_v)] <- -exp(-sum(theta_v * dLam_v))
  
  # --- term 1
  int_t1_val <- sum((u_int <= U_i) * (theta_v - 1) * dLam_v)
  term1 <- Y_i * (theta_i^Delta_i) * exp(-int_t1_val)

  # --- term 2
  inner_v <- cumsum((theta_v - 1) * deLam_v)
  int_t2_val <- inner_v
  idx_ui <- max(which(u_int <= U_i))
  if(idx_ui < length(u_int)) {
    int_t2_val[(idx_ui + 1):length(u_int)] <- inner_v[idx_ui]
  }
  term2 <- sum(mu_v * int_t2_val * dfLam_v)

  # --- term 3
  denom_t3 <- exp(-Lam_v[idx_ui])
  int_t3_val <- sum((u_int > U_i) * mu_v * dfLam_v)
  term3 <- ((theta_i - 1) / denom_t3) * int_t3_val
  
  return(term1 - term2 + term3)
}

# -- psi_hat(theta): fit the two (semi)parametric nuisances, then average phi_i
#    Lambda: Cox proportional hazards;  mu: linear regression.  No cross-fitting.
psi_hat <- function(hz_ratio, data, tau) {
  n <- nrow(data)
  fit <- coxph(Surv(time = U, event = Delta) ~ L, data = data)  # nuisance Lambda
  mu_model <- lm(Y ~ U + L, data = data)                        # nuisance mu
  bh <- basehaz(fit, centered = FALSE)
  u_int <- bh$time
  
  phi_i_v <- sapply(1:n, phi_i, hz_ratio = hz_ratio, data = data, mu_model = mu_model, 
                    tau = tau, fit = fit, bh = bh, u_int = u_int)
  
  return(phi_i_v)
}
#------------------------------------------------

#------------------------------------------------
# Bootstrap estimates for psi(theta)
psi_boot_W <- function(hz_ratio, data, tau, B.seed) {
  set.seed(B.seed)
  boot_weights <- rexp(N)
  boot_weights <- boot_weights / mean(boot_weights)
  n <- nrow(data)
  fit <- coxph(Surv(time = U, event = Delta) ~ L, data = data, weights = boot_weights)
  mu_model <- lm(Y ~ U + L, data = data, weights = boot_weights)
  bh <- basehaz(fit, centered = FALSE)
  u_int <- bh$time
  
  boot_est <- mean(boot_weights * 
                     sapply(1:n, phi_i, hz_ratio = hz_ratio, data = data, mu_model = mu_model,
                            tau = tau, fit = fit, bh = bh, u_int = u_int))
  
  return(boot_est)
}

psi_ci_W <- function(hz_ratio, data, tau, B) {
  psi_boots_W <- sapply(1:B, psi_boot_W, hz_ratio = hz_ratio, data = data, tau = tau)
  
  return(sd(psi_boots_W))
}
#------------------------------------------------

#------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
N    <- as.integer(args[1])
seed <- as.integer(args[2])

#B <- 200
tau <- 2

set.seed(seed)
df <- data_gen(N, para_set, c(0, 0, 0))
est_mat <- sapply(cand_set, psi_hat, data = df, tau = tau)
est <- colMeans(est_mat)

se <- apply(est_mat, 2, sd) / sqrt(N)
CP <- truth_set <= est + qnorm(0.975) * se & truth_set >= est + qnorm(0.025) * se

#boot_se <- sapply(cand_set, psi_ci_W, data = df, tau = tau, B = B)
#boot_CP <- truth_set <= est + qnorm(0.975) * boot_se & truth_set >= est + qnorm(0.025) * boot_se

#result <- list(est_set = est, se_set = se, CP_set = CP, boot_se_set = boot_se, boot_CP_set = boot_CP)

result <- list(est_set = est, se_set = se, CP_set = CP)
#------------------------------------------------

#------------------------------------------------
dir_name <- paste0("eif_", N)

dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)

saveRDS(result, file = file.path(getwd(), dir_name, paste0("seed_", seed, ".rds")))
#------------------------------------------------

