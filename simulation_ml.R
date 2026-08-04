#!/usr/bin/env Rscript

# =============================================================================
# Simulation: EIF/AIPW estimator with MACHINE-LEARNING nuisances + CROSS-FITTING
#
# Paper: "Efficient Inference for Incremental Causal Effects of Time to
#         Treatment" (Zhao, Ying & Xu), Section 5.
#
# Same estimand and same EIF (three terms term1 - term2 + term3) as
# simulation_eif.R, but here the nuisances are estimated with flexible ML and
# combined via K-fold cross-fitting (Section 4 of the paper):
#   - Lambda (cumulative hazard): spline hazard regression, hare() from polspline
#   - mu (outcome regression):    random forest, ranger()
# Cross-fitting (K = 10): for each fold k, nuisances are fit on the other K-1
# folds and phi is evaluated on fold k, so no observation is scored by a model
# trained on it. This is what allows sqrt(n)-inference with ML nuisances.
#
# USAGE (command line; N = sample size, seed = RNG seed):
#   Rscript simulation_ml.R <N> <seed>
# The per-seed result is written to  ml_<N>/seed_<seed>.rds
# =============================================================================

#------------------------------------------------
library(polspline)  # hare(): spline-based hazard regression for Lambda
library(ranger)     # random forest for the outcome regression mu
#------------------------------------------------

#------------------------------------------------
para_set <- list(beta = 0.2, 
                 mu_cY = 3, mu_LY = 0.6, mu_UY = 1, sigma_Y = 0.5)

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

set.seed(12345)
N <- 10^7
cand_set <- 
  list(c(0.9, 0.3, -0.7), c(0.9, 0.5, -0.7), c(0.7, 0.3, -0.5), c(0.7, 0.5, -0.5), 
       c(0.5, 0.1, -0.1), c(0.5, 0.1, -0.2), c(0.3, 0.1,  0.4), c(0.3, 0.1,  0.6))
truth_set <- c()

for (hz_ratio in cand_set) {
  df <- data_gen(N, para_set, hz_ratio)
  truth_set <- c(truth_set, mean(df$Y))
}
#------------------------------------------------

#------------------------------------------------
# -- theta_fun(t, l)
theta_fun <- function(t, l, hz_ratio) {
  theta_a <- hz_ratio[1]
  theta_b <- hz_ratio[2]
  theta_beta <- hz_ratio[3]
  
  val <- (theta_a * t + theta_b) * exp((l * theta_beta))
  
  return(val)
}

# -- phi_i(theta) for the first and second integrals
phi <- function(i, hz_ratio, data, mu_model, tau, Lam_mat, dLam_mat, deLam_mat, u_int) {
  L_i <- data[i, "L"]
  U_i <- data$U[i]
  Delta_i <- data$Delta[i]
  Y_i     <- data$Y[i]
  
  theta_i <- theta_fun(U_i, L_i, hz_ratio)
  
  Lam_v <- Lam_mat[i, ]
  dLam_v <- dLam_mat[i, ]
  deLam_v <- deLam_mat[i, ]
    
  new_data <- data.frame(U = u_int, L = L_i)
  mu_v <- predict(mu_model, data = new_data)$predictions
  
  theta_v <- theta_fun(u_int, L_i, hz_ratio)
  dfLam_v <- diff(c(1, exp(-cumsum(theta_v * dLam_v))))
  dfLam_v[length(dfLam_v)] <- -exp(-sum(theta_v * dLam_v))
  
  # --- term 1
  int_t1_val <- sum((u_int <= U_i) * (theta_v - 1) * dLam_v)
  term1 <- Y_i * (theta_i^Delta_i) * exp(-int_t1_val)
  
  # --- term 2
  inner_v <- cumsum((theta_v - 1) * deLam_v)
  int_t2_val <- inner_v
  
  indices <- which(u_int <= U_i)
  idx_ui <- if(length(indices) > 0) max(indices) else 0
  
  if (idx_ui == 0) {
    int_t2_val <- rep(0, length(u_int))
  } else if (idx_ui < length(u_int)) {
    int_t2_val[(idx_ui + 1):length(u_int)] <- inner_v[idx_ui]
  }
  
  term2 <- sum(mu_v * int_t2_val * dfLam_v)
  
  # --- term 3
  denom_t3 <- if(idx_ui != 0) exp(-Lam_v[idx_ui]) else 1
  int_t3_val <- sum((u_int > U_i) * mu_v * dfLam_v)
  term3 <- ((theta_i - 1) / denom_t3) * int_t3_val
  
  return(term1 - term2 + term3)
}

# -- Cross-fitted EIF: fit ML nuisances out-of-fold, evaluate phi on each fold
#    K-fold cross-fitting (Section 4). For each fold k:
#      train nuisances (RF for mu, HARE spline for Lambda) on folds != k,
#      then compute phi for the held-out fold k. Returns phi_i for every i.
cross_fitted_phi <- function(hz_ratio, data, tau, K = 10, global_seed = 12345) {
  n <- nrow(data)
  set.seed(global_seed)
  folds <- sample(rep(1:K, length.out = n))  # random fold assignment
  phi_i_v <- numeric(n)
  
  for (k in 1:K) {
    idx_train <- which(folds != k)
    idx_valid <- which(folds == k)
    train_data <- data[idx_train, ]
    valid_data <- data[idx_valid, ]
    
    # Nuisance mu: random forest for E[Y | U, L], trained on the other folds
    mu_model <- ranger(Y ~ U + L, 
                       data = train_data, 
                       num.trees = 100, 
                       min.node.size = 15, 
                       mtry = 2, 
                       seed = global_seed + 100 * k)

    # Fit spline for Lambda
    hare_fit <- hare(data = train_data$U,
                     delta = train_data$Delta,
                     cov = as.matrix(train_data[, "L"]), 
                     prophaz = TRUE)
    
    # Predict Lambda for validation set
    L_valid <- as.matrix(valid_data[, "L"])
    u_int <- sort(unique(train_data$U[train_data$Delta == 1]))
    F_mat <- sapply(u_int, function(t) {phare(q = t, cov = L_valid, fit = hare_fit)})
    
    Lam_mat <- -log(1 - F_mat)
    Lam_mat <- cbind(Lam_mat, Lam_mat[, ncol(Lam_mat)])
    dLam_mat <- t(apply(Lam_mat, 1, function(row) {diff(c(0, row))}))
    
    eLam_mat <- 1 / (1 - F_mat)
    eLam_mat <- cbind(eLam_mat, eLam_mat[, ncol(eLam_mat)])
    deLam_mat <- t(apply(eLam_mat, 1, function(row) {diff(c(1, row))}))
    
    u_int <- c(u_int, tau)
    
    # Compute phi for each obs in the validation fold (local index j maps to global idx_valid[j])
    phi_vals <- sapply(1:length(idx_valid), phi, hz_ratio = hz_ratio, data = valid_data, 
                      mu_model = mu_model, tau = tau, Lam_mat = Lam_mat, dLam_mat = dLam_mat, 
                      deLam_mat = deLam_mat, u_int = u_int)

    phi_i_v[idx_valid] <- phi_vals
  }
  
  return(phi_i_v)
}
#------------------------------------------------

#------------------------------------------------
N    <- 200
seed <- 12345
tau  <- 2

set.seed(seed)
df <- data_gen(N, para_set, c(0, 0, 0))
est_mat <- sapply(cand_set, cross_fitted_phi, data = df, tau = tau)
est <- colMeans(est_mat)

se <- apply(est_mat, 2, sd) / sqrt(N)
CP <- truth_set <= est + qnorm(0.975) * se & truth_set >= est + qnorm(0.025) * se

result <- list(est_set = est, se_set = se, CP_set = CP)
#------------------------------------------------


