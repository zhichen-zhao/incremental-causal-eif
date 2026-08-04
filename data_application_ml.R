# =============================================================================
# Data application: HPV testing and CIN2+ detection (Section 6)
#
# Paper: "Efficient Inference for Incremental Causal Effects of Time to
#         Treatment" (Zhao, Ying & Xu).
#
# Applies the cross-fitted, ML-nuisance EIF/AIPW estimator to the Norwegian
# cervical-cancer-screening example. For the PreTectProofer group we estimate
# psi(theta) over a grid of hazard ratios theta, scaling the subsequent-testing
# hazard, and report pointwise 95% CIs together with a uniform 95% confidence
# band (multiplier bootstrap). The estimator mirrors simulation_ml.R
# (RF for mu, spline hazard for Lambda, K-fold cross-fitting).
#
# DATA. The original screening data are not publicly available. As described in
# Section 6, we use the SIMULATED data set from the companion repository of
# Roysland et al. (2025):
#     https://github.com/palryalen/paper-code
# Download their simulated data and save it as "sim_data.RData" in this
# working directory before running. It has 1428 individuals (1147 Amplicor/HC2,
# 281 PreTectProofer); "technology" = DNA (Amplicor/HC2) or RNA (PreTectProofer),
# and the outcome Y indicates CIN2+ detection.
# =============================================================================

#-----------------------------------------------------
library(dplyr)
library(survival)
library(ggplot2)
library(polspline)
library(ranger)
library(flexsurv)   # flexsurvspline: spline hazard for Lambda
#-----------------------------------------------------

#-----------------------------------------------------
# Simulated screening data from Roysland et al. (2025); see header note above.
load("sim_data.RData")
#-----------------------------------------------------

#-----------------------------------------------------
id_list <- unique(sim_data$id)
tau <- 4
df <- data.frame(id = id_list,
                 U = rep(NA, length(id_list)), 
                 Delta = rep(NA, length(id_list)), 
                 Y = rep(NA, length(id_list)), 
                 technology = rep(NA, length(id_list)))

for (i in id_list) {
  data <- sim_data %>% filter(id == i)
  time_vec <- matrix(as.matrix(data[, c("from", "to")]), nrow = 1, byrow = T)
  state_vec <- matrix(as.matrix(data[, c("from.state", "to.state")]), nrow = 1, byrow = T)
  T.time <- ifelse(length(unique(time_vec[state_vec == "follow-up"])) == 0, 
                   tau, unique(time_vec[state_vec == "follow-up"]))
  df[i, "U"] <- min(T.time, tau)
  df[i, "Delta"] <- as.numeric(T.time < tau)
  df[i, "Y"] <- sum(state_vec == "cin2+")
  df[i, "technology"] <- unique(data$technology)
  rm(data); rm(time_vec); rm(state_vec); rm(T.time)
}

df <- df[df$U!=0, ]
#-----------------------------------------------------

#-----------------------------------------------------
dfDNA <- df[df$technology == "DNA", ]
dfRNA <- df[df$technology == "RNA", ]
#-----------------------------------------------------

#-----------------------------------------------------
# -- phi_i(theta) for the first and second integrals
phi <- function(i, theta, data, mu_model, tau, Lambda, d_neg_S, u_int, surv_fit) {
  U_i <- data$U[i]
  Y_i <- data$Y[i]
  Delta_i <- data$Delta[i]
  Lambda_i <- summary(surv_fit, t = U_i, type = "cumhaz", tidy = TRUE)[, "est"]
  
  Lambda_vals <- Lambda
  d_neg_S_vals <- d_neg_S
  mu_vals <- predict(mu_model, data = data.frame(U = u_int))$predictions[, 2]
  
  term1 <- Y_i * theta^Delta_i * exp(-(theta - 1) * Lambda_i)
  term2 <- sum(mu_vals * theta^(u_int < tau) * exp(-(theta - 2) * Lambda_vals) * (u_int <= U_i) * d_neg_S_vals)
  term3 <- sum(mu_vals * theta^(u_int < tau) * exp(-(theta - 1) * Lambda_vals) * d_neg_S_vals)
  return(term1 + (theta - 1) * term2 - (theta - 1) * term3)
}

# -- give phi_i for each observation
cross_fitted_phi <- function(theta, data, tau, K = 5, global_seed = 12345) {
  n <- nrow(data)
  set.seed(global_seed)
  folds <- sample(rep(1:K, length.out = n))  # Randomly assign group IDs
  phi_i_v <- numeric(n)
  
  for (k in 1:K) {
    idx_train <- which(folds != k)
    idx_valid <- which(folds == k)
    train_data <- data[idx_train, ]
    valid_data <- data[idx_valid, ]

    # Fit RF for mu
    train_data$Y <- as.factor(train_data$Y)
    mu_model <- ranger(Y ~ U, probability = TRUE,
                       data = train_data, 
                       num.trees = 500, 
                       min.node.size = 50, 
                       mtry = 1,
                       seed = global_seed + 100 * k)
    
    # Fit spline for Lambda
    surv_fit <- flexsurvspline(Surv(U, Delta) ~ 1, data = train_data, k = 1)
    
    # Predict Lambda for validation set
    u_int <- sort(unique(train_data$U[train_data$Delta == 1]))
    Lambda <- summary(surv_fit, t = u_int, type = "cumhaz", tidy = TRUE)[, "est"]
    Lambda <- c(Lambda, Lambda[length(Lambda)])
    neg_S <- -exp(-Lambda)
    d_neg_S <- diff(c(-1, neg_S))
    d_neg_S[length(d_neg_S)] <- -neg_S[length(neg_S)]
    
    u_int <- c(u_int, tau)
    
    # Compute phi for each obs in the validation fold (local index j maps to global idx_valid[j])
    phi_vals <- sapply(1:length(idx_valid), phi, theta = theta, data = valid_data, 
                       mu_model = mu_model, tau = tau, Lambda = Lambda, d_neg_S = d_neg_S, u_int = u_int, surv_fit = surv_fit)
    
    phi_i_v[idx_valid] <- phi_vals
  }
  
  return(phi_i_v)
}
#------------------------------------------------

#------------------------------------------------
# Uniform confidence band via Rademacher multiplier bootstrap (Theorem 4).
# Approximates the (1-alpha) quantile of sup_theta | multiplier process | to
# get the critical value c_alpha for a band valid uniformly over theta.
ucb <- function(est_mat, est, sd, B, alpha = 0.05) {
  n <- nrow(est_mat)
  sup_vals <- replicate(B, {
    xi <- sample(c(-1, 1), n, replace = TRUE)
    Z <- scale(est_mat, center = est, scale = sd)
    process_b <- (t(xi) %*% Z) / sqrt(n)
    max(abs(process_b))
  })
  c_alpha <- quantile(sup_vals, 1 - alpha)
  
  return(c_alpha)
}
#-----------------------------------------------------

#-----------------------------------------------------
tau <- 4
theta <- seq(0.5, 2, 0.01)
est_mat <- sapply(theta, cross_fitted_phi, data = dfRNA, tau = tau)
est <- colMeans(est_mat)

N <- nrow(dfRNA)
sd <- apply(est_mat, 2, sd)
se <- sd / sqrt(N)
CI_L <- est + qnorm(0.025) * se; CI_U <- est + qnorm(0.975) * se

B <- 10000
c_alpha <- ucb(est_mat, est, sd, B, alpha = 0.05)
CB_L <- est - c_alpha * se; CB_U <- est + c_alpha * se

result <- list(est_set = est, se_set = se, CI_L = CI_L, CI_U = CI_U, 
               c_alpha_set = c_alpha, CB_L = CB_L, CB_U = CB_U)
#-----------------------------------------------------



