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
#
# NOTE. This script also reads intermediate .rds files (eif.result.rds,
# ipw.result.rds, dfDNA.rds, haz_ratio.rds) produced by the companion runs when
# assembling the comparison plot; generate those first, or comment out the
# plotting blocks to run the estimator alone.
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
# data set includes follow-up at 4 which we treat as Delta=0, but Hines treat that case still as Delta=1
#fit <- coxph(Surv(from, to, to.state == "follow-up") ~ technology, data = sim_data[sim_data$from.state == "",])
#summary(fit)
#zph <- cox.zph(fit)
#plot(zph)

#sim_data[sim_data$to.state == "follow-up" & sim_data$to == 4,]
#sim_data[sim_data$from.state == "follow-up" & sim_data$from == 0,]
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
  df[i, "Delta"] <- as.numeric(T.time < tau) # as.numeric(length(unique(time_vec[state_vec == "follow-up"])) != 0)
  df[i, "Y"] <- sum(state_vec == "cin2+")
  df[i, "technology"] <- unique(data$technology)
  rm(data); rm(time_vec); rm(state_vec); rm(T.time)
}

mean(df$Y[df$technology=="DNA"]); mean(df$Y[df$technology=="RNA"])
#df$technology <- as.factor(df$technology)
#df0 <- df
df <- df[df$U!=0, ]
#-----------------------------------------------------

#-----------------------------------------------------
km_fit <- survfit(Surv(time = U, event = Delta) ~ technology, data = df, conf.type = "log-log", conf.int = 0.95)

km_df <- data.frame(time = km_fit$time,
                    surv = 1 - km_fit$surv,
                    lower = 1 - km_fit$lower,
                    upper = 1 - km_fit$upper,
                    technology = rep(names(km_fit$strata), km_fit$strata))

km_df$technology <- gsub("technology=DNA", "Amplicor/HC2", km_df$technology)
km_df$technology <- gsub("technology=RNA", "PreTectProofer", km_df$technology)

# km_df$surv[642]
# tail(km_df$surv)
ggplot(km_df, aes(x = time, y = surv, color = technology, fill = technology)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Time since 2005", y = "Probability of receiving subsequent test") +
  scale_color_manual(values = c("orange", "#1F77B4")) + 
  scale_fill_manual(values = c("orange", "#1F77B4")) + 
  ylim(0, 0.7) +  # Set y-axis limits to 0–0.75
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Add black axis lines
    panel.border = element_blank(),  # Remove default border (if any)
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15)
  )
#-----------------------------------------------------

#-----------------------------------------------------
dfDNA <- df[df$technology == "DNA", ]
dfRNA <- df[df$technology == "RNA", ]

# Coxph fit
df$technology <- relevel(factor(df$technology), ref = "RNA")
cox_fit <- coxph(Surv(U, Delta) ~ technology, data = df)
summary(cox_fit)
#saveRDS(summary(cox_fit)$conf.int[1], "haz_ratio.rds")
#ggsurvplot(survfit(cox_fit), data = df)

cox_zph <- cox.zph(cox_fit)
print(cox_zph) # coxph assumption is violated
#plot(cox_zph)

summary(cox_fit)$conf.int
#saveRDS(dfRNA, "dfRNA.rds");saveRDS(dfDNA, "dfDNA.rds")
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
theta <- sort(c(seq(0.5, 2, 0.01), exp(cox_fit$coefficients)))
est_mat <- sapply(theta, cross_fitted_phi, data = dfRNA, tau = tau)
est <- colMeans(est_mat)

N <- nrow(dfRNA)
sd <- apply(est_mat, 2, sd)
se <- sd / sqrt(N)
CI_L <- est + qnorm(0.025) * se; CI_U <- est + qnorm(0.975) * se

B <- 10000
c_alpha <- ucb(est_mat, est, sd, B, alpha = 0.05)
CB_L <- est - c_alpha * se; CB_U <- est + c_alpha * se

ml_res <- list(est_set = est, se_set = se, CI_L = CI_L, CI_U = CI_U, 
               c_alpha_set = c_alpha, CB_L = CB_L, CB_U = CB_U)

saveRDS(ml_res, "ml.result.rds")
#-----------------------------------------------------

#-----------------------------------------------------
# plot
dfDNA <-  readRDS("dfDNA.rds")
cf_res <- readRDS("ml.result.rds")
eif <- readRDS("eif.result.rds")
ipw <- readRDS("ipw.result.rds")
haz_ratio <- readRDS("haz_ratio.rds")
ybar <- round(eif$est[eif$theta==1], 4)
ybarDNA <- round(mean(dfDNA$Y), 4)

cf <- data.frame(theta = eif$theta, est = cf_res[[1]], 
                 CI_U = cf_res[["CI_U"]], CI_L = cf_res[["CI_L"]], 
                 CB_U = cf_res[["CB_U"]], CB_L = cf_res[["CB_L"]])

library(ggplot2)
library(latex2exp)

label_pw_cf  <- TeX("Pointwise 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_un_cf  <- TeX("Uniform 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_pw_eif <- TeX("Pointwise 95% CI for $\\hat{\\psi}(\\theta)$")
label_est_cf <- TeX("$\\hat{\\psi}_{cf}(\\theta)$")
label_est_eif <- TeX("$\\hat{\\psi}(\\theta)$")

ggplot() +
  # --- Uniform CI ---
  geom_ribbon(data = cf, aes(x = theta, ymin = CB_L, ymax = CB_U, color = "un_cf"), 
              fill = "lightblue", alpha = 0.3, linetype = 0) +
  
  # --- Pointwise CIs ---
  # cf Pointwise
  geom_line(data = cf, aes(x = theta, y = CI_L, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 1) +
  geom_line(data = cf, aes(x = theta, y = CI_U, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 1) +
  # eif Pointwise
  geom_line(data = eif, aes(x = theta, y = CI_L, color = "pw_eif"), 
            linetype = "dotted", alpha = 0.9, linewidth = 1.2) +
  geom_line(data = eif, aes(x = theta, y = CI_U, color = "pw_eif"), 
            linetype = "dotted", alpha = 0.9, linewidth = 1.2) +
  
  # --- Estimates ---
  geom_line(data = cf, aes(x = theta, y = est, color = "est_cf"), 
            linewidth = 2, alpha = 0.7) +
  geom_line(data = eif, aes(x = theta, y = est, color = "est_eif"), 
            linewidth = 0.7, alpha = 0.9) +
  
  # --- Others ---
  annotate("segment", x = 0.5, xend = 2, y = ybar, yend = ybar, 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = 0.5, xend = 2, y = ybarDNA, yend = ybarDNA, 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = haz_ratio, linetype = "dotted", color = "gray70") +
  
  # --- Axis ---
  labs(y = "Estimated proportion of CIN2+ detected", x = expression(paste("Hazards ratio ", theta))) +
  scale_x_log10(breaks = c(0.5, 1, round(haz_ratio, 4), 2, 5),
                labels = c("0.5", "1", round(haz_ratio, 4), "2", "5")) +
  scale_y_continuous(breaks = c(0.000, seq(0.025, 0.125, 0.025), ybar, ybarDNA),
                     labels = c("0.000", "0.025", "0.050", "0.075", "0.100", "0.125",
                                as.character(ybar), as.character(ybarDNA))) +
  
  # --- Label ---
  scale_color_manual(name = NULL,
    values = c("est_eif" = "pink", "est_cf"  = "steelblue", "pw_eif"  = "pink",
               "pw_cf"   = "steelblue", "un_cf"   = "lightblue"),
    labels = c("est_eif" = label_est_eif, "est_cf"  = label_est_cf, "pw_eif"  = label_pw_eif, 
               "pw_cf"   = label_pw_cf, "un_cf"   = label_un_cf),
    breaks = c("est_eif", "est_cf", "pw_eif", "pw_cf", "un_cf")) +

  guides(color = guide_legend(ncol = 1, 
    override.aes = list(
      linetype = c("solid", "solid", "dotted", "longdash", "blank"), 
      fill     = c(NA, NA, NA, NA, "lightblue"), linewidth = c(0.7, 2, 1.2, 1, 0),
      alpha    = c(0.9, 0.7, 0.9, 0.5, 0.3)))) +
  
  theme_classic() +
  theme(legend.position = c(0.02, 1), 
    legend.justification = c("left", "top"), 
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),     
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line = element_line(color = "black", linewidth = 0.3)
)
#-----------------------------------------------------

#-----------------------------------------------------
result <- cf
cand_set <- result$theta

round(100*(result$est[cand_set==2]/result$est[cand_set==1]-1), 1)
round(100*(1-result$est[cand_set==0.5]/result$est[cand_set==1]), 1)

round(result[c(152, 1, 19), ], 4)

#> round(summary(cox_fit)$conf.int, 4)
#exp(coef) exp(-coef) lower .95 upper .95
#technologyDNA    0.6768     1.4776    0.5727    0.7998
#-----------------------------------------------------

#-----------------------------------------------------
# plot
dfDNA <-  readRDS("dfDNA.rds")
cf_res <- readRDS("ml.result.rds")
eif <- readRDS("eif.result.rds")
ipw <- readRDS("ipw.result.rds")
haz_ratio <- readRDS("haz_ratio.rds")
ybar <- round(eif$est[eif$theta==1], 4)
ybarDNA <- round(mean(dfDNA$Y), 4)

cf <- data.frame(theta = eif$theta, est = cf_res[[1]], 
                 CI_U = cf_res[["CI_U"]], CI_L = cf_res[["CI_L"]], 
                 CB_U = cf_res[["CB_U"]], CB_L = cf_res[["CB_L"]])

library(ggplot2)
library(latex2exp)

label_pw_cf  <- TeX("Pointwise 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_un_cf  <- TeX("Uniform 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_pw_eif <- TeX("Pointwise 95% CI for $\\hat{\\psi}(\\theta)$")
label_est_cf <- TeX("$\\hat{\\psi}_{cf}(\\theta)$")
label_est_eif <- TeX("$\\hat{\\psi}(\\theta)$")
label_pw_ipw <- TeX("Pointwise 95% CI for $\\hat{\\psi}_{ipw}(\\theta)$")
label_est_ipw <- TeX("$\\hat{\\psi}_{ipw}(\\theta)$")

ggplot() +
  # --- Estimates ---
  geom_line(data = eif, aes(x = theta, y = est, color = "est_eif"), linetype = "dotted", 
            linewidth = 0.7, alpha = 0.9) +
  geom_line(data = cf, aes(x = theta, y = est, color = "est_cf"), 
            linewidth = 0.7, alpha = 0.5) +
  geom_line(data = ipw, aes(x = theta, y = est, color = "est_ipw"), linetype = "dashed", 
            linewidth = 0.7, alpha = 0.9) +
  
  # --- Pointwise CIs ---
  # eif Pointwise
  geom_line(data = eif, aes(x = theta, y = CI_L, color = "pw_eif"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.5) +
  geom_line(data = eif, aes(x = theta, y = CI_U, color = "pw_eif"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.5) +
  # cf Pointwise
  geom_line(data = cf, aes(x = theta, y = CI_L, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 0.5) +
  geom_line(data = cf, aes(x = theta, y = CI_U, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 0.5) +
  # ipw Pointwise
  geom_line(data = ipw, aes(x = theta, y = CI_L, color = "pw_ipw"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.5) +
  geom_line(data = ipw, aes(x = theta, y = CI_U, color = "pw_ipw"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.5) +
  
  # --- Uniform CI ---
  geom_ribbon(data = cf, aes(x = theta, ymin = CB_L, ymax = CB_U, color = "un_cf"), 
              fill = "lightblue", alpha = 0.3, linetype = 0) +
  
  # --- Others ---
  annotate("segment", x = 0.5, xend = 2, y = ybar, yend = ybar, 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = 0.5, xend = 2, y = ybarDNA, yend = ybarDNA, 
           color = "red", linetype = "dotted", linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray70") +
  geom_vline(xintercept = haz_ratio, linetype = "dotted", color = "gray70") +
  
  # --- Axis ---
  labs(y = "Estimated proportion of CIN2+ detected", x = expression(paste("Hazards ratio ", theta))) +
  scale_x_log10(breaks = c(0.5, 1, round(haz_ratio, 4), 2, 5, 0.5727, 0.7998),
                labels = c("0.5", "1", round(haz_ratio, 4), "2", "5", "0.5727", "0.7998")) +
  scale_y_continuous(limits = c(0, 0.14),
                     breaks = c(0.000, seq(0.025, 0.125, 0.025), ybar, ybarDNA),
                     labels = c("0.000", "0.025", "0.050", "0.075", "0.100", "0.125",
                                as.character(ybar), as.character(ybarDNA))) +
  
  # --- Label ---
  scale_color_manual(
    name = NULL,
    values = c("est_eif" = "slateblue", "est_cf" = "steelblue", "est_ipw" = "pink",
               "pw_eif" = "slateblue", "pw_cf" = "steelblue", "pw_ipw" = "pink",
               "un_cf" = "lightblue"),
    labels = c("est_eif" = label_est_eif, "est_cf" = label_est_cf, "est_ipw" = label_est_ipw, 
               "pw_eif" = label_pw_eif, "pw_cf" = label_pw_cf, "pw_ipw" = label_pw_ipw, 
               "un_cf" = label_un_cf),
    breaks = c("est_eif", "est_cf", "est_ipw", "pw_eif", "pw_cf", "pw_ipw", "un_cf")) +
  
  guides(
    color = guide_legend(
      ncol = 1, 
      override.aes = list(
        linetype = c("dotted", "solid", "dashed",
                     "solid", "longdash", "solid", 
                     "blank"), 
        fill     = c(NA, NA, NA, NA, NA, NA, "lightblue"), 
        linewidth = c(0.7, 0.7, 0.7, 
                      0.5, 0.5, 0.5, 
                      0),
        alpha    = c(0.9, 0.5, 0.9, 
                     0.9, 0.5, 0.9, 
                     0.3)))) +
  
  theme_classic() +
  theme(legend.position = c(0.02, 1), 
        legend.justification = c("left", "top"), 
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),     
        legend.text = element_text(size = 9),
        legend.key.width = unit(1.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        axis.line = element_line(color = "black", linewidth = 0.3)
  )
#-----------------------------------------------------




#-----------------------------------------------------
# plot
dfDNA <- readRDS("dfDNA.rds")
dfRNA <- readRDS("dfRNA.rds")
cf_res <- readRDS("ml.result.rds")
eif <- readRDS("eif.result.rds")
ipw <- readRDS("ipw.result.rds")
haz_ratio <- readRDS("haz_ratio.rds")
ybar <- round(eif$est[eif$theta==1], 3)
ybarDNA <- round(mean(dfDNA$Y), 3)

N <- nrow(dfRNA)
eif$CI_L <- eif$est + qnorm(0.025) * eif$sigma / sqrt(N)
eif$CI_U <- eif$est + qnorm(0.975) * eif$sigma / sqrt(N)

cf <- data.frame(theta = eif$theta, est = cf_res[[1]], 
                 CI_U = cf_res[["CI_U"]], CI_L = cf_res[["CI_L"]], 
                 CB_U = cf_res[["CB_U"]], CB_L = cf_res[["CB_L"]])
ext.idx <- 62
cf <- cf[1:ext.idx, ]
eif <- eif[1:ext.idx, ]
ipw <- ipw[1:ext.idx, ]

library(ggplot2)
library(latex2exp)
library(patchwork)

label_pw_cf  <- TeX("Pointwise 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_un_cf  <- TeX("Uniform 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_pw_eif <- TeX("Pointwise 95% CI for $\\hat{\\psi}(\\theta)$")
label_est_cf <- TeX("$\\hat{\\psi}_{cf}(\\theta)$")
label_est_eif <- TeX("$\\hat{\\psi}(\\theta)$")
label_pw_ipw <- TeX("Pointwise 95% CI for $\\hat{\\psi}_{ipw}(\\theta)$")
label_est_ipw <- TeX("$\\hat{\\psi}_{ipw}(\\theta)$")


p <- ggplot() +
  # --- Uniform CI ---
  geom_ribbon(data = cf, aes(x = theta, ymin = CB_L, ymax = CB_U, color = "un_cf"), 
              fill = "lightblue", alpha = 0.3, linetype = 0) +
  
  # --- Estimates ---
  geom_line(data = cf, aes(x = theta, y = est, color = "est_cf"), 
            linewidth = 0.7, alpha = 0.7) +
  geom_line(data = ipw, aes(x = theta, y = est, color = "est_ipw"), linetype = "dashed", 
            linewidth = 0.7, alpha = 0.9) +
  geom_line(data = eif, aes(x = theta, y = est, color = "est_eif"), linetype = "dotted", 
            linewidth = 0.7, alpha = 0.9) +
  
  # --- Pointwise CIs ---
  # eif Pointwise
  geom_line(data = eif, aes(x = theta, y = CI_L, color = "pw_eif"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.3) +
  geom_line(data = eif, aes(x = theta, y = CI_U, color = "pw_eif"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.3) +
  # cf Pointwise
  geom_line(data = cf, aes(x = theta, y = CI_L, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 0.3) +
  geom_line(data = cf, aes(x = theta, y = CI_U, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 0.3) +
  # ipw Pointwise
  geom_line(data = ipw, aes(x = theta, y = CI_L, color = "pw_ipw"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.3) +
  geom_line(data = ipw, aes(x = theta, y = CI_U, color = "pw_ipw"), 
            linetype = "solid", alpha = 0.9, linewidth = 0.3) +
  
  # --- Others ---
  annotate("segment", x = 0.5, xend = 1.1, y = ybarDNA, yend = ybarDNA, 
           color = "red", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = 0.5, xend = 1, y = ybar, yend = ybar, 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = 0.5, xend = round(haz_ratio, 3), 
           y = cf$est[19], yend = cf$est[19], 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  
  annotate("segment", x = 1, xend = 1, y = -Inf, yend = ybar, 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = haz_ratio, xend = haz_ratio, y = -Inf, yend = cf$est[19], 
           color = "gray70", linetype = "dotted", linewidth = 0.6) + 
  
  # --- Axis ---
  labs(y = "Estimated proportion of CIN2+ detected", x = expression(paste("Hazards ratio ", theta))) +
  scale_x_log10(breaks = c(0.5, 1, 1.1, round(haz_ratio, 3), 0.573, 0.8),
                labels = c("0.5", "1", "1.1", round(haz_ratio, 3), "0.573", "0.8")) +
  scale_y_continuous(#limits = c(0, 0.13),
    breaks = c(0.025, 0.075,
               ybar, round(cf$est[19], 3), ybarDNA),
    labels = c("0.025", "0.075",
               as.character(ybar), 
               as.character(round(cf$est[19], 3)),
               as.character(ybarDNA))) +
  
  # --- Label ---
  scale_color_manual(
    name = NULL,
    values = c("est_eif" = "green3", "est_cf" = "blue", "est_ipw" = "gold",
               "pw_eif" = "green3", "pw_cf" = "blue", "pw_ipw" = "gold",
               "un_cf" = "lightblue"),
    labels = c("est_eif" = label_est_eif, "est_cf" = label_est_cf, "est_ipw" = label_est_ipw, 
               "pw_eif" = label_pw_eif, "pw_cf" = label_pw_cf, "pw_ipw" = label_pw_ipw, 
               "un_cf" = label_un_cf),
    breaks = c("est_eif", "est_cf", "est_ipw", "pw_eif", "pw_cf", "pw_ipw", "un_cf")) +
  
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(
        linetype = c("dotted", "solid", "dashed",
                     "solid", "longdash", "solid", 
                     "blank"), 
        fill     = c(NA, NA, NA, NA, NA, NA, "lightblue"), 
        linewidth = c(0.7, 0.7, 0.7, 
                      0.3, 0.3, 0.3, 
                      0),
        alpha    = c(0.9, 0.7, 0.9, 
                     0.9, 0.5, 0.9, 
                     0.3)))) +
  
  theme_classic() +
  theme(legend.position = "right", #c(0.02, 1), 
        legend.justification = "top", #c("left", "top"), 
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),     
        legend.text = element_text(size = 10),
        legend.key.width = unit(1.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        #axis.line = element_line(color = "black", linewidth = 0.3),
        axis.line = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
  )

plot_spacer() + p + plot_spacer() +
  plot_layout(widths = c(2, 7, 2))
#-----------------------------------------------------



#-----------------------------------------------------
result <- cf
cand_set <- result$theta

round(100*(result$est[cand_set==2]/result$est[cand_set==1]-1), 1)
round(100*(1-result$est[cand_set==0.5]/result$est[cand_set==1]), 1)

round(result[c(152, 1, 19), ], 3)
#-----------------------------------------------------



#-----------------------------------------------------
# plot
dfDNA <- readRDS("dfDNA.rds")
dfRNA <- readRDS("dfRNA.rds")
cf_res <- readRDS("ml.result.rds")
eif <- readRDS("eif.result.rds")
ipw <- readRDS("ipw.result.rds")
haz_ratio <- readRDS("haz_ratio.rds")
ybar <- round(eif$est[eif$theta==1], 3)
ybarDNA <- round(mean(dfDNA$Y), 3)

N <- nrow(dfRNA)
eif$CI_L <- eif$est + qnorm(0.025) * eif$sigma / sqrt(N)
eif$CI_U <- eif$est + qnorm(0.975) * eif$sigma / sqrt(N)

cf <- data.frame(theta = eif$theta, est = cf_res[[1]], 
                 CI_U = cf_res[["CI_U"]], CI_L = cf_res[["CI_L"]], 
                 CB_U = cf_res[["CB_U"]], CB_L = cf_res[["CB_L"]])
ext.idx <- 62
cf <- cf[1:ext.idx, ]
eif <- eif[1:ext.idx, ]
ipw <- ipw[1:ext.idx, ]

library(ggplot2)
library(latex2exp)
library(patchwork)

label_pw_cf  <- TeX("Pointwise 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_un_cf  <- TeX("Uniform 95% CI for $\\hat{\\psi}_{cf}(\\theta)$")
label_pw_eif <- TeX("Pointwise 95% CI for $\\hat{\\psi}(\\theta)$")
label_est_cf <- TeX("$\\hat{\\psi}_{cf}(\\theta)$")
label_est_eif <- TeX("$\\hat{\\psi}(\\theta)$")
label_pw_ipw <- TeX("Pointwise 95% CI for $\\hat{\\psi}_{ipw}(\\theta)$")
label_est_ipw <- TeX("$\\hat{\\psi}_{ipw}(\\theta)$")


p <- ggplot() +
  # --- Uniform CI ---
  geom_ribbon(data = cf, aes(x = theta, ymin = CB_L, ymax = CB_U, color = "un_cf"), 
              fill = "lightblue", alpha = 0.3, linetype = 0) +
  
  # --- Estimates ---
  geom_line(data = cf, aes(x = theta, y = est, color = "est_cf"), 
            linewidth = 0.7, alpha = 0.7) +
  
  # --- Pointwise CIs ---
  # cf Pointwise
  geom_line(data = cf, aes(x = theta, y = CI_L, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 0.3) +
  geom_line(data = cf, aes(x = theta, y = CI_U, color = "pw_cf"), 
            linetype = "longdash", alpha = 0.5, linewidth = 0.3) +
  
  # --- Others ---
  annotate("segment", x = 0.5, xend = 1.1, y = ybarDNA, yend = ybarDNA, 
           color = "red", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = 0.5, xend = 1, y = ybar, yend = ybar, 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = 0.5, xend = round(haz_ratio, 3), 
           y = cf$est[19], yend = cf$est[19], 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  
  annotate("segment", x = 1, xend = 1, y = -Inf, yend = ybar, 
           color = "gray70", linetype = "dotted", linewidth = 0.6) +
  annotate("segment", x = haz_ratio, xend = haz_ratio, y = -Inf, yend = cf$est[19], 
           color = "gray70", linetype = "dotted", linewidth = 0.6) + 
  
  # --- Axis ---
  labs(y = "Estimated proportion of CIN2+ detected", x = expression(paste("Hazards ratio ", theta))) +
  scale_x_log10(breaks = c(0.5, 1, 1.1, round(haz_ratio, 3), 0.573, 0.8),
                labels = c("0.5", "1", "1.1", round(haz_ratio, 3), "0.573", "0.8")) +
  scale_y_continuous(#limits = c(0, 0.13),
    breaks = c(0.025, 0.075,
               ybar, round(cf$est[19], 3), ybarDNA),
    labels = c("0.025", "0.075",
               as.character(ybar), 
               as.character(round(cf$est[19], 3)),
               as.character(ybarDNA))) +
  
  # --- Label ---
  scale_color_manual(
    name = NULL,
    values = c("est_cf" = "blue", "pw_cf" = "blue", "un_cf" = "lightblue"),
    labels = c("est_cf" = label_est_cf, "pw_cf" = label_pw_cf, "un_cf" = label_un_cf),
    breaks = c("est_cf", "pw_cf", "un_cf")) +
  
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(
        linetype = c("solid", "longdash", "blank"), 
        fill     = c(NA, NA, "lightblue"), 
        linewidth = c(0.7, 0.3, 0),
        alpha    = c(0.7, 0.5, 0.3)))) +
  
  theme_classic() +
  theme(legend.position = c(0.02, 1), 
        legend.justification = c("left", "top"), 
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),     
        legend.text = element_text(size = 10),
        legend.key.width = unit(1.2, "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        #axis.line = element_line(color = "black", linewidth = 0.3),
        axis.line = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
  )

plot_spacer() + p + plot_spacer() +
  plot_layout(widths = c(2, 3, 2))
#-----------------------------------------------------

