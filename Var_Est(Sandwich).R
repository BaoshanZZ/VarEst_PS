#' ---
#' title: "Asymptotic Sandwich Variance Estimators for WATE"
#' description: >
#'   R functions to calculate and compare asymptotic sandwich variance
#'   estimators for the Weighted Average Treatment Effect (WATE).
#'   Includes implementations for standard IPW and Overlap weighting,
#'   as well as Augmented WATE (A-WATE) estimators.
#'   The script compares manually implemented estimators (PSE, L&D, Li)
#'   with the numerical methods from the PSweight package.
#' ---

#' Calculate Standard Errors for WATE (IPW and Overlap)
#'
#' Computes various asymptotic sandwich variance estimators for the WATE.
#' It contrasts estimators assuming fixed (known) PS
#' versus those that account for the estimation uncertainty of the PS model.
#'
#' @param sample data.frame. The input dataset. *Must* contain a
#'   pre-calculated propensity score column named `e_X`, treatment indicator `Z`,
#'   and observed outcome `Y_obs`.
#' @param ps.formula formula. The formula used for the propensity score model
#'   (e.g., `Z ~ X1 + X2`). This is used by `PSweight` and to identify
#'   covariates.
#' @param p integer. The number of covariates (parameters) in the PS model.
#'   *Note: This is used to subset columns (e.g., `sample[i, 1:p]`), which is
#'   a fragile practice. See refinement suggestions.*
#' @param weight character. The weighting method: "IPW" or "overlap".
#'
#' @return A data.frame summarizing the Standard Errors (SE) from
#'   different methods (e.g., `PSw.Fixed`, `PES.Est`, `PSw.Est`).
#'
#' @details
#'   - "Fixed" methods treat the `e_X` propensity score as known.
#'   - "Est" (Estimated) methods account for the uncertainty from
#'     *estimating* the PS model parameters via a sandwich estimator.
#'   - `PES`: Purely Empirical Sandwich estimator (manual implementation).
#'   - `PES_LD`: Empirical Sandwich estimator from  L&D (2017), Eq(19) (manual).
#'   - `PES_Li`: Empirical Sandwich estimator from Li et al. (2018) for Overlap (manual).
#'   - `PSw`: Estimator from the `PSweight` package.
#'

Calc.WATE_SE <- function(sample, ps.formula, p, weight = "overlap") {
  ### ps.formula, p: # of para of PS model 
  n <- nrow(sample)
  ## Propensity Model (Logistics)
  # ps_model <- glm(ps.formula, data = sample, family = binomial) # PS model formula
  # sample$e_X <- predict(ps_model, type = "response")
  ####################POINT ESTIMATOR L&D IPW2############################
  Z <- sample$Z
  e <- sample$e_X
  Y <- sample$Y_obs
  # Calculate the treated and control weights
  if (weight == "IPW") {
    # Calculate the treated and control weights for IPW
    w_trt <- Z / e 
    w_ctrl <- (1 - Z) / (1 - e)
    # Calculate the weighted ATE
    mu_1 <- sum(w_trt * Y) / sum(w_trt)
    mu_0 <- sum(w_ctrl * Y) / sum(w_ctrl)
    out.res_trt <- Z * (Y - mu_1)/e
    out.res_ctrl <- (1 - Z) * (Y - mu_0) / (1 - e) 
  } else if (weight == "overlap") {
    # Calculate the treated and control weights for overlap weighting
    w_trt <- Z * (1 - e) 
    w_ctrl <- (1 - Z) * e
    # Calculate the weighted ATE
    mu_1 <- sum(w_trt * Y) / sum(w_trt)
    mu_0 <- sum(w_ctrl * Y) / sum(w_ctrl)
    out.res_trt <- Z * (Y - mu_1)*(1-e)
    out.res_ctrl <- (1 - Z) * (Y - mu_0)* e 
  } else {
    stop("Invalid weight option. Choose 'IPW' or 'overlap'.")
  }
  w_trt_Avg<- mean(w_trt)
  w_ctrl_Avg <- mean(w_ctrl)

  ###################################################################################################
  ################################## VAR ESTIMATOR of WATE (PSE, L&D);##############################
  ###### CASE A: PS Model Correct, Fixed PS Model Coefficient ######
  # ### A1: Derived Purely empirical sandwich estimator (ONLY IPW)#####
  # Var_PES.fixed <- w_trt_Avg^{-2} * mean(out.res_trt^2)  +  w_ctrl_Avg^{-2} * mean(out.res_ctrl^2)
  # SE_PES.fixed <- sqrt(Var_PES.fixed/n)
  
  #### A2: PSweight() function with fixed coefficient ###
  m.fixed <- PSweight(ps.estimate = sample$e_X, yname = "Y_obs", data = sample, weight = weight, zname = "Z")
  SE_PSw.fixed <- summary(m.fixed)$estimates[, "Std.Error"]
  #########################################################################################
  #########################################################################################
  ########## CASE B: PS Model Correct, Estimated PS Model Coefficient ####################
  if(weight == "IPW"){
    out.res_trt <- Z * (Y-mu_1)/e
    out.res_ctrl <- (1 - Z) * (Y - mu_0) / (1 - e) 
    ### B1 (Only for ATE): Purely empirical sandwich estimator (PSE)  #####
    E_bb <- Reduce("+", lapply(1:n, function(i) {
      # E_bb.inv are the same for B1 and B2
      e_X_i <- sample$e_X[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      (e_X_i * (1 - e_X_i)) * (X_i %*% t(X_i))
    })) / n
    E_bb.inv <- solve(E_bb)
    H_beta_1 <- Reduce("+", lapply(1:n, function(i) {
      e_X_i <- sample$e_X[i]
      Z_i <- sample$Z[i]
      Y_i <- sample$Y_obs[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      term1 <- Z_i * (Y_i - mu_1) * (1 - e_X_i) / e_X_i
      term2 <- (1 - Z_i) * (Y_i - mu_0) * e_X_i / (1 - e_X_i)
      (w_trt_Avg^{-1} * term1 + w_ctrl_Avg^{-1} * term2) * X_i
    })) / n
    Var_PES.est <-  Reduce("+", lapply(1:n, function(i) {
      e_X_i <- sample$e_X[i]
      Z_i <- sample$Z[i]
      out.res_trt_i <- out.res_trt[i]
      out.res_ctrl_i <- out.res_ctrl[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      (w_trt_Avg^{-1} * out.res_trt_i - w_ctrl_Avg^{-1} * out.res_ctrl_i - (Z_i-e_X_i)%*%t(H_beta_1)%*%E_bb.inv%*%X_i)^2
    })) / n
    SE_PES.est <- sqrt(Var_PES.est/n)
    colnames(SE_PES.est) <- "PES.Est"
    
    ### B2: Empirical sandwich estimator from L&D (19) #####
    H_beta_2 <- Reduce("+", lapply(1:n, function(i) {
      e_X_i <- sample$e_X[i]
      Z_i <- sample$Z[i]
      Y_i <- sample$Y_obs[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      term1 <- Z_i * (Y_i - mu_1) * (1 - e_X_i) / e_X_i
      term2 <- (1 - Z_i) * (Y_i - mu_0) * e_X_i / (1 - e_X_i)
      (term1 + term2) * X_i
    })) / n
    Var_PES.estLD <-  Reduce("+", lapply(1:n, function(i) {
      e_X_i <- sample$e_X[i]
      Z_i <- sample$Z[i]
      out.res_trt_i <- out.res_trt[i]
      out.res_ctrl_i <- out.res_ctrl[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      (out.res_trt_i -out.res_ctrl_i - (Z_i-e_X_i)%*%t(H_beta_2)%*%E_bb.inv%*%X_i)^2
    })) / n
    SE_PES.estLD <- sqrt(Var_PES.estLD/n)
    colnames(SE_PES.estLD) <- "PES_LD.Est"
  } else if (weight == "overlap"){
    ########### B1: Li's Sandwich Estimator  ###########
    E_bb <- Reduce("+", lapply(1:n, function(i) {
      e_X_i <- sample$e_X[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      (e_X_i * (1 - e_X_i)) * (X_i %*% t(X_i))
    })) / n
    E_bb.inv <- solve(E_bb)
    H_beta <- Reduce("+", lapply(1:n, function(i) {
      e_X_i <- sample$e_X[i]
      Z_i <- sample$Z[i]
      Y_i <- sample$Y_obs[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      term1 <- Z_i * (Y_i - mu_1)     
      term2 <- (1 - Z_i) * (Y_i - mu_0) 
      (term1 + term2) * e_X_i * (1 - e_X_i) * X_i
    })) / n
    Var_Li.est <-  Reduce("+", lapply(1:n, function(i) {
      e_X_i <- sample$e_X[i]
      Z_i <- sample$Z[i]
      out.res_trt_i <- out.res_trt[i]
      out.res_ctrl_i <- out.res_ctrl[i]
      X_i <- t(as.matrix(sample[i, 1:p]))
      (out.res_trt_i - out.res_ctrl_i - (Z_i-e_X_i)%*%t(H_beta)%*%E_bb.inv%*%X_i)^2
    })) / n
    theta.hat <- mean(e * (1-e))
    SE_Li.est <- sqrt(Var_Li.est/n)/theta.hat
    colnames(SE_Li.est) <- "PES_Li.Est"
  }
  ### B3: PSweight() function with estimated coefficient (w/o Augmented by Outcome Model)###
  m.estimated <- PSweight(ps.formula = ps.formula, yname = "Y_obs", data = sample, weight = weight, zname = "Z")
  SE_PSw.est <- summary(m.estimated)$estimates[, "Std.Error"]
  ### Summary of Estimated SE  ###
  if (weight == "IPW") {
    SE_Summary <- data.frame(PSw.Fixed = SE_PSw.fixed, 
                             PES.Est = SE_PES.est, PES_LD.Est= SE_PES.estLD, PSw.Est = SE_PSw.est)
  }else if (weight == "overlap")
  {
    SE_Summary <- data.frame(PSw.Fixed = SE_PSw.fixed, 
                             PES_Li.Est= SE_Li.est, PSw.Est = SE_PSw.est) 
  }
  return(SE_Summary)
}


#' Calculate Standard Errors for Augmented WATE (TRUE Outcome Model)
Calc.AugT_WATE_SE <- function(sample, ps.formula, out.TRUE.formula, weight = "IPW") {
  ### ps.formula, p: # of para of PS model 
  n <- nrow(sample)
  ## Propensity Model (Logistics)
  # ps_model <- glm(ps.formula, data = sample, family = binomial) # PS model formula
  # sample$e_X <- predict(ps_model, type = "response")
  # Sample$e_X contains in Sample
  ###################################################################################################
  ##### ########################### VAR ESTIMATOR of Augmented WATE #################
  ###### CASE A: PS Model Correct, Fixed PS Model Coefficient ######
  #### A3: Augmented by Outcome Model (Consider TURE/Misspecifed Outcome model). Fixed PS Model #####
  ## Augmented by TRUE Outcome MODEL
  Aug.T.Fixed <- PSweight(ps.estimate = sample$e_X, out.formula = out.TRUE.formula, 
                          augmentation = TRUE, yname = "Y_obs", data = sample, weight = weight, zname = "Z", family = "binomial")
  SE_TRUE.Aug.PSw.Fixed <- summary(Aug.T.Fixed)$estimates[, "Std.Error"]
  #########################################################################################
  #########################################################################################
  ########## CASE B: PS Model Correct, Estimated PS Model Coefficient #####################
  #### B4: Augmented by Outcome Model (Consider TURE/Mis-specified Outcome model) with Est. PS Model #####
  ## Augmented by TRUE Outcome MODEL
  Aug.T.Est <- PSweight(ps.formula = ps.formula, out.formula = out.TRUE.formula, 
                        augmentation = TRUE, yname = "Y_obs", data = sample, weight = weight, zname = "Z", family = "binomial")
  SE_TRUE.Aug.PSw.Est <- summary(Aug.T.Est)$estimates[, "Std.Error"]
  ### Summary of Estimated SE  ###
  SE_Summary <- data.frame( PSw.Fixed.AugT = SE_TRUE.Aug.PSw.Fixed, PSw.Est.AugT = SE_TRUE.Aug.PSw.Est) 
  return(SE_Summary)
}

#' Calculate Standard Errors for Augmented WATE (False Outcome Model)
Calc.AugF_WATE_SE <- function(sample, ps.formula, out.FALSE.formula, weight = "IPW") {
  ### ps.formula, p: # of para of PS model 
  n <- nrow(sample)
  ## Propensity Model (Logistics Regression)
  # ps_model <- glm(ps.formula, data = sample, family = binomial) # PS model formula
  # sample$e_X <- predict(ps_model, type = "response")
  ############################ VAR ESTIMATOR of Augmented WATE ##############################
  ###### CASE A: PS Model Correct, Fixed PS Model Coefficient ################################
  #### A3: Augmented by Outcome Model (Consider TURE/Mis-specified Outcome model). Fixed PS Model #####
  Aug.F.Fixed <- PSweight(ps.estimate = sample$e_X, out.formula = out.FALSE.formula, 
                          augmentation = TRUE, yname = "Y_obs", data = sample, weight = weight, zname = "Z", family = "binomial")
  SE_FALSE.Aug.PSw.fixed <- summary(Aug.F.Fixed)$estimates[, "Std.Error"]
  #########################################################################################
  ########## CASE B: PS Model Correct, Estimated PS Model Coefficient #####################
  #### B4: Augmented by Outcome Model (Consider TURE/Mis-specified Outcome model) with Est. PS Model #####
  Aug.F.Est <- PSweight(ps.formula = ps.formula, out.formula = out.FALSE.formula, 
                        augmentation = TRUE, yname = "Y_obs", data = sample, weight = weight, zname = "Z", family = "binomial")
  SE_FALSE.Aug.PSw.Est <- summary(Aug.F.Est)$estimates[, "Std.Error"]
  ### Summary of Estimated SE  ###
  SE_Summary <- data.frame(PSw.Fixed.AugF = SE_FALSE.Aug.PSw.fixed, 
                           PSw.Est.AugF = SE_FALSE.Aug.PSw.Est) 
  return(SE_Summary)
}

#' Calculate Parametric (Wald) Confidence Intervals
#' Calculates 95% (or other alpha) confidence intervals based on the
#' point estimates and the various standard errors computed.
#'
Calc.CI_Par <- function(SE_Summary, alpha, sample, ps.formula, weight, out.TRUE.formula, out.FALSE.formula) {
  # ps_model <- glm(ps.formula, data = sample, family = binomial)
  # sample$e_X <- predict(ps_model, type = "response")
  # Point Estimate Summary based on Random Sample
  Point.Est_WATE <- Calc.Point_WATE(sample, sample$e_X, weight)
  Point.Est_AugT.WATE <- Calc.Point_AugWATE(sample, sample$e_X, weight, out.TRUE.formula)
  Point.Est_AugF.WATE <- Calc.Point_AugWATE(sample, sample$e_X, weight, out.FALSE.formula)
  z_value <- qnorm(1 - alpha / 2)
  CIs <- data.frame(Method = character(), Lower = numeric(), Upper = numeric())
  for (col in names(SE_Summary)) {
    if (grepl("AugT", col)) { Point.Est <- Point.Est_AugT.WATE
    } else if (grepl("AugF", col)) { Point.Est <- Point.Est_AugF.WATE
    } else { Point.Est <- Point.Est_WATE }
    SE_value <- SE_Summary[[col]][1]
    lower_bound <- Point.Est - z_value * SE_value
    upper_bound <- Point.Est + z_value * SE_value
    CIs <- rbind(CIs, data.frame(Method = paste0(col, "_Par"), Lower = lower_bound, Upper = upper_bound))
  }
  return(CIs)
}
