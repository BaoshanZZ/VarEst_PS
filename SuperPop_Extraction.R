extract_superpop <- function(trt_prob_ind, out_prob_ind, Full_super_pop) {
  Z <- Full_super_pop[[paste0("Z_", trt_prob_ind, "0")]]
  popPS <- Full_super_pop[[paste0("popPS_", trt_prob_ind, "0")]]
  Y0 <- Full_super_pop[[paste0("Y0_", out_prob_ind, "0")]]
  Y1 <- Full_super_pop[[paste0("Y1_", out_prob_ind, "0")]]
  Y_obs <- Z * Y1 + (1-Z) * Y0
  features <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  super_pop <- data.frame(Full_super_pop[features], Z = Z, popPS = popPS, Y0 = Y0, Y1 = Y1, Y_obs = Y_obs)
  return(super_pop)
}

estimand <- function(super_pop, weight){
  if (weight == "IPW") {
    Estimand <- mean(super_pop$Y1-super_pop$Y0)
  }else if (weight == "overlap"){
    mu_1 <- sum(super_pop$Z*super_pop$Y_obs*(1-super_pop$popPS))/sum(super_pop$Z*(1-super_pop$popPS))
    mu_2 <- sum((1-super_pop$Z) * super_pop$Y_obs * super_pop$popPS)/sum((1-super_pop$Z) * super_pop$popPS)
    Estimand <- mu_1 - mu_2
  }
  return(Estimand)
}

# Stratified Random Sample method by proportion of treatment
Strtf.sampling <- function(n, superpopDATA, trt_prob_ind){
  # trt_prob_ind is an integral, e.g., 20% we input 2
  trt_prob <- trt_prob_ind / 10
  N.super_pop <- nrow(superpopDATA)
  # Calculate the number of treated and control samples
  n_treated <- floor(n * trt_prob)
  n_control <- n - n_treated
  
  # Stratified sampling
  treated_indices <- sample(which(superpopDATA$Z == 1), n_treated, replace = FALSE)
  control_indices <- sample(which(superpopDATA$Z == 0), n_control, replace = FALSE)
  
  # Combine sampled indices
  sample_indices <- c(treated_indices, control_indices)
  sample <- superpopDATA[sample_indices, ]
  
  return(sample)
}
