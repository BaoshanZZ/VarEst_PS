####### Point Estimate of WATE from Sample (Option: IPW and overlap) #######
Calc.Point_WATE <- function(sample, ps.estimate, weight = "IPW") {
  N.sample <- nrow(sample) 
  e_X <- ps.estimate
  if (weight == "IPW") {
    omega <- rep(1, N.sample)
  }  else if  (weight == "overlap") {
    omega <- e_X * (1-e_X)
  }
  Z <- sample$Z; Y <- sample$Y_obs
  W <- omega / (Z * e_X + (1 - Z) * (1 - e_X))
  mu_1 <- sum(W * Z * Y) / sum(W * Z)
  mu_2 <- sum(W * (1 - Z) * Y) / sum(W * (1 - Z))
  return(mu_1 - mu_2)
}

Calc.Point_AugWATE <- function(sample, ps.estimate, weight = "IPW", out.formula) {
  N.sample <- nrow(sample) 
  e_X <- ps.estimate
  if (weight == "IPW") {
    omega <- rep(1, N.sample)
  } else if  (weight == "overlap") {
    omega <- e_X * (1 - e_X)}
  Z <- sample$Z; Y <- sample$Y_obs
  W <- omega / (Z * e_X + (1 - Z) * (1 - e_X))
  m1 <- glm(out.formula, data = sample[Z == 1,], family = binomial)
  m0 <- glm(out.formula, data = sample[Z == 0,], family = binomial)
  pred_m1 <- predict(m1, newdata = sample, type = "response")
  pred_m0 <- predict(m0, newdata = sample, type = "response")

  Term1 <- sum(omega * (pred_m1 - pred_m0)) / sum(omega)
  Term2 <- sum(W * Z * (Y - pred_m1)) / sum(W * Z)
  Term3 <- sum(W * (1 - Z) * (Y - pred_m0)) / sum(W * (1 - Z))
  WATE_Aug <- Term1 +  Term2 - Term3
  return(WATE_Aug)
}

