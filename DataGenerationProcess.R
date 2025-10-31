##### Data Generating Process (Austin PC, 2022) #####
rm(list = ls())
library(MASS)
library(dplyr)
getwd()
#setwd("/Users/baoshanzhang/Library/CloudStorage/OneDrive-DukeUniversity/Research/DCRI/Code/Super_Pop/")
setwd("/hpc/home/bz91/DCRI_VarEst")
set.seed(111)

n <- 1000000  # Population size
mean_vector <- rep(0, 10)
cov_matrix <- matrix(0.2, nrow = 10, ncol = 10)
diag(cov_matrix) <- 1

baseline_covariates <- mvrnorm(n, mu = mean_vector, Sigma = cov_matrix)

x1 <- baseline_covariates[, 1]
x2 <- baseline_covariates[, 2]
x3 <- baseline_covariates[, 3]
x4 <- baseline_covariates[, 4]
x5 <- baseline_covariates[, 5]
x6 <- as.numeric(baseline_covariates[, 6] < quantile(baseline_covariates[, 6], 0.1))
x7 <- as.numeric(baseline_covariates[, 7] < quantile(baseline_covariates[, 7], 0.2))
x8 <- as.numeric(baseline_covariates[, 8] < quantile(baseline_covariates[, 8], 0.3))
x9 <- as.numeric(baseline_covariates[, 9] < quantile(baseline_covariates[, 9], 0.4))
x10 <- as.numeric(baseline_covariates[, 10] < quantile(baseline_covariates[, 10], 0.5))

logit <- function(x) {exp(x) / (1 + exp(x))}

# Function to generate population_level PS
Pop_PS <- function(alpha0) {
  log_odds <- alpha0 + log(1.1) * x1 + log(1.2) * x2 + log(1.5) * x3 + log(1.75) * x4 + log(2) * x5 + log(1.25) * x6 + 
    log(1.5) * x7 + log(2) * x8 + log(0.8) * x9 + log(0.5) * x10
  p_treat <- logit(log_odds)
  return(p_treat)
}

# Function to generate treatment variable
generate_treatment <- function(alpha0) {
  Z <- rbinom(n, 1, Pop_PS(alpha0))
  return(Z)
}

pop_data <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
for (i in 1:9) {
  assign(paste0("pop_data", i), data.frame(pop_data))
}


########## PS MODEL ###########
# Prevalence of treatment 10% ~ 90%
prevalences <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
alpha0_treat <- c(-2.69, -1.69, -1.00, -0.42, 0.13, 0.67, 1.25, 1.95, 2.93)


for (i in 1:9) {
  alpha0 <- alpha0_treat[i]
  pop_data[[paste0("Z_",i,"0")]] <- generate_treatment(alpha0)
  pop_data[[paste0("popPS_",i,"0")]] <- Pop_PS(alpha0)
}
head(pop_data)

########## OUTCOME MODEL ###########
# Define the logistic model for outcome
generate_outcome <- function(alpha0, alpha_treat, Z) {
  log_odds <- alpha0 + alpha_treat * Z + log(2) * x1 + log(1.75) * x2 + log(1.1) * x3 + log(1.5) * x4 + log(1.2) * x5 + log(2) * x6 
  + log(1.5) * x7 + log(1.1) * x8 + log(1.25) * x9 + log(2) * x10
  p_outcome <- logit(log_odds)
  Y <- rbinom(n, 1, p_outcome)
  return(Y)
}

#Adjust based on desired prevalence of outcome to be 0.1~0.5 (untreated arm) and ATE as -0.02
alpha0_outcome <- c(-2.79, -1.83, -1.16, -0.59, -0.07)
alpha_treat <- c(-0.28, -0.16, -0.125, -0.115, -0.113)
Z_0 <- rep(0, n); Z_1 <- rep(1, n)

for (i in 1:5) {
  alpha0_out <- alpha0_outcome[i]
  alpha_trt <- alpha_treat[i]
  pop_data[[paste0("Y0_",i,"0")]] <- generate_outcome(alpha0_out, alpha_trt, Z_0)
  pop_data[[paste0("Y1_",i,"0")]] <- generate_outcome(alpha0_out, alpha_trt, Z_1)
}

head(pop_data)

# table(pop_data$Y0_20)
mean(pop_data$Y1_50-pop_data$Y0_50)

write.csv(pop_data, file = "SuperPop.csv")
