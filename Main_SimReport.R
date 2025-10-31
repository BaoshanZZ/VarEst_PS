rm(list = ls())
getwd()
setwd("/hpc/home/bz91/DCRI_VarEst")
library(dplyr)
library(PSweight)
library(rms)
library(coxed)
library(parallel)
library(pbapply)  # for progress bar in parallel
source("Point_Est.R")
source("Var_Est(Sandwich).R")
source("Var&CI_BootstrapEst.R")
source("Extract_Super.R")
source("Sim_Function.R")

Full_SuperPop <- read.csv("Super_Pop/SuperPop.csv")
N.super_pop <- nrow(Full_SuperPop)[1]

ps.formula <- Z ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
p <- 10
out.TRUE.formula <- Y_obs ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 
out.FALSE.formula <- Y_obs ~ x1 + x2 + x3 +  x6 + x10

Running <- 10000
B.sample <- 1000
alpha <- 0.05
Quasi_Separa <- 0.1
requested_cores <- 32

# Demo Code for this main report
for(weight in c("IPW")){
  for (n.sample in c(200)) {
    for (out_prob in c(5:2)) {
      for (trt_prob in c(1)) {
        sub_super_pop <- extract_superpop(out_prob_ind = out_prob, trt_prob_ind = trt_prob , Full_super_pop = Full_SuperPop)
        Estimand <- estimand(super_pop = sub_super_pop, weight)
        rst_name <- paste0("RstY", out_prob, "0Z", trt_prob, "0n", n.sample)
        rep_name <- paste0("RepY", out_prob, "0Z",trt_prob, "0n", n.sample)
        Sim.results <- Sim_Parallel.Process(Running, B.sample, n.sample, 
                                            ps.formula, p, sub_super_pop, trt_prob, alpha, Quasi_Separa, weight, requested_cores,
                                            out.TRUE.formula, out.FALSE.formula)
        assign(rst_name, Sim.results)
        assign(rep_name, Sim_Report(Sim.results, Estimand))
        save(list = c(rst_name, rep_name), file = paste0("Rdata_WATE_Augment/", weight,"/Report_Y", out_prob, "0Z", trt_prob, "0n", n.sample, ".RData"))
        print(c(n.sample, out_prob, trt_prob))}
    }
  }
}
