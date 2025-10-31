################################################################################
##################### Simulation & Outcome Report ##############################

Sim_Process <- function(Running, B.sample, n.sample, ps.formula, p, sub_super_pop, trt_prob, out.TRUE.formula, our.FALSE.formula, alpha, Quasi_Separa){
  ATE.est <- rep(NULL, Running)
  SE.SummaryList <- vector("list", Running)
  CI.SummaryList <- vector("list", Running)
  Running_Skip <- data.frame(Std_QuasiSepara = 0, Strtf_QuasiSepara = 0)
  Warning.SummaryList <- vector("list", Running)
  pb <- txtProgressBar(min = 0, max = Running, style = 3)
  for (i in 1:Running) {
    sample_data <- Strtf.sampling(n.sample, superpopDATA = sub_super_pop,
                                  trt_prob_ind = trt_prob) 
    ATE.est[i] <- ATE_IPW2(sample_data, ps.formula, p)
    ##### SE Estimation (Bootstrap + Var Formula) ######
    # Var Formula
    VF.repSE <- Calc.SE_ATE(sample_data, ps.formula, p)
    # Standardized Bootstrap
    B.Std.results <- Std.Bootstrap(B.sample, ps.formula, p, sample_data)
    if(B.Std.results$Summary.Warning$Prop.Quasi_Separa > Quasi_Separa){
      Running_Skip$Std_QuasiSepara <-  Running_Skip$Std_QuasiSepara  + 1
      next
    }
    B.Std.rep <- Summary.Bootstrap(B.Std.results$B.ATE_results, ps.formula, p, sample_data, alpha)
    Warning.SummaryList[[i]] <- B.Std.results$Summary.Warning
    # Stratified Bootstrap
    B.Strtf.results <- Strtf.Bootstrap(B.sample, ps.formula, p, sample_data)
    if(B.Strtf.results$Summary.Warning$Prop.Quasi_Separa > Quasi_Separa){
      Running_Skip$Strtf_QuasiSepara <- Running_Skip$Strtf_QuasiSepara + 1
      next
    }
    B.Strtf.rep <- Summary.Bootstrap(B.Strtf.results$B.ATE_results, ps.formula, p, sample_data, alpha)
    SE.Summary.i <- cbind(Calc.SE_ATE(sample_data, ps.formula, p), B.Std.rep$SEs, B.Strtf.rep$SEs)
    SE.SummaryList[[i]] <- SE.Summary.i
    ##### CI Estimation  (Bootstrap + Var Formula) #########
    # Var Formula
    CF.repCI <- Calc.CI_ATE(SE_Summary = VF.repSE, alpha, sample_data, ps.formula, p) # Para CI (Austin's) based on Var Formula
    # Add Standardized Bootstrap to CIs
    CI.Summary.i <- rbind(CF.repCI, B.Std.rep$CIs, B.Strtf.rep$CIs)
    CI.SummaryList[[i]] <- CI.Summary.i
    setTxtProgressBar(pb, i)
  }
  close(pb)
  SE.results <- do.call(rbind, SE.SummaryList)
  CI.results <- do.call(rbind, CI.SummaryList)
  Warning.Summary <- do.call(rbind, Warning.SummaryList)
  return(list(ATE.Results = ATE.est, SE.Results = SE.results, CI.Results = CI.results, Warning.Results = Warning.Summary, 
              SkipRunning = Running_Skip/Running))
}


####### Parallel Processing 
Sim_Parallel.Process <- function(Running, B.sample, n.sample, ps.formula, p, sub_super_pop, 
                                 trt_prob, alpha, Quasi_Separa, weight, requested_cores, out.TRUE.formula, out.FALSE.formula) {
  cl <- makeCluster(requested_cores - 2)  
  clusterExport(cl, c( "trt_prob", "sub_super_pop", "n.sample", "p", "alpha", "ps.formula",
                      "B.sample", "Quasi_Separa", "weight", "out.TRUE.formula", "out.FALSE.formula"))
  clusterEvalQ(cl, {
    library(PSweight)
    library(rms)
    library(coxed)
    source("Point_Est.R")
    source("Var_Est(Sandwich).R")
    source("Var&CI_BootstrapEst.R")
    source("Extract_Super.R")
  })
  clusterEvalQ(cl, {
    print(ls()) 
  })
  results <- pblapply(1:Running, function(i) {
    sample_data <- Strtf.sampling(n.sample, superpopDATA = sub_super_pop, trt_prob_ind = trt_prob) 
    #### Point Estimate of WATE
    # table(sample_data$Z)
    ps_model <- glm(ps.formula, data = sample_data, family = binomial) # PS model formula
    sample_data$e_X <- predict(ps_model, type = "response")
    Point.Est_WATE <- Calc.Point_WATE(sample_data, sample_data$e_X, weight)
    Point.Est_AugT.WATE <- Calc.Point_AugWATE(sample_data, sample_data$e_X, weight, out.TRUE.formula)
    Point.Est_AugF.WATE <- Calc.Point_AugWATE(sample_data, sample_data$e_X, weight, out.FALSE.formula)
    Point_Est_i <- as.data.frame(cbind(Point.Est_WATE, Point.Est_AugT.WATE, Point.Est_AugF.WATE))
    #### Var Formula Est of SE
    VF_repSE <- cbind(Calc.WATE_SE(sample_data, ps.formula, p, weight), 
                  Calc.AugT_WATE_SE(sample_data, ps.formula, out.TRUE.formula, weight), 
                  Calc.AugF_WATE_SE(sample_data, ps.formula, out.FALSE.formula, weight))
    #### Bootstrap Est of SE
    # Standardized Bootstrap
    B_Std_results <- Std.Bootstrap(B.sample, ps.formula, p, sample_data, weight, out.TRUE.formula, out.FALSE.formula)
    Warning_i <- B_Std_results$Summary.Warning
    if (Warning_i$Prop.Quasi_Separa > Quasi_Separa) {
      return(list(Skip = TRUE, Point_Est = NA, SE_Summary = NA, CI_Summary = NA, Warning_Summary = Warning_i))
    }
    B_Std_rep <- Summary.Bootstrap(B_Std_results$B.Point_results,  ps.formula, p, sample_data, alpha, weight, out.TRUE.formula, out.FALSE.formula)
    SE_Summary_i <- cbind(VF_repSE, B_Std_rep$SEs)
    #### CI Report ####
    # Var Formula CI
    CF_repCI <- Calc.CI_Par(SE_Summary = SE_Summary_i, alpha, sample_data, ps.formula, weight, out.TRUE.formula, out.FALSE.formula)
    # Combine CIs
    CI_Summary_i <- rbind(CF_repCI, B_Std_rep$CIs)
    list(Skip = FALSE, Point_Est = Point_Est_i, SE_Summary = SE_Summary_i, CI_Summary = CI_Summary_i, Warning_Summary = Warning_i)
  }, cl = cl)
  stopCluster(cl)
  
  # Initialize summary variables
  Running_Skip <- data.frame(Std_QuasiSepara = 0)
  Point_resultsList <- vector("list", Running)
  SE_SummaryList <- vector("list", Running)
  CI_SummaryList <- vector("list", Running)
  Warning_SummaryList <- vector("list", Running)
  
  # Process results
  for (i in 1:Running) {
    result <- results[[i]]
    if (result$Skip) {
      Running_Skip$Std_QuasiSepara <- Running_Skip$Std_QuasiSepara + 1
      Point_resultsList[[i]] <- NA
      SE_SummaryList[[i]] <- NA
      CI_SummaryList[[i]] <- NA
      Warning_SummaryList[[i]] <- result$Warning_Summary
    } else {
      Point_resultsList[[i]] <- result$Point_Est
      SE_SummaryList[[i]] <- result$SE_Summary
      CI_SummaryList[[i]] <- result$CI_Summary
      Warning_SummaryList[[i]] <- result$Warning_Summary
    }
  }
  Point_results <- do.call(rbind, Point_resultsList)
  SE_results <- do.call(rbind, SE_SummaryList)
  CI_results <- do.call(rbind, CI_SummaryList)
  Warning_Summary <- do.call(rbind, Warning_SummaryList)
  
  return(list(Point.Results = Point_results, SE.Results = SE_results, 
              CI.Results = CI_results, Warning.Results = Warning_Summary, 
              SkipRunning.Prop = Running_Skip / Running))
}

# RST<- Sim_Parallel.Process(Running = 100, B.sample = 100, n.sample = n.sample, ps.formula = ps.formula, p = p, sub_super_pop = sub_super_pop, 
#                            trt_prob = trt_prob, alpha = 0.05, Quasi_Separa = 0.1, weight = weight, requested_cores = 8, 
#                            out.TRUE.formula = out.TRUE.formula, out.FALSE.formula = out.FALSE.formula)
# #CI <- Sim_Report(RST, Ture_Estimate)

# Sim.Results<-RST
######### Input Sim-Process results Sim.Results of a list(ATE, SE, and CI), output summary
Sim_Report <- function(Sim.Results, Ture_Estimate){
  Point.Est <- Sim.Results$Point.Results
  ##### ATE Summary Report  ##### 
  Empirical.SE <- sd(Point.Est$Point.Est_WATE, na.rm = TRUE) ### true SE of Estimated by WATE (not Augmented)
  average_results <- data.frame(t(lapply(Point.Est, mean, na.rm = TRUE)))
  colnames(average_results) <- paste0("Mean_", colnames(Point.Est))
  Point.Report <- data.frame(TRUE_Point = Ture_Estimate, Avg_Point.Est = average_results, Empirical.SE = Empirical.SE)
  ##### SE Summary Report  #####
  SE.Report <- data.frame(
    Mean.Est.SE = apply(Sim.Results$SE.Results, 2, function(x) mean(x, na.rm = TRUE)),
    Ratio.SE = apply(Sim.Results$SE.Results, 2, function(x) mean(x, na.rm = TRUE)) / Empirical.SE,
    SE.SE = apply(Sim.Results$SE.Results, 2, function(x) sd(x, na.rm = TRUE))/sqrt(nrow(Sim.Results$SE.Results))
  )
  SE.Report <- cbind(Method = rownames(SE.Report), SE.Report)
  rownames(SE.Report) <- NULL
  ##### CI Summary Report ######
  CI.Report <- Sim.Results$CI.Results %>%
    filter(!is.na(Method) & Method != "") %>%
    group_by(Method) %>%
    summarise(
      Coverage_Prob = mean(Ture_Estimate >= Lower & Ture_Estimate <= Upper, na.rm = TRUE),
      Width_Avg = mean(Upper - Lower, na.rm = TRUE),
      LMiss_Prob = mean(Ture_Estimate < Lower, na.rm = TRUE),
      RMiss_Prob = mean(Ture_Estimate > Upper, na.rm = TRUE),
      Upper_Avg = mean(Upper, na.rm = TRUE),
      Lower_Avg = mean(Lower, na.rm = TRUE)
    )
  Warning.Report <- Sim.Results$Warning.Results %>%
    summarise(
      Mean_Prop.QuaSp = mean(Prop.Quasi_Separa, na.rm = TRUE),
      Max_Prop.QuaSp = max(Prop.Quasi_Separa, na.rm = TRUE),
      Mean_Prop.NonConvg = mean(Prop.Non_Convg, na.rm = TRUE),
      Max_Prop.NonConvg = max(Prop.Non_Convg, na.rm = TRUE),
      Mean_Prop.validB = mean(Prop.validBoostrap, na.rm = TRUE)
    )
  Skip.Report <- Sim.Results$SkipRunning
  return(list(Point.Report = Point.Report, SE.Report = SE.Report, CI.Report= CI.Report, 
              Boost.Warning.Report = Warning.Report, Skip.Report = Skip.Report))
}

