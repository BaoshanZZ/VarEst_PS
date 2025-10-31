
handle_warning <- function(w, b) {
  if (grepl("glm.fit: algorithm did not converge", w$message)) {
    count.non_convg <<- count.non_convg + 1
  }
  if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", w$message)) {
    count.quasi_separa <<- count.quasi_separa + 1
  }
  return(NULL)
}

########### Bootstrap Variance Estimation Method ################
## Generate Standardized Bootstrap Sample and get summary for each B.sample
## Fixed and Estimated Propensity Weighting are applied
Std.Bootstrap <- function(B.sample, ps.formula, p, sample, weight = "IPW", out.TRUE.formula, out.FALSE.formula){
    n.sample <- nrow(sample)
    # # PS model for random sample
    # ps_model <- glm(ps.formula, data = sample, family = binomial)
    # sample$e_X <- predict(ps_model, type = "response")
    B.Point <- data.frame(StdB.Fixed = numeric(), StdB.Est = numeric(),
                          StdB.Fixed.AugT = numeric(), StdB.Est.AugT = numeric(),
                          StdB.Fixed.AugF = numeric(), StdB.Est.AugF = numeric())
    # Initialize warning counter, global variable
    count.non_convg <<- 0
    count.quasi_separa <<- 0
    count.validBoostrap <- 0
    b <- 0
    #pb <- txtProgressBar(min = 0, max = B.sample, style = 3)
    while (count.validBoostrap < B.sample)  {
      B_sample.ind <- sample(1:n.sample, size = n.sample, replace = TRUE)
      B_sample <- sample[B_sample.ind, ]
      if (nrow(B_sample[B_sample$Z == 1, ]) == 0 || nrow(B_sample[B_sample$Z == 0, ]) == 0) {
        next
      }
      b <- b + 1
      B.ps_model <- tryCatch(
        {
          glm(formula = ps.formula, data = B_sample, family = binomial, control = list(maxit = 30))
        },
        error = function(e) {
          if (grepl("system is exactly singular", e$message)) {
            print(paste("Singular matrix encountered in bootstrap sample", b))
            return(NULL)
          } else {
            stop(e) 
          }
        },
        warning = function(w) {
          handle_warning(w, b)
          return(NULL)
        }
      )
      
      if (is.null(B.ps_model)) {
        next  
      }
      count.validBoostrap <- count.validBoostrap + 1
      #setTxtProgressBar(pb, count.validBoostrap)
      # Calculate the propensity score predictions
      B_sample$B.e_X <- predict(B.ps_model, type = "response")
      ### Fixed and Est for WATE
      B.Point_Fixed <- Calc.Point_WATE(sample = B_sample, ps.estimate = B_sample$e_X, weight = weight)
      B.Point_Est <- Calc.Point_WATE(sample = B_sample, ps.estimate = B_sample$B.e_X, weight = weight)
      ### Fixed and Est for Aug-WATE (TRUE)
      B.Point_Fixed.AugT <- Calc.Point_AugWATE(sample = B_sample, ps.estimate = B_sample$e_X, 
                                          weight = weight, out.formula = out.TRUE.formula)
      B.Point_Est.AugT <- Calc.Point_AugWATE(sample = B_sample, ps.estimate = B_sample$B.e_X, 
                                          weight = weight, out.formula = out.TRUE.formula)
      ### Fixed and Est for Aug-WATE (FALSE)
      B.Point_Fixed.AugF <- Calc.Point_AugWATE(sample = B_sample, ps.estimate = B_sample$e_X,
                                          weight = weight, out.formula = out.FALSE.formula)
      B.Point_Est.AugF <- Calc.Point_AugWATE(sample = B_sample, ps.estimate = B_sample$B.e_X,
                                          weight = weight, out.formula = out.FALSE.formula)
      B.Point[count.validBoostrap,] <- c(B.Point_Fixed, B.Point_Est, 
                                         B.Point_Fixed.AugT, B.Point_Est.AugT,
                                         B.Point_Fixed.AugF, B.Point_Est.AugF)
    }
    #close(pb)
    count.Summary <- data.frame(Prop.Quasi_Separa = count.quasi_separa/b,
                                Prop.Non_Convg = count.non_convg/b,
                                Prop.validBoostrap = count.validBoostrap/b)
    return(list(B.Point_results = B.Point, 
                Summary.Warning = count.Summary))
}

# Std.Bootstrap(B.sample = 1000, ps.formula, p, sample, weight = "IPW", out.TRUE.formula, out.FALSE.formula)

Summary.Bootstrap <- function(B.Point_results, ps.formula, p, sample, alpha, weight, out.TRUE.formula, out.FALSE.formula) {
  ps_model <- glm(ps.formula, data = sample, family = binomial)
  ps.estimate <- predict(ps_model, type = "response")
  # Point Estimate Summary based on Random Sample
  Point.Est_WATE <- Calc.Point_WATE(sample, ps.estimate, weight)
  Point.Est_AugT.WATE <- Calc.Point_AugWATE(sample, ps.estimate, weight, out.TRUE.formula)
  Point.Est_AugF.WATE <- Calc.Point_AugWATE(sample, ps.estimate, weight, out.FALSE.formula)
  ########### SE Estimation based on Bootstrap Sample ##########
  SE_Summary <- as.data.frame(t(apply(B.Point_results, MARGIN = 2, sd)))
  #################### CI Estimation #########################
  z_value <- qnorm(1 - alpha / 2)
  colnames <- names(SE_Summary)
  # Mapping point estimates to column names
  Point.Ests <- setNames(c(Point.Est_WATE, Point.Est_AugT.WATE, Point.Est_AugF.WATE),
                         c("WATE", "AugT", "AugF"))
  CIs <- data.frame(Method = character(), Lower = numeric(), Upper = numeric(), stringsAsFactors = FALSE)
  for (col in colnames) {
    Point.Est <- Point.Est_WATE  # Default to Point.Est_WATE
    if (grepl("AugT", col)) {
      Point.Est <- Point.Ests["AugT"]
    } else if (grepl("AugF", col)) {
      Point.Est <- Point.Ests["AugF"]
    }
    SE_value <- SE_Summary[[col]][1]
    # Basic Bootstrap or Parametric Bootstrap
    lower_bound <- Point.Est - z_value * SE_value
    upper_bound <- Point.Est + z_value * SE_value
    CIs <- rbind(CIs, data.frame(Method = paste0(col, "_Par"), Lower = lower_bound, Upper = upper_bound))
    # Empirical Bootstrap based on Percentile
    lower_bound <- quantile(B.Point_results[[col]], probs = alpha / 2)
    upper_bound <- quantile(B.Point_results[[col]], probs = 1 - alpha / 2)
    CIs <- rbind(CIs, data.frame(Method = paste0(col, "_Pct"), Lower = lower_bound, Upper = upper_bound))
    # Double Bootstrap Method
    lower_bound <- 2 * Point.Est - quantile(B.Point_results[[col]], probs = 1 - alpha / 2)
    upper_bound <- 2 * Point.Est - quantile(B.Point_results[[col]], probs = alpha / 2)
    CIs <- rbind(CIs, data.frame(Method = paste0(col, "_Double"), Lower = lower_bound, Upper = upper_bound))
    # Bias-Corrected and Accelerated Interval (BCa)
    bca_Results <- bca(B.Point_results[[col]], conf.level = 1 - alpha)
    lower_bound <- bca_Results[1]
    upper_bound <- bca_Results[2]
    CIs <- rbind(CIs, data.frame(Method = paste0(col, "_BCa"), Lower = lower_bound, Upper = upper_bound))
  }
  
  row.names(CIs) <- NULL
  return(list(SEs = SE_Summary, CIs = CIs))
}
