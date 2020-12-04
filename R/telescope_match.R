#' Perform telescope matching to estimate the controlled
#' direct effect of a binary treatment net the effect of a binary mediator.
#'
#' @param outcome character indicating the name of the
#' variable in \code{data} being used as the outcome. This variable must appear
#' on the left-hand side of both \code{s1.formula} and \code{s2.formula}.
#' @param treatment character indicating the name of the
#' variable in \code{data} being used as the treatment indicator. This variable
#' must be a binary integer (either 0 or 1). This variable must appear
#' on the right-hand side of  \code{s1.formula}. It should also appear in \code{s2.formula}.
#' @param mediator character indicating the name of the
#' variable in \code{data} being used as the mediator indicator. This variable
#' must be a binary integer (either 0 or 1). This method returns controlled direct effects
#' of setting the mediator to 0. Recode the indicator if the controlled direct effect of setting the 
#' mediator to 1 is of interest. This variable cannot appear in either
#' \code{s1.formula} or \code{s2.formula}.
#' @param s1.formula  A formula object denoting the stage 1 formula. The outcome should appear on the left-hand side.
#' Treatment, pre- and post-treatment covariates should appear on the right-hand side. The mediator is omitted as the model
#' is fit only within the subset of observations with \code{mediator} equal to 0.
#' @param s2.formula A formula object denoting the stage 2 formula estimating the ACDE of treatment using the demediated
#' matches from the first stage. The outcome should appear on the left-hand side.
#' Treatment and pre-treatment covariates should appear on the right-hand side.
#' @param data A dataframe containing columns referenced by
#' \code{outcome}, \code{treatment} and \code{mediator} along with any variables
#' referenced in \code{s1.formula} and \code{s2.formula}.
#' @param caliper A scalar denoting the caliper to be used in matching in the treatment stage (calipers cannot be used for matching
#' on the mediator). Observations outside of the caliper are dropped. Calipers are specified in standard deviations of the covariates.
#' NULL by default (no caliper).
#' @param L Number of matches to use for each unit. Must be a numeric vector of eitehr length 1 or 2. If length 1, L 
#' sets the number of matches used in both the first stage (matching on mediator) and in the second stage (matching on treatment).
#' If length 2, the first element sets the number of matches used in the first stage (matching on mediator) and the second element
#' sets the number of matches used in the second stage (matching on treatment) Default is 5.
#' @param boot logical indicating whether to conduct inference using the weighted bootstrap 
#' for matching estimators extended from Otsu and Rai (2017) (\code{TRUE}) or the asymptotic variance 
#' estimator from Blackwell and Strezhnev (2019) (\code{FALSE}). Defaults to \code{FALSE}.
#' @param nBoot If \code{boot} is \code{TRUE}, number of bootstrap iterations to use. Default is \code{5000}.
#' @param ci percent level of confidence interval to return. If \code{boot} is \code{FALSE}, returns symmetric asymptotic 
#' interval centered on the estimated treatment effect. If \code{boot} is \code{TRUE} returns the \code{(100 - ci)/2} and
#' \code{100 - (100 - ci)/2} percentiles of the bootstrap distribution. Must be in the interval \code{(0, 100)}. Defaults to 95.
#' @param verbose logical indicating whether to display progress information. Default is \code{TRUE}.
#' 
#' @details The \code{telescope_match} function implements the two-stage
#' "telescope matching" procedure developed by Blackwell and Strezhnev 
#' (2019).
#' 
#'  The procedure first estimates a demediated outcome using a combination
#'  of matching and a regression bias-correction. The first stage formula 
#'  \code{s1.formula} specifies the pre- and post-treatment covariates to be used
#'  in matching along with the specification for the bias-correction regression. In this stage,
#'  all units with M = 1 to units with M = 0 with identical treatment values and 
#'  similar pre- and post-treatment covariates. The potential outcomes under M = 0
#'  are imputed using the average of the matches + the bias correction from the regression
#'  model. The second stage estimates the controlled direct effect of treatment
#'  on this demediated outcome using a similar matching/bias-correction procedure with
#'  the formula specified in \code{s2.formula} indicating the pre-treatment covariates used
#'  along with the treatment.
#'  
#'  Matching is performed using the \code{Match()} routine from the \code{Matching} package.
#'  By default, matching is L-to-1 nearest neighbor with replacement using Mahalanobis distance.
#'  
#'
#' See the references below for more details.
#'
#' @return Returns an object of \code{class} \code{tmatch}. Contains the following components
#' \itemize{
#' \item estimate: Estimated ACDE fixing M=0
#' \item std.err: Estimated asymptotic standard error. \code{NULL} if \code{boot} is \code{TRUE}
#' \item boot.dist: Bootstrap distribution of \code{estimate}. \code{NULL} if \code{boot} is \code{FALSE}
#' \item conf.low: Lower bound of \code{ci} confidence interval for the estimate
#' \item conf.high: Upper bound of \code{ci} confidence interval for the estimate
#' \item ci.level: Level of the confidence interval
#' \item outcome: Name of outcome variable
#' \item treatment: Name of treatment variable
#' \item mediator: Name of mediator variable
#' \item pre.treatment: Vector of names of pre-treatment confounders (appear in both stage 1 and 2)
#' \item post.treatment: Vector of names of post-treatment confounders (appear only in stage 1)
#' \item s1.formula: Stage 1 bias-correction regression formula (pre-/post-treatment covariates)
#' \item s2.formula: Stage 2 bias-correction regression formula (pre-treatment covariates)
#' \item outcome.vec: Vector of outcomes used in estimation
#' \item treatment.vec: Vector of treatment indicators used in estimation
#' \item mediator.vec: Vector of mediator indicators used in estimation
#' \item L_m: Number of matches found for each unit in the first stage mediator matching procedure
#' \item L_a: Number of matches found for each unit in the second stage mediator matching procedure
#' \item KLm: Number of times unit is used as a match in the first stage mediator matching procedure
#' \item KLa: Number of times unit is used as a match in the second stage treatment matching procedure
#' \item N: Number of observations
#' \item N_summary: Number of observations in each treatment/mediator combination.
#' \item caliper: Caliper (if any) used in the treatment stage matching to drop distant observations.
#' }

#' @references Blackwell, Matthew, and Strezhnev, Anton (2020)
#' "Telescope Matching: Reducing Model Dependence 
#' in the Estimation of Direct Effects." Working Paper.
#' 
#' @examples
#' data(jobcorps)
#' 
#' ## Split male/female
#' jobcorps_female <- jobcorps %>% filter(female == 1)
#' 
#' ## Telescope matching formula - First stage (X and Z)
#' tm_stage1 <- exhealth30 ~ treat*(schobef + trainyrbef + jobeverbef + jobyrbef + health012 + health0mis +  pe_prb0 + 
#'                                    everalc + alc12 + everilldrugs + age_cat +  eduhigh + rwhite + everarr + hhsize + hhsizemis +  hhinc12 + hhinc8 + fdstamp +
#'                                    welf1 + welf2 + publicass + emplq4 + emplq4full + pemplq4 + pemplq4mis + vocq4 + vocq4mis + 
#'                                    health1212 + health123 + pe_prb12 + pe_prb12mis  + 
#'                                    narry1 + numkidhhf1zero + numkidhhf1onetwo + pubhse12 + h_ins12a + h_ins12amis)
#' 
#' ## Telescope matching formula - second stage (X)
#' tm_stage2 <- exhealth30 ~ treat*(schobef + trainyrbef + jobeverbef + jobyrbef + health012 + health0mis +  pe_prb0 + 
#'                                    everalc + alc12 + everilldrugs + age_cat +  eduhigh +  rwhite + everarr + hhsize + hhsizemis + hhinc12 + hhinc8 + fdstamp +
#'                                    welf1 + welf2 + publicass)
#' 
#' 
#' ### Estimate ACDE for women holding employment at 0
#' telescopeMatch.result.0 <-  telescope_match(outcome = "exhealth30", treatment = "treat", mediator = "work2year2q", 
#'                                            s1.formula = tm_stage1, 
#'                                            s2.formula = tm_stage2, data=jobcorps_female, L=3, boot=F, nBoot=1000, verbose=T, ci=95)
#' 
#' @export
#' @importFrom Matching Match
#'
telescope_match <- function(outcome, treatment, mediator, s1.formula, s2.formula, data,  caliper = NULL, L=5, 
                            boot = F, nBoot=5000, ci = 95, verbose=T){

  ########################
  ### Quick pre-processing
  ########################
  
  ### If length of L is 1, L_a and L_m are the same 
  if (!is.numeric(L)|!is.vector(L)){
    stop("Error: L must be a numeric vector of length 1 or 2")
  }else if (length(L) > 2|length(L) < 1){
    stop("Error: L must be a numeric vector of length 1 or 2")
  }else if (length(L) == 2){
    L_m <- L[1]
    L_a <- L[2]
  }else if (length(L) == 1){
    L_m <- L
    L_a <- L
  }

  ### Force data to be a data frame. Tidy tibbles crash the existing code
  data <- as.data.frame(data)
  
  ### Get first order terms from s1.formula and s2.formula
  s1.terms <- rownames(attributes(terms(s1.formula))$factors)
  s2.terms <- rownames(attributes(terms(s2.formula))$factors)
  
  ### Number of observations
  N <- nrow(data)

  ####################
  #### Sanity checks 
  ####################

  ### Check -2 - Is the ci within (0, 100).
  if (ci <= 0 | ci >= 100){
    warning("Warning: 'ci' must be within the interval (0, 100). Defaulting to 95% confidence intervals.")
    ci <- 95
  }
  ci.alpha <- 1 - ci/100
  
  ### Check -1 - Outcome, Treatment, Mediator are in data?
  if (!is.character(outcome)| !is.character(treatment)| !is.character(mediator)){
    stop("Error: 'outcome', 'treatment', and 'mediator' must be characters.")
  }else if ((length(outcome) != 1)|(length(treatment) != 1)|(length(mediator) != 1)){
    stop("Error: Cannot have a vector greater than length 1 for 'outcome', 'treatment', and 'mediator'.")
  }else if (!(outcome %in% colnames(data))|!(treatment %in% colnames(data))|!(mediator %in% colnames(data))){
    stop("Error: 'outcome', 'treatment', or 'mediator' not in 'data'.")
  }
  
  ### Check 0 - Outcome in both s1.terms and s2.terms
  if (!(outcome %in% s1.terms)|!(outcome %in% s2.terms)){
    stop("Error: 'outcome' not in s1.formula or s2.formula.")
  }
  
  ## Check 1 - Treatment is in s1.terms and s2.terms
  if (!(treatment %in% s1.terms)|!(treatment %in% s2.terms)){
    stop("Error: 'treatment' is not in s1.formula or s2.formula")
  }
  
  ## Check 2 - Mediator isn't in either s1.formula, s2.formula
  if ((mediator %in% s1.terms)|(mediator %in% s2.terms)){
    stop("Error: 'mediator' in s1.formula or s2.formula.")
  }
  
  ## Check 3 - Treatment is binary 0-1
  if (!(is.numeric(data[[treatment]]))){
    stop("Error: 'treatment' variable not numeric.")
  }else if(!isTRUE(all.equal(unique(data[[treatment]])[order(unique(data[[treatment]]))], c(0,1)))){
    stop("Error: 'treatment' must have only two levels: 0,1")
  }
  
  ## Check 4 - Mediator is binary 0-1
  if (!(is.numeric(data[[mediator]]))){
    stop("Error: 'mediator' variable not numeric.")
  }else if(!isTRUE(all.equal(unique(data[[mediator]])[order(unique(data[[mediator]]))], c(0,1)))){
    stop("Error: 'mediator' must have only two levels: 0,1")
  }
  
  ### Pre-treatment covariates are s2.terms
  pre.treatment <- s2.terms[!(s2.terms %in% c(outcome, treatment, mediator))]
  s1.covariates <- s1.terms[!(s1.terms %in% c(outcome, treatment, mediator))]
  
  ### Check 5 - Are there any s2.terms that aren't in s1?
  if (any(!(pre.treatment %in% s1.covariates))){
    stop("Error: some covariates in second stage 's2.formula' not in the first stage 's1.formula'")
  }  
  
  ######################
  ### Pre-estimation set-up
  ######################
  
  ### Post-treatment covariates are those in first stage but not in second stage
  post.treatment <- s1.covariates[!(s1.covariates %in% pre.treatment)]
  
  ## Combine all covariates into a single vector
  all.covariates <- c(pre.treatment, post.treatment)
  
  ## Diagnostic, print implied pre-treatment and post-treatment covariates
  if (verbose){
    cat("Telescope matching setup:\n")
    cat(paste("Outcome: ", outcome, "\n", sep=""))
    cat(paste("Treatment: ", treatment, "\n", sep=""))
    cat(paste("Mediator: ", mediator, "\n", sep=""))
    cat(paste("Pre-treatment Covariates: ", paste(pre.treatment, collapse = ", "), "\n", sep=""))
    cat(paste("Post-treatment Covariates: ", paste(post.treatment, collapse = ", "), "\n", sep=""))
    cat(paste("Number of matches in first stage (mediator): ", L_m, "\n", sep=""))
    cat(paste("Number of matches in second stage (treatment): ", L_a, "\n", sep=""))
    cat("\n")
  }
  
  ######################
  ### First-stage: Estimate Y(a, 0) conditional on pre-/post-treatment covariates
  ######################
  
  ### First stage matching - inexact on all pre.treatment/post.treatment covariates, exact on A"
  tm.first <- Match(Y = data[[outcome]], Tr = data[[mediator]], X = data[,c(all.covariates, treatment)], 
                    exact = c(rep(FALSE, length(all.covariates)), TRUE), M=L_m, BiasAdjust=FALSE, estimand = "ATT", ties=F)

  ### Summarize input - First Stage
  if (verbose){
    cat("First-stage matching: Mediator on pre-treatment, post-treatment\n")
    cat("Number of observations with 'mediator' = 0 matched to each observation with 'mediator' = 1: ", L_m, "\n", sep="")
    cat("Number of observations with 'mediator' = 1: ", sum(data[[mediator]]), "\n", sep = "" )
    cat("Number of observations with 'mediator' = 0: ", sum(1 - data[[mediator]]), "\n", sep = "")
    cat("\n")
  }
  
  ### Summarize sample sizes with mediator/treatment
  n_summary = data.frame(c(0,1,0,1), c(0,0,1,1), as.vector(table(data[[treatment]],data[[mediator]])))
  colnames(n_summary) <- c(treatment, mediator, "N")
  
  ### Count of # of matches for each M=0 unit in the first stage (K_L^m)
  KLm <- table(tm.first$index.control)
  ### Which M=1 units did each M=0 match to?
  first_stage_matchedTo <- tapply(tm.first$index.treated, tm.first$index.control, function(x) x)
  ## Save K_L^m in data
  data$KLm <- 0
  data[as.numeric(names(KLm)), "KLm"] <- KLm
  
  ### Regression imputation in the first stage
  s1.reg.0 <- lm(s1.formula, data=data[data[[mediator]]==0,]) ## Fit a regression in the set with mediator = 0
  
  ## Predicted value of regression for each X under control
  data$pred.Y.m0 <- NA
  data$pred.Y.m0 <- predict(s1.reg.0, newdata = data[,c(treatment, pre.treatment, post.treatment)]) 
  
  #### Average of matched values for each unit
  data$y.m0.m.imp <- data[[outcome]] ## Matching imputation
  data[unique(tm.first$index.treated), "y.m0.m.imp"] <- tapply(tm.first$index.control, tm.first$index.treated, function(x) mean(data[x, outcome]))
  #### Imputed value under M=0 using regression-imputed value of the matches
  data$y.m0.r.imp <- data$pred.Y.m0 ## Regression imputation
  data[unique(tm.first$index.treated), "y.m0.r.imp"] <- tapply(tm.first$index.control, tm.first$index.treated, function(x) mean(data[x, "pred.Y.m0"]))
  
  ### Impute the outcome under M=0 for units with M=1
  data$Ytilde <- data[[outcome]] ## Bias-corrected matching imputation for the M=1 units
  ### Average of imputations + (Regression prediction for X_i - regression prediction for all imputed units)
  data$Ytilde[data[[mediator]] == 1] <- data[data[[mediator]] == 1,]$y.m0.m.imp + data[data[[mediator]] == 1,]$pred.Y.m0 - data[data[[mediator]] == 1,]$y.m0.r.imp 
  
  ################################################
  ### Second stage - Estimate E[Y(1,0) - Y(0,0)] - conditional on pre-treatment covariates
  ################################################
  
  ### Summarize input - Second stage
  if (verbose){
    cat("Second-stage matching: Treatment on pre-treatment\n")
    cat("Number of observations matched to each unit with opposite treatment: ", L_a,  "\n", sep="")
    cat("Number of observations with 'treatment' = 1: ", sum(data[[treatment]]), "\n", sep = "" )
    cat("Number of observations with 'treatment' = 0: ", sum(1 - data[[treatment]]), "\n", sep = "")
    cat("\n")
  }
  
  ### Fix the outcome in s2.formula (replace with Ytilde)
  s2.form.character <- as.character(s2.formula)
  s2.form.character[2] <- "Ytilde"
  s2.formula.fixed <- formula(paste(s2.form.character[2], " ~ ", s2.form.character[3], sep=""))
  
  ##### Second stage regression on blipped-down outcome
  ## Ytilde = Y_i(a, 0)
  s2.reg <- lm(s2.formula.fixed, data = data) ## Fully interacted regression of Ytilde with treatment A
  data$pred.Y.a1 <- NA
  newdata.treat <- data[,c(treatment, pre.treatment)]
  newdata.treat[[treatment]] <- 1
  data$pred.Y.a1 <- predict(s2.reg, newdata = newdata.treat) ## Predicted outcome under treatment
  newdata.ctrl <- data[,c(treatment, pre.treatment)]
  newdata.ctrl[[treatment]] <- 0
  data$pred.Y.a0 <- NA
  data$pred.Y.a0 <-  predict(s2.reg, newdata = newdata.ctrl) ## Predicted outcome under control
  
  data$pred.Y.A <- data[[treatment]] * data$pred.Y.a1 + (1 - data[[treatment]]) * data$pred.Y.a0 ### Predict factual
  data$pred.Y.1.A <- data[[treatment]] * data$pred.Y.a0 + (1 - data[[treatment]]) * data$pred.Y.a1 ## Predict counterfactual
  
  ### Match controls to treated
  tm.second.a1 <- Match(Y = data$Ytilde, Tr = data[[treatment]], X = data[,c(pre.treatment)], estimand = "ATT", caliper = caliper, M = L_a, ties=F)
  KLa0 <- table(tm.second.a1$index.control) ## Count of matched controls - stage 2
  ### Match treateds to control
  tm.second.a0 <- Match(Y = data$Ytilde, Tr = data[[treatment]], X = data[,c(pre.treatment)], estimand = "ATC", caliper = caliper, M = L_a, ties=F)
  KLa1 <- table(tm.second.a0$index.treated) ## Count of matched treated - stage 2
  
  ## Total match counts - stage 2
  data$KLa <- 0
  data[as.numeric(names(KLa0)), "KLa"] <- KLa0
  data[as.numeric(names(KLa1)), "KLa"] <- KLa1
  
  SLmatch <- lapply(first_stage_matchedTo, function(x) sum(data$KLa[x])) ## Sum of KLa of matched units
  data$SLm <- 0
  data$SLm[as.integer(names(SLmatch))] <- unlist(SLmatch)
  
  ## Imputed Y_i(0,0) of matches - regression
  data$Yhat00.r <- data$pred.Y.a0
  data[unique(tm.second.a1$index.treated), "Yhat00.r"] <- tapply(tm.second.a1$index.control, tm.second.a1$index.treated, function(x) mean(data[x, "pred.Y.a0"]))
  ## Imputed Y_i(1,0) of matches - regression
  data$Yhat10.r <- data$pred.Y.a1
  data[unique(tm.second.a0$index.control), "Yhat10.r"] <- tapply(tm.second.a0$index.treated, tm.second.a0$index.control, function(x) mean(data[x, "pred.Y.a1"]))
  
  ## regression imputations of the CEF under 1-A_i
  ## (1/L) \sum_{j \in J^a(i)} mu_{1-A_i,0}(X_j, 1-A_i)
  data$pred.Y.1.A.r <- data[[treatment]] * data$Yhat00.r + (1 - data[[treatment]]) * data$Yhat10.r
  
  ## Linearized form for bootstrapping
  sm.part <- (1 - data[[mediator]]) * (1 + data$KLa/L_a + data$KLm/L_m + data$SLm/(L_a*L_m)) * data[[outcome]]
  bm.part <- ((1 - data[[mediator]])*(data$KLm/L_a + data$SLm/(L_a*L_m)) - data[[mediator]]*(1 + data$KLa/L_a))*data$pred.Y.m0
  ba.part <- data$pred.Y.1.A + (data$KLa/L_a)*data$pred.Y.A
  
  data$tau.i <- (2 * data[[treatment]] - 1) * (sm.part - bm.part - ba.part)
  
  ## Point estimate
  tau <- mean(data$tau.i)

  ########################################
  ## Variance estimation
  ########################################
  
  if (boot == F){
    
    ###############
    ### Asymptotic variance
    ###################
    
    ### Calculate weights
    data$ww <- (1 - data[[mediator]]) * (1 + data$KLa/L_a + data$KLm/L_m + data$SLm/(L_a*L_m))
    
    ### Number of units with M=0
    N0 <- sum(1 - data[[mediator]])
    ### Parameters in the first-stage regression
    P1 <- length(coef(s1.reg.0))
    ### Parameters in the second-stage regression
    P2 <- length(coef(s2.reg))
    
    ### Variance component 1
    data$em.var <- (N0/(N0 - P1)) * (data[[outcome]] - data$pred.Y.m0)^2
    
    ### Variance component 2
    data$ea.var <- (N/(N-P2)) * (data$pred.Y.m0 - data$pred.Y.A)^2
    
    ### Variance component 3
    tau.var <- mean((data$pred.Y.a1 - data$pred.Y.a0 - tau)^2)
    
    ### combine all three components
    se.est2 <- sqrt(tau.var/N + (mean((1-data[[mediator]]) * data$ww^2 * data$em.var) + mean((1+data$KLa/L_a)^2*data$ea.var))/N)
    
    ### No bootstrap
    Tstar <- NULL
    
    ### CI bounds
    ci.low = tau - abs(qnorm(ci.alpha/2))*se.est2
    ci.high = tau + abs(qnorm(ci.alpha/2))*se.est2
    
  }else{
    
    ################
    ### Wild Bootstrap
    ################
    
    ## De-mean
    tau.norm <- data$tau.i - tau
    
    ## Bootstrap iterations
    W.bern <- sapply(1:nBoot, function(x) rbinom(N, 1, prob = (sqrt(5) - 1)/(2*sqrt(5))))
    Wstar <- (((sqrt(5) + 1)/2)*W.bern + (-(sqrt(5) - 1)/2)*(1-W.bern))/N
    
    ## Apply bootstrap weights to each "observation"
    Tstar <- (t(Wstar) %*% tau.norm)
    
    ## Add in tau
    Tstar <- Tstar + tau
    
    ### Asypmtotic estimate is null
    se.est2 <- NULL
    
    ### Get quantiles for the CI
    ci.low = quantile(Tstar, ci.alpha/2)
    ci.high = quantile(Tstar, 1 - ci.alpha/2)
  }
  
  ### Return output
  output <- list(outcome = outcome, treatment = treatment, mediator = mediator, s1.formula = s1.formula, s2.formula = s2.formula,
                 N = N, L_m = L_m, L_a = L_a, N_summary = n_summary, 
                 estimate=tau, std.err = se.est2, boot.dist=as.vector(Tstar), KLm = data$KLm,
                 KLa = data$KLa, outcome.vec = data[[outcome]], treatment.vec = data[[treatment]], mediator.vec = data[[mediator]],
                 pre.treatment = pre.treatment, post.treatment = post.treatment,
                 conf.low = ci.low, conf.high = ci.high, ci.level = ci)
  
  class(output) <- "tmatch"
  return(output)
  
}

#' Summarize telescope match objects
#' 
#' @details \code{summary} method for \code{tmatch} objects returned by
#' \code{telescope_match}
#'
#' @param object an object of class \code{tmatch} -- results from a call to \code{telescope_match} 
#' 
#' @details Returns a summary data frame containing the estimate, standard error and confidence
#' interval from the `telescope_match` object.
#'
#' @return Returns an object of \code{class} \code{summary.tmatch}. Contains the following components
#' \itemize{
#' \item outcome: Name of outcome variable
#' \item treatment: Name of treatment variable
#' \item mediator: Name of mediator variable
#' \item pre.treatment: Vector of names of pre-treatment confounders (appear in both stage 1 and 2)
#' \item post.treatment: Vector of names of post-treatment confounders (appear only in stage 1)
#' \item sizes: Number of observations in each treatment-mediator combination
#' \item L_m: Number of units matched in the mediator (first) stage
#' \item L_a: Number of units matched in the treatment (second) stage
#' \item estimate: Point estimate of the ACDE
#' \item se.type: Character indicating the type of standard error estimator used in the 
#' telescope_matching routine, either 'Asymptotic" or 'Bootstrap'
#' \item std.err: Estimated standard error of the ACD
#' \item ci.level: Confidence interval level
#' \item conf.low: Lower bound of the `ci.level` asymptotically normal confidence interval
#' \item conf.high: Upper bound of the `ci.level` asymptotically normal confidence interval
#' }
#' @export
summary.tmatch <- function(object){
  
  summary_obj <- NULL
  
  ### Names of treatments/outcome
  summary_obj$outcome <- object$outcome
  summary_obj$treatment <- object$treatment
  summary_obj$mediator <- object$mediator
  
  ### Covariates
  summary_obj$pre.treatment <- object$pre.treatment
  summary_obj$post.treatment <- object$post.treatment
  
  ### Sample size
  summary_obj$sizes <- object$N_summary
  
  ### Number of units matched in each stage
  summary_obj$L_m <- object$L_m
  summary_obj$L_a <- object$L_a
  
  ### Estimates
  if (!is.null(object$std.err)){
    summary_obj$se.type = "Asymptotic"
    summary_obj$std.err = object$std.err
  }else{
    summary_obj$se.type = "Bootstrap"
    summary_obj$std.err = sd(object$boot.dist)
  }
  
  summary_obj$estimate = object$estimate
  summary_obj$ci.level = object$ci.level
  summary_obj$conf.low = object$conf.low
  summary_obj$conf.high = object$conf.high
  
  class(summary_obj) <- "summary.tmatch"
  
  return(summary_obj)
  
}

print.summary.tmatch <- function(object, digits = max(3, getOption("digits") - 3)){
  cat("Telescope matching results\n")
  cat("----------------------------\n")
  cat("Variables:\n")
  cat(paste("Outcome: ", object$outcome, "\n", sep=""))
  cat(paste("Treatment: ", object$treatment, "\n", sep=""))
  cat(paste("Mediator: ", object$mediator, "\n", sep=""))
  cat("----------------------------\n")
  cat(paste("Pre-treatment Covariates: ", paste(object$pre.treatment, collapse = ", "), "\n", sep=""))
  cat(paste("Post-treatment Covariates: ", paste(object$post.treatment, collapse = ", "), "\n", sep=""))
  cat("----------------------------\n")
  cat(paste("Number of units matched to each observation in first stage (mediator): ", object$L_m, "\n", sep=""))
  cat(paste("Number of units matched to each observation in second stage (treatment): ", object$L_a, "\n", sep=""))
  cat("----------------------------\n")
  cat("Results:\n")
  if (object$se.type == "Asymptotic"){
    cat(paste("Estimated ACDE: ", format(object$estimate, digits=digits), "\n", sep=""))
    cat(paste("Asymptotic Standard Error: ", format(object$std.err, digits=digits), "\n", sep=""))
    cat(paste("Asymptotic ", format(object$ci.level, digits=digits), "% confidence interval: [", format(object$conf.low, digits=digits), ", ", format(object$conf.high, digits=digits), "]\n", sep=""))
  }else if (object$se.type == "Bootstrap"){
    cat(paste("Estimated ACDE: ", format(object$estimate, digits=digits), "\n", sep=""))
    cat(paste("Bootstrap Standard Error: ", format(object$std.err, digits=digits), "\n", sep=""))
    cat(paste("Bootstrap ", format(object$ci.level, digits=digits), "% confidence interval (percentile): [", format(object$conf.low, digits=digits), ", ", format(object$conf.low, digits=digits), "]\n", sep=""))
  }else{
    cat("Something went wrong with the summary - SEs not Asymptotic or Bootstrap\n")
  }
  cat("----------------------------\n")
  cat("Sample sizes (N):\n")
  cat(paste("Treatment = 0, Mediator = 0: ", object$sizes$N[object$sizes[[object$treatment]] == 0&object$sizes[[object$mediator]] == 0], "\n", sep=""))
  cat(paste("Treatment = 0, Mediator = 1: ", object$sizes$N[object$sizes[[object$treatment]] == 0&object$sizes[[object$mediator]] == 1], "\n", sep=""))
  cat(paste("Treatment = 1, Mediator = 0: ", object$sizes$N[object$sizes[[object$treatment]] == 1&object$sizes[[object$mediator]] == 0], "\n", sep=""))
  cat(paste("Treatment = 1, Mediator = 1: ", object$sizes$N[object$sizes[[object$treatment]] == 1&object$sizes[[object$mediator]] == 1], "\n", sep=""))
  cat("----------------------------\n")

}

#' Balance diagnostics for Telescope Match objects
#' 
#' @details Provides matching balance diagnostics for \code{tmatch} objects returned by
#' \code{telescope_match}
#'
#' @param object an object of class \code{tmatch} -- results from a call to \code{telescope_match} 
#' @param vars a formula object containing either the treatment or the mediator as the dependent variable (which denotes whether first-stage or
#' second-stage balance diagnostics are returned) and the covariates for which balance diagnostics are requested as the independent variables. Each covariate
#' or function of covariates (e.g. higher-order polynomials or interactions) should be separated by a +.
#' @param data the data frame used in the call to \code{telescope_match}
#'
#' @return Returns a data frame with the following columns.
#' \itemize{
#' \item variable: Name of covariate
#' \item Before_M0/Before_A0: Pre-matching average of the covariate in the mediator == 0 (if first stage balance) or 
#' treatment == 0 (if second stage balance) condition
#' \item Before_M1/Before_A1: Pre-matching average of the covariate in the mediator == 1 (if first stage balance) or 
#' treatment == 1 (if second stage balance) condition
#' \item After_M0/After_A0: Post-matching average of the covariate in the mediator == 0 (if first stage balance) or 
#' treatment == 0 (if second stage balance) condition
#' \item After_M1/After_A1: Post-matching average of the covariate in the mediator == 1 (if first stage balance) or 
#' treatment == 1 (if second stage balance) condition
#' \item SD: standard deviation of the outcome (pre-Matching)
#' \item Before_Diff: Pre-matching covariate difference between mediator arms (if first stage balance) or treatment arms
#' (if second stage balance).
#' \item Before_Std_Diff: Pre-matching standardized covariate difference between mediator arms (if first stage balance) or
#' treatment arms (if second stage balance), Equal to Before_Diff/SD.
#' \item After_Diff: Post--matching covariate difference between mediator arms (if first stage balance) or treatment arms
#' (if second stage balance).
#' \item After_Std_Diff: Post-matching standardized covariate difference between mediator arms (if first stage balance) or
#' treatment arms (if second stage balance), Equal to Before_Diff/SD.
#' }
#' 
#' @export


balance.tmatch <- function(object, vars, data){
  
  ##################
  ### Sanity checks
  ##################
  
  ### Is the class a 'tmatch'
  if (class(object) != "tmatch"){
    stop("Error: 'object' not of class 'tmatch'")
  }
  
  ### Does N of data match number of obs 
  if (nrow(data) != object$N){
    stop("Error: number of rows in data not equal to 'N' parameter in object")
  }
  
  ##########################
  ### Validating the formula
  if (object$treatment == as.character(vars[2])){
    ### Balance diagnostics for treatment
    get.balance = "treatment"
  }else if (object$mediator == as.character(vars[2])){
    ### Balance diagnostics for mediator
    get.balance = "mediator"
  }else{
    stop("Error: Left-hand side of 'vars' must be either the 'treatment' or the 'mediator' from 'object'")
  }
  
  ###########################
  ### Validating the data frame
  
  ## Do model.frame first to parse any functions
  covariate.frame = tryCatch({model.frame(vars, data[,-1])},error = function(e) { stop("Could not extract all variables in 'formula' from data")})
  #print(head(covariate.frame))
  ## Do model.matrix to get interactions
  covariate.frame = tryCatch({model.matrix(vars, covariate.frame)},error = function(e) { stop("Could not extract all variables in 'formula' from data")})

  ## Strip out the "(Intercept) column
  covariate.frame <- covariate.frame[,!(colnames(covariate.frame) %in% c("(Intercept)"))]
  
  ### Generate weights on each observation
  if (get.balance == "mediator"){
    ### If it's the mediator, all M=1 units get KLm+1, all M=0 get KLm
    obs.weights <- (object$KLm/object$L_m + 1)*(data[[object$mediator]] == 1) + (object$KLm/object$L_m)*(data[[object$mediator]] == 0)
    augmented.frame <- data.frame(cbind(data[[object$mediator]], obs.weights, covariate.frame))
    colnames(augmented.frame)[1] <- "mediator"
    ## Fix column names
    colnames(augmented.frame)[-c(1,2)] <- colnames(covariate.frame)

  }else if (get.balance == "treatment"){
    ### If it's the treatment, all A=1 units get KLa+1, all A=0 get KLa+1
    obs.weights <- (object$KLa/object$L_a + 1)
    augmented.frame <- data.frame(cbind(data[[object$treatment]], obs.weights, covariate.frame))
    colnames(augmented.frame)[1] <- "treatment"
    ## Fix column names
    colnames(augmented.frame)[-c(1,2)] <- colnames(covariate.frame)
    
  }

  ## Make a data-frame containing balance results
  if (get.balance == "mediator"){
    balance.frame <- data.frame(variable = colnames(augmented.frame)[-c(1,2)], Before_M0 = NA, Before_M1 = NA, After_M0 = NA, After_M1= NA, SD= NA)
  }else if (get.balance == "treatment"){
    balance.frame <- data.frame(variable = colnames(augmented.frame)[-c(1,2)], Before_A0 = NA, Before_A1 = NA, After_A0 = NA, After_A1= NA, SD = NA)
  }

  
  ### For each term in vars (aside from the treatment and the weights)
  for (variable in colnames(augmented.frame)[-c(1,2)]){

    ## Run the regression to get difference-in-means
    if (get.balance == "mediator"){
      ## Unweighted difference
      un.weight.reg <- lm(as.formula(paste(variable, "~", "mediator")), data=augmented.frame)
      balance.frame[balance.frame$variable == variable,]$Before_M0 = coef(un.weight.reg)[1] ## Intercept is M = 0
      balance.frame[balance.frame$variable == variable,]$Before_M1 = coef(un.weight.reg)[1] + coef(un.weight.reg)[2] ## Intercept + Beta1 is M = 1
      balance.frame[balance.frame$variable == variable,]$SD <- sd(augmented.frame[[variable]]) ## Get the marginal SD of the variable in the unweighted data.
        
      ## Weighted difference
      balance.reg <- lm(as.formula(paste(variable, "~", "mediator")), data=augmented.frame, weights=obs.weights)
      balance.frame[balance.frame$variable == variable,]$After_M0 = coef(balance.reg)[1] ## Intercept is M = 0
      balance.frame[balance.frame$variable == variable,]$After_M1 = coef(balance.reg)[1] + coef(balance.reg)[2] ## Intercept + Beta1 is M = 1
      
    }else if (get.balance == "treatment"){
      ## Unweighted difference
      un.weight.reg <- lm(as.formula(paste(variable, "~", "treatment")), data=augmented.frame)
      balance.frame[balance.frame$variable == variable,]$Before_A0 = coef(un.weight.reg)[1] ## Intercept is A = 0
      balance.frame[balance.frame$variable == variable,]$Before_A1 = coef(un.weight.reg)[1] + coef(un.weight.reg)[2] ## Intercept + Beta1 is A = 1
      balance.frame[balance.frame$variable == variable,]$SD <- sd(augmented.frame[[variable]])  ## Get the marginal SD of the variable in the unweighted data.
        
      ## Weighted difference
      balance.reg <- lm(as.formula(paste(variable, "~", "treatment")), data=augmented.frame, weights=obs.weights)
      balance.frame[balance.frame$variable == variable,]$After_A0 = coef(balance.reg)[1] ## Intercept is M = 0
      balance.frame[balance.frame$variable == variable,]$After_A1 = coef(balance.reg)[1] + coef(balance.reg)[2] ## Intercept + Beta1 is M = 1
      
    }
    
  }
  
  ## Save the standardized differences in means
  if (get.balance == "mediator"){
    balance.frame$Before_Diff <- balance.frame$Before_M1 - balance.frame$Before_M0
    balance.frame$Before_Std_Diff <- balance.frame$Before_Diff/balance.frame$SD
    
    balance.frame$After_Diff <- balance.frame$After_M1 - balance.frame$After_M0
    balance.frame$After_Std_Diff <- balance.frame$After_Diff/balance.frame$SD    
  }else if (get.balance == "treatment"){
    balance.frame$Before_Diff <- balance.frame$Before_A1 - balance.frame$Before_A0
    balance.frame$Before_Std_Diff <- balance.frame$Before_Diff/balance.frame$SD
    
    balance.frame$After_Diff <- balance.frame$After_A1 - balance.frame$After_A0
    balance.frame$After_Std_Diff <- balance.frame$After_Diff/balance.frame$SD    
  }
  
  
  ## Output
  return(balance.frame)
  
}

#' Histograms of matching weights
#' 
#' @details Provides histograms of the number of times each unit is used as a match given a  \code{tmatch} object returned by
#' \code{telescope_match}
#'
#' @param object an object of class \code{tmatch} -- results from a call to \code{telescope_match} 
#' @param stage a character vector equal to either 'mediator' or 'treatment'. If equal to 'mediator', returns a histogram
#' of matching weights for units with mediator = 0. If equal to 'treatment', returns a histogram of matching weights for
#' all units.
#' 
#' @returns Outputs a `plot()` object containing the histogram of match counts
#' 
#' @export 

plotDiag.tmatch <- function(object, stage = "mediator"){
  
  ##################
  ### Sanity checks
  ##################
  
  ### Is the class a 'tmatch'
  if (class(object) != "tmatch"){
    stop("Error: 'object' not of class 'tmatch'")
  }
  
  if (stage != "mediator" & stage != "treatment"){
    stop("Error: 'stage' must be either 'mediator' or 'treatment'.")
  }
  
  #######################
  ### Main output
  #######################
  
  ### If it's mediator, plot the histogram of K_Lm for mediator == 0
  if (stage == "mediator"){
    plot_title = paste("Matching weights for first stage (mediator) among M = 0")
    hist(object$KLm[object$mediator.vec == 0]/object$L_m, main=plot_title, xlab="Number of times unit is matched")
    abline(v = 1, col="red", lty=2, lwd=2)
  }else if (stage == "treatment"){
    plot_title = paste("Matching weights for second stage (treatment)")
    hist(object$KLa/object$L_a, main=plot_title, xlab="Number of times unit is matched")
    abline(v = 1, col="red", lty=2, lwd=2)
  }
  
  return()
  
}