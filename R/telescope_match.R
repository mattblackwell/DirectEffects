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
#' @param s2.formula A formula object denoting the stage 2 formula.
#' @param data A dataframe containing columns referenced by
#' \code{outcome}, \code{treatment} and \code{mediator} along with any variables
#' referenced in \code{s1.formula} and \code{s2.formula}.
#' @param subset A vector of logicals indicating which rows of \code{data} to keep.
#' @param verbose logical indicating whether to display progress information. Default is TRUE.
#' 
#' 
#' 
#' @details The \code{telescope_match} function implements the linear
#' sequential g-estimator developed by Vansteelandt (2009) with the
#' consistent variance estimator developed by Acharya, Blackwell, and
#' Sen (2016).
#'
#'  The formula specifies specifies the full first-stage model
#'   including treatment, baseline confounders, intermediate
#'   confounders, and the mediators. The user places \code{|} bars to
#'   separate out these different components of the model. For
#'   example, the formula should have the form \code{y ~ tr + x1 + x2
#'   | z1 + z2 | m1 + m2}. where \code{tr} is the name of the
#'   treatment variable, \code{x1} and \code{x2} are baseline
#'   covariates, \code{z1} and \code{z2} are intermediate covariates,
#'   and \code{m1} and \code{m2} are the names of the mediator
#'   variables. This last set of variables specify the 'blip-down' or
#'   'demediation' function that is used to remove the average effect
#'   of the mediator (possibly interacted) from the outcome to create
#'   the blipped-down outcome. This blipped-down outcome is the passed
#'   to a standard linear model with the covariates as specified for
#'   the direct effects model.
#'
#' See the references below for more details.
#'
#' @return Returns an object of \code{class} A \code{telescopeM}. Contains the following components
#' \itemize{
#' \item estimate: Estimated ACDE fixing M=0
#' \item se: Standard deviation of bootstrapped estimates
#' \item ci_low: Lower bound of `ci` confidence interval for the estimate - percentile of the bootstrapped distribution
#' \item ci_high: Upper bound of `ci` confidence interval for the estimate - percentile of the bootstrapped distribution
#' \item match_fs: Matching object returned from the first matching stage (on the mediator)
#' \item match.ss.c: Matching object returned from the second matching stage - matching treated units to control
#' \item match.ss.t: Matching object returned from the second matching stage - matching control units to treated
#' }
#' 
#' @references Blackwell, Matthew, and Strezhnev, Anton (2019)
#' "Telescope Matching: Reducing Model Dependence 
#' in the Estimation of Direct Effects." Working Paper.
#' 
#' @examples
#' data(jobcorps)
#'
#' @export
#' @importFrom Matching Match
#'
telescope_match <- function(outcome, treatment, mediator, s1.formula, s2.formula, data, subset=NULL,
                            boot = F, nBoot=5000, L=5, verbose=T){

  ########################
  ### Quick pre-processing
  ########################
  
  ### Force data to be a data frame. Tidy tibbles crash the existing code
  data <- as.data.frame(data)
  
  ### Get first order terms from s1.formula and s2.formula
  s1.terms <- rownames(attributes(terms(s1.formula))$factors)
  s2.terms <- rownames(attributes(terms(s2.formula))$factors)
  
  ### Number of observations
  N <- nrow(data)

  ### If subset isn't null
  if !is.null(subset){
    if (!is.logical(subset)){
      stop("Error: 'subset' must be a logical vector equal in length to the number of rows in 'data'")
    }else if(length(subset) != N){
      stop("Error: Length of 'subset' must equal the number of rows in 'data'")
    }else{
      data <- data[subset,]
    }
  }
  
  ####################
  #### Sanity checks 
  ####################
  
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
  
  ## Check 1 - Treatment is in s1.formula
  if (!(treatment %in% s1.terms)){
    stop("Error: 'treatment' is not in s1.formula")
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
  
  ### Are there any s2.terms that aren't in s1?
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
    cat("\n")
  }
  
  ######################
  ### First-stage: Estimate Y(a, 0) conditional on pre-/post-treatment covariates
  ######################
  
  ### First stage matching - inexact on all pre.treatment/post.treatment covariates, exact on A
  tm.first <- Match(Y = data[[outcome]], Tr = data[[mediator]], X = data[,c(all.covariates, treatment)], 
                    exact = c(rep(FALSE, length(all.covariates)), TRUE), M=L, BiasAdjust=FALSE, estimand = "ATT", ties=F)

  ### Summarize input - First Stage
  if (verbose){
    cat("First-stage matching: Mediator on pre-treatment, post-treatment\n")
    cat("Number of observations with 'mediator' = 1: ", sum(data[[mediator]]), "\n", sep = "" )
    cat("Number of observations with 'mediator' = 0: ", sum(1 - data[[mediator]]), "\n", sep = "")
    cat("\n")
  }
  
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
  tm.second.a1 <- Match(Y = data$Ytilde, Tr = data[[treatment]], X = data[,c(pre.treatment)], estimand = "ATT", M = L, ties=F)
  KLa0 <- table(tm.second.a1$index.control) ## Count of matched controls - stage 2
  ### Match treateds to control
  tm.second.a0 <- Match(Y = data$Ytilde, Tr = data[[treatment]], X = data[,c(pre.treatment)], estimand = "ATC", M = L, ties=F)
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
  sm.part <- (1 - data[[mediator]]) * (1 + data$KLa/L + data$KLm/L + data$SLm/(L^2)) * data[[outcome]]
  bm.part <- ((1 - data[[mediator]])*(data$KLm/L + data$SLm/(L^2)) - data[[mediator]]*(1 + data$KLa/L))*data$pred.Y.m0
  ba.part <- data$pred.Y.1.A + (data$KLa/L)*data$pred.Y.A
  
  data$tau.i <- (2 * data[[treatment]] - 1) * (sm.part - bm.part - ba.part)
  
  ## Point estimate
  tau <- mean(data$tau.i)
  
  ###############
  ### Asymptotic variance
  ###################
  
  if (boot == F){
    ### Calculate weights
    data$ww <- (1 - data[[mediator]]) * (1 + data$KLa/L + data$KLm/L + data$SLm/(L^2))
    
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
    se.est2 <- sqrt(tau.var/N + (mean((1-data[[mediator]]) * data$ww^2 * data$em.var) + mean((1+data$KLa/L)^2*data$ea.var))/N)
    
    ### No bootstrap
    Tstar <- NULL
    
  }else{
    
    ################
    ### Bootstrap
    ################
    
    ## De-mean
    tau.norm <- data$tau.i - tau
    
    ## Bootstrap iterations
    W.bern <- sapply(1:nBoot, function(x) rbinom(N, 1, prob = (sqrt(5) - 1)/(2*sqrt(5))))
    Wstar <- (((sqrt(5) + 1)/2)*W.bern + (-(sqrt(5) - 1)/2)*(1-W.bern))/N
    
    ## Apply bootstrap weights to each "observation"
    Tstar <- (t(Wstar) %*% tau.norm)
    
    ### Asypmtotic estimate is null
    se.est2 <- NULL
    
  }
  
  ### Return output
  return(list(estimate=tau, boot.dist=as.vector(Tstar), se.asymp=se.est2))
  
}