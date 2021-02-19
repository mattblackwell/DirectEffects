#' Perform telescope matching to estimate the controlled
#' direct effect of a binary treatment net the effect of a binary mediator.
#'
#' @param formula A formula object that specifies the outcome,
#' baseline covariates, treatment, intermediate covariate, and
#' mediator to be used in the matching. Each of the last four variable
#' groups should be separated by \code{|}. See below for more details.
#' @param data A dataframe containing columns referenced by
#' \code{formula}.
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
#' @param verbose logical indicating whether to display progress
#' information. Default is \code{TRUE}.
#' @param subset A vector of logicals indicating which rows of
#' \code{data} to keep.
#' @param separate_bc logical indicating whether or not bias
#' correction regressions should be run separately within levels of
#' the treatment and mediator. Defaults to \code{TRUE}. If
#' \code{TRUE}, any interactions between treatment/mediator and
#' covariates in the specification should be omitted.
#' @inheritParams stats::lm
#' 
#' @details The \code{telescope_match} function implements the two-stage
#' "telescope matching" procedure developed by Blackwell and Strezhnev 
#' (2019).
#' 
#'  The procedure first estimates a demediated outcome using a combination
#'  of matching and a regression bias-correction based on the
#' specification in \code{formula}, which should be specified as
#' \code{Y ~ X | A | Z | M}, where \code{Y} is the outcome, \code{X}
#' is a formula of baseline covariates, \code{A} is a single variable
#' name indicating the binary treatment, \code{Z} is a formula of
#' intermediate covariates, and \code{M} is a single variable name
#' indicating the mediator.
#'
#' Under the default \code{separate_bc == TRUE}, the  first stage will
#' match for \code{M} on \code{X} and \code{Z} within levels of
#' \code{A}. The bias correction regressions will be linear with a
#' specification of \code{Y ~ X + Z} within levels of \code{A} and
#' \code{M}. In the second stage, we match for \code{A} within levels
#' of \code{X} and run bias correction regressions of \code{Ytilde ~ X}
#' within levels of \code{A}, where \code{Ytilde} is the imputed value
#' of \code{Y(A, 0)} from the first stage. This second step estimates
#' the controlled direct effect of treatment when the mediator is 0.
#'
#' When \code{separate_bc} is \code{FALSE}, the bias correction
#' regressions are not broken out by the treatment/mediator and
#' \code{A} and \code{M} are simply included as separate terms in the
#' linear regression. In this setting, interactions between the
#' treatment/mediator and covariates can be added on a selective basis
#' to the \code{X} and \code{Z} specifications. 
#'  
#'  Matching is performed using the \code{Match()} routine from the
#' \code{Matching} package. By default, matching is L-to-1 nearest
#' neighbor with replacement using Mahalanobis distance. 
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
#' jobcorps_female <- subset(jobcorps, female == 1)
#' 
#' ## Telescope matching formula - First stage (X and Z)
#' tm_form <- exhealth30 ~  schobef + trainyrbef + jobeverbef  |
#' treat | emplq4 + emplq4full | work2year2q
#'  
#' 
#' ### Estimate ACDE for women holding employment at 0
#' tm_out <-  telescope_match(
#'   tm_form,
#'   data = jobcorps_female,
#'   L = 3,
#'   boot = FALSE,
#'   verbose = TRUE
#' )
#' 
#' @export
#' @importFrom Matching Match
#' @importFrom stats as.formula model.frame predict rbinom var
#'
telescope_match <- function(formula, data, caliper = NULL, L = 5,
                            boot = FALSE, nBoot = 5000, ci = 95,
                            verbose = TRUE, subset, contrasts = NULL,
                            separate_bc = TRUE, ...) {
  
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

  cl <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)

  m <- match(
    x = c("formula", "data", "subset", "na.action"),
    table = names(mf),
    nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE
  
  ## must be valid formula
  formula <- Formula::as.Formula(formula)
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:4)
  if (inherits(try(terms(formula), silent = TRUE), "try-error")) {
    stop("cannot use dot '.' in formulas")
  }

  mt_x <- terms(formula, data = data, rhs = 1)
  mt_a <- terms(formula, data = data, rhs = 2)
  mt_xz <- terms(formula, data = data, rhs = c(1,3))
  mt_m <- terms(formula, data = data, rhs = 4)
  aname <- attr(mt_a, "term.labels")
  mname <- attr(mt_m, "term.labels")
  yname <- rownames(mt_x$factors)[1]
  
  if (length(mname) > 1) stop("only a single mediator allowed.")
  if (length(aname) > 1) stop("only a single treatment allowed.")
  
  ## add to mf call
  mf$formula <- formula
  #  finally evaluate model.frame, create data matrix
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # terms object
  included <- rownames(data) %in% rownames(mf)

  all_data <- model.matrix(mt, mf, contrasts)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt_x, mf, contrasts)
  XZ <- model.matrix(mt_xz, mf, contrasts)
  A <- all_data[, aname, drop = FALSE]
  M <- all_data[, mname, drop = FALSE]


  ## pick which columns need to be dropped from X and XZ for matching.
  ## any column that is a function A or M and the intercept
  a_in_x <- which(rownames(attributes(mt_x)$factors) == aname)
  a_in_xz <- which(rownames(attributes(mt_xz)$factors) == aname)
  m_in_xz <- which(rownames(attributes(mt_xz)$factors) == mname)

  if (separate_bc & length(c(a_in_x, a_in_xz, m_in_xz))) {
    stop(
      "Treatment and/or mediator terms in covariate ",
      "specifications are not permitted with separate_bc = TRUE.",
      call. = FALSE
    )
  }
  x_drop <- which(attributes(mt_x)$factors[a_in_x,] == 1)
  if (colnames(X)[1] == "(Intercept)") x_drop <- c(1, x_drop + 1)

  xz_drop <- c(
    which(attributes(mt_xz)$factors[a_in_xz,] == 1),
    which(attributes(mt_xz)$factors[m_in_xz,] == 1)
  )
  xz_drop <- unique(xz_drop)
  if (colnames(XZ)[1] == "(Intercept)") xz_drop <- c(1, xz_drop + 1)

  ### Number of observations
  ### we don't use nrow(data) since it might be subsetted
  N <- length(Y)

  ####################
  #### Sanity checks 
  ####################

  ### Check -2 - Is the ci within (0, 100).
  if (ci <= 0 | ci >= 100){
    warning("Warning: 'ci' must be within the interval (0, 100). Defaulting to 95% confidence intervals.")
    ci <- 95
  }
  ci.alpha <- 1 - ci/100
  
  
  
  ## Check 3 - Treatment is binary 0-1
  if (!(is.numeric(A))) {
    stop("Error: 'treatment' variable not numeric.")
  } else if(!isTRUE(all.equal(unique(A)[order(unique(A))], c(0, 1)))) {
    stop("Error: 'treatment' must have only two levels: 0,1")
  }
  
  ## Check 4 - Mediator is binary 0-1
  if (!(is.numeric(M))){
    stop("Error: 'mediator' variable not numeric.")
  } else if(!isTRUE(all.equal(unique(M)[order(unique(M))], c(0, 1)))) {
    stop("Error: 'mediator' must have only two levels: 0,1")
  }
  
  
  ######################
  ### Pre-estimation set-up
  ######################

  pre.treatment <- colnames(X[, -x_drop])
  ### Post-treatment covariates are those in first stage but not in second stage
  post.treatment <- setdiff(colnames(XZ[, -xz_drop]), colnames(X[, -x_drop]))
  
  ## Combine all covariates into a single vector
  all.covariates <- c(pre.treatment, post.treatment)
  
  ## Diagnostic, print implied pre-treatment and post-treatment covariates
  if (verbose){
    cat("Telescope matching setup:\n")
    cat(paste("Outcome: ", yname, "\n", sep=""))
    cat(paste("Treatment: ", aname, "\n", sep=""))
    cat(paste("Mediator: ", mname, "\n", sep=""))
    cat(paste("Pre-treatment Covariates: ", paste(pre.treatment, collapse = ", "), "\n", sep=""))
    cat(paste("Post-treatment Covariates: ", paste(post.treatment, collapse = ", "), "\n", sep=""))
    cat(paste("Number of matches in first stage (mediator): ", L_m, "\n", sep=""))
    cat(paste("Number of matches in second stage (treatment): ", L_a, "\n", sep=""))
    cat("\n")
  }
  
  ######################
  ### First-stage: Estimate Y(a, 0) conditional on pre-/post-treatment covariates
  ######################
  
  ### First stage matching - inexact on X and Z covariates, exact on A
  tm.first <- Match(Y = Y, Tr = M, X = cbind(XZ[, -xz_drop], A), M = L_m, 
                    exact = c(rep(FALSE, ncol(XZ[, -xz_drop])), TRUE),
                    BiasAdjust = FALSE, estimand = "ATT", ties = FALSE,
                    Weight = 2)

  ### Summarize input - First Stage
  if (verbose){
    cat("First-stage matching: Mediator on pre-treatment, post-treatment\n")
    cat("Number of observations with 'mediator' = 0 matched to each observation with 'mediator' = 1: ", L_m, "\n", sep="")
    cat("Number of observations with 'mediator' = 1: ", sum(M), "\n", sep = "" )
    cat("Number of observations with 'mediator' = 0: ", sum(1 - M), "\n", sep = "")
    cat("\n")
  }
  
  ### Summarize sample sizes with mediator/treatment
  n_summary <- data.frame(
    c(0,1,0,1),
    c(0,0,1,1),
    as.vector(table(A, M))
  )
  colnames(n_summary) <- c(aname, mname, "N")
  
  ### Count of # of matches for each M=0 unit in the first stage (K_L^m)
  KLm_tab <- table(tm.first$index.control)
  ### Which M=1 units did each M=0 match to?
  first_stage_matchedTo <- tapply(tm.first$index.treated, tm.first$index.control, function(x) x)
  ## Save K_L^m in data
  KLm <- rep(0, times = N)
  KLm[as.numeric(names(KLm_tab))] <- KLm_tab
  
  ### Regression imputation in the first stage
  ## Fit a regression in the set with mediator = 0
  if (separate_bc) {
    s1.reg.10 <- lm.fit(y = Y[A == 1 & M == 0], x = XZ[A == 1 & M == 0,]) 
    s1.reg.00 <- lm.fit(y = Y[A == 0 & M == 0], x = XZ[A == 0 & M == 0,])
    pred.Y.10 <- XZ[, names(s1.reg.10$coefficients)] %*% s1.reg.10$coefficients
    pred.Y.00 <- XZ[, names(s1.reg.00$coefficients)] %*% s1.reg.00$coefficients    
    pred.Y.m0 <- A * pred.Y.10 + (1 - A) * pred.Y.00
  } else {
    s1.reg <- lm(formula(mt), data = mf)
    mf.m0 <- mf
    mf[[mname]] <- 0
    pred.Y.m0 <- predict(s1.reg, data = mf.m0)
  }
  
  ## Predicted value of regression for each X under control
  
  #### Average of matched values for each unit
  y.m0.m.imp <- Y ## Matching imputation
  y.m0.m.imp[unique(tm.first$index.treated)] <- tapply(tm.first$index.control, tm.first$index.treated, function(x) mean(Y[x]))
  #### Imputed value under M=0 using regression-imputed value of the matches
  y.m0.r.imp <- pred.Y.m0 ## Regression imputation
  y.m0.r.imp[unique(tm.first$index.treated)] <- tapply(tm.first$index.control, tm.first$index.treated, function(x) mean(pred.Y.m0[x]))
  
  ### Impute the outcome under M=0 for units with M=1
  Ytilde <- Y ## Bias-corrected matching imputation for the M=1 units
  ### Average of imputations +
  ### (Regression prediction for X_i - regression prediction for all imputed units)
  Ytilde[M == 1] <- y.m0.m.imp[M == 1] + pred.Y.m0[M == 1] - y.m0.r.imp[M == 1]
  
  ################################################
  ### Second stage - Estimate E[Y(1,0) - Y(0,0)] - conditional on pre-treatment covariates
  ################################################
  
  ### Summarize input - Second stage
  if (verbose){
    cat("Second-stage matching: Treatment on pre-treatment\n")
    cat("Number of observations matched to each unit with opposite treatment: ", L_a,  "\n", sep="")
    cat("Number of observations with 'treatment' = 1: ", sum(A), "\n", sep = "" )
    cat("Number of observations with 'treatment' = 0: ", sum(1 - A), "\n", sep = "")
    cat("\n")
  }
    
  ##### Second stage regression on blipped-down outcome
  ## Ytilde = Y_i(a, 0)
  if (separate_bc) {
    badcols.a0 <- which(apply(X[A == 0,], 2, var) < 10 * .Machine$double.eps)
    badcols.a0 <- badcols.a0[colnames(X)[badcols.a0] != "(Intercept)"]
    badcols.a1 <- which(apply(X[A == 1,], 2, var) < 10 * .Machine$double.eps)
    badcols.a1 <- badcols.a1[colnames(X)[badcols.a1] != "(Intercept)"]
    
    if (!setequal(badcols.a0, badcols.a1)) {
      warning("Some baseline covariates do not vary in one treatment arm.")
    }
    
    s2.reg.a0 <- lm.fit(y = Ytilde[A == 0], x = X[A == 0,])
    pred.Y.a0 <- X[, names(s2.reg.a0$coefficients)] %*% s2.reg.a0$coefficients
    s2.reg.a1 <- lm.fit(y = Ytilde[A == 1], x = X[A == 1,])
    pred.Y.a1 <- X[, names(s2.reg.a1$coefficients)] %*% s2.reg.a1$coefficients
  } else {
    s2.reg <- lm(formula(terms(formula, rhs = 1:2)), data = mf)
    mf.a1 <- mf.a0 <- mf
    mf.a1[[aname]] <- 1
    mf.a0[[aname]] <- 0
    pred.Y.a1 <- predict(s2.reg, newdata = mf.a1)
    pred.Y.a0 <- predict(s2.reg, newdata = mf.a0)
  }
  pred.Y.A <- A * pred.Y.a1 + (1 - A) * pred.Y.a0 ### Predict factual
  pred.Y.1.A <- A * pred.Y.a0 + (1 - A) * pred.Y.a1 ## Predict counterfactual
  
  ### Match controls to treated
  tm.second.a1 <- Match(Y = Ytilde, Tr = A, X = X[, -x_drop],
                        estimand = "ATT", caliper = caliper,
                        M = L_a, ties = FALSE,
                    Weight = 2)
  
  KLa0 <- table(tm.second.a1$index.control) ## Count of matched controls - stage 2
  ### Match treateds to control
  tm.second.a0 <- Match(Y = Ytilde, Tr = A, X = X[, -x_drop],
                        estimand = "ATC", caliper = caliper,
                        M = L_a, ties = FALSE,
                        Weight = 2)
  KLa1 <- table(tm.second.a0$index.treated) ## Count of matched treated - stage 2
  
  ## Total match counts - stage 2
  KLa <- rep(0, times = N)
  KLa[as.numeric(names(KLa0))] <- KLa0
  KLa[as.numeric(names(KLa1))] <- KLa1
  
  SLmatch <- lapply(first_stage_matchedTo, function(x) sum(KLa[x])) ## Sum of KLa of matched units
  SLm <- rep(0, times = N)
  SLm[as.integer(names(SLmatch))] <- unlist(SLmatch)
  
  ## Imputed Y_i(0,0) of matches - regression
  Yhat00.r <- pred.Y.a0
  Yhat00.r[unique(tm.second.a1$index.treated)] <- tapply(tm.second.a1$index.control, tm.second.a1$index.treated, function(x) mean(pred.Y.a0[x]))

  ## Imputed Y_i(1,0) of matches - regression
  Yhat10.r <- pred.Y.a1
  Yhat10.r[unique(tm.second.a0$index.control)] <- tapply(tm.second.a0$index.treated, tm.second.a0$index.control, function(x) mean(pred.Y.a1[x]))
  
  ## regression imputations of the CEF under 1-A_i
  ## (1/L) \sum_{j \in J^a(i)} mu_{1-A_i,0}(X_j, 1-A_i)
  pred.Y.1.A.r <- A * Yhat00.r + (1 - A) * Yhat10.r
  
  ## Linearized form for bootstrapping
  sm.part <- (1 - M) * (1 + KLa/L_a + KLm/L_m + SLm/(L_a * L_m)) * Y
  bm.part <- ((1 - M) * (KLm/L_a + SLm/(L_a * L_m)) -
                M * (1 + KLa/L_a)) * pred.Y.m0
  ba.part <- pred.Y.1.A + (KLa/L_a) * pred.Y.A
  
  tau.i <- (2 * A - 1) * (sm.part - bm.part - ba.part)
  
  ## Point estimate
  tau <- mean(tau.i)

  ########################################
  ## Variance estimation
  ########################################
  
  if (boot == FALSE){
    
    ###############
    ### Asymptotic variance
    ###################
    
    ### Calculate weights
    ww <- (1 - M) * (1 + KLa/L_a + KLm/L_m + SLm/(L_a * L_m))
    
    ### Number of units with M=0
    N0 <- sum(1 - M)
    ### Parameters in the first-stage regression
    if (separate_bc) {
      P1 <- length(coef(s1.reg.10)) + length(coef(s1.reg.00))
      P2 <- length(coef(s2.reg.a1)) + length(coef(s2.reg.a0)) 
    } else {
      P1 <- length(coef(s1.reg))
      P2 <- length(coef(s2.reg))

    }
    
    ### Variance component 1
    em.var <- (N0/(N0 - P1)) * (Y - pred.Y.m0)^2
    
    ### Variance component 2
    ea.var <- (N/(N-P2)) * (pred.Y.m0 - pred.Y.A)^2
    
    ### Variance component 3
    tau.var <- mean((pred.Y.a1 - pred.Y.a0 - tau)^2)
    
    ### combine all three components
    se.est2 <- sqrt(tau.var/N +
                      (mean((1 - M) * ww^2 * em.var) +
                         mean((1 + KLa / L_a)^2 * ea.var))/N)
    
    ### No bootstrap
    Tstar <- NULL
    
    ### CI bounds
    ci.low <- tau - abs(qnorm(ci.alpha/2)) * se.est2
    ci.high <- tau + abs(qnorm(ci.alpha/2)) * se.est2
    
  }else{
    
    ################
    ### Wild Bootstrap
    ################
    
    ## De-mean
    tau.norm <- tau.i - tau
    
    ## Bootstrap iterations
    W.bern <- sapply(1:nBoot, function(x) rbinom(N, 1, prob = (sqrt(5) - 1)/(2*sqrt(5))))
    Wstar <- (((sqrt(5) + 1)/2) * W.bern + (-(sqrt(5) - 1)/2) * (1-W.bern))/N
    
    ## Apply bootstrap weights to each "observation"
    Tstar <- (t(Wstar) %*% tau.norm)
    
    ## Add in tau
    Tstar <- Tstar + tau
    
    ### Asypmtotic estimate is null
    se.est2 <- NULL
    
    ### Get quantiles for the CI
    ci.low <- quantile(Tstar, ci.alpha/2)
    ci.high <- quantile(Tstar, 1 - ci.alpha/2)
  }
  
  ### Return output
  output <- list(formula = formula, outcome = yname, treatment = aname,
                 mediator = mname, included = included,
                 N = N, L_m = L_m, L_a = L_a, N_summary = n_summary, 
                 estimate=tau, std.err = se.est2, boot.dist=as.vector(Tstar), KLm = KLm,
                 KLa = KLa, outcome.vec = Y, treatment.vec = A, mediator.vec = M,
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
#' @param object an object of class \code{tmatch} -- results from a
#' call to \code{telescope_match}
#' @param ... additional arguments affecting the summary produced.
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
summary.tmatch <- function(object, ...){
  
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
#' @importFrom stats as.formula model.frame


balance.tmatch <- function(object, vars, data){
  
  ##################
  ### Sanity checks
  ##################
  
  ### Is the class a 'tmatch'
  if (class(object) != "tmatch"){
    stop("Error: 'object' not of class 'tmatch'")
  }
  
  ### Does N of data match number of obs 
  if (nrow(data) != length(object$included)){
    stop("Error: number of rows in data not consistent with object")
  }

  data_tr <- data[[object$treatment]][object$included]
  if (!all.equal(as.numeric(data_tr),
                 as.numeric(object$treatment.vec))) {
    stop("Error: treatment vector in data and object don't match")
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

  ## subset to the selected rows used in matching
  data <- data[object$included,]
  
  ## Do model.frame first to parse any functions
  covariate.frame = tryCatch({model.frame(vars, data)},error = function(e) { stop("Could not extract all variables in 'formula' from data")})
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
#' @importFrom graphics hist

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
    plot_title = paste("Matching weights for first stage (mediator)\namong M = 0")
    hist(object$KLm[object$mediator.vec == 0]/object$L_m, main=plot_title, xlab="Number of times unit is matched")
    abline(v = 1, col="red", lty=2, lwd=2)
  }else if (stage == "treatment"){
    plot_title = paste("Matching weights for second stage (treatment)")
    hist(object$KLa/object$L_a, main=plot_title, xlab="Number of times unit is matched")
    abline(v = 1, col="red", lty=2, lwd=2)
  }
  
  return()
  
}
