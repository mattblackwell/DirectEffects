#' Perform linear sequential g-estimation to estimate the controlled
#' direct effect of a treatment net the effect of a mediator.
#'
#' @inheritParams stats::lm
#' @param formula formula specification of the first-stage,
#'   second-stage, and blip-down models. The right-hand side of the
#'   formula should have three components separated by the \code{|},
#'   with the first component specifying the first-stage model with
#'   treatment and any baseline covariates, the second component
#'   specifying the intermediate covariates for the first-stage, and
#'   the third component specifying the blip-down model. See Details
#'   below for more information.
#' @param data A dataframe to apply \code{formula} on.
#' @param subset A vector of logicals indicating which rows of \code{data} to keep.
#' @param verbose logical indicating whether to suppress progress bar. Default is FALSE.
#' @details The \code{sequential_g} function implements the linear
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
#' @return Returns an object of \code{class} A \code{"seqg"}. Similar
#' to the output of a call to \code{lm}. Contains the following
#' components:
#' \itemize{
#'   \item coefficients: a vector of named coefficients for the direct
#' effects model.
#'   \item residuals: the residuals, that is the blipped-down outcome
#' minus the fitted values.
#'   \item rank: the numeric rank of the fitted linear direct effects
#' model.
#'   \item fitted.values: the fitted mean values of the direct effects
#' model.
#'   \item weights: (only for weighted fits) the specified weights.
#'   \item df.residual: the residual degrees of freedom for the direct
#' effects model.
#'   \item aliased: logical vector indicating if any of the terms were
#'   dropped or aliased due to perfect collinearity.
#'   \item terms: the list of \code{terms} object used. One for the
#'   baseline covariates and treatment (\code{X}) and one for the
#'   variables in the blip-down model (\code{M}).
#'   \item formula: the \code{formula} object used, possibly modified
#' to drop a constant in the blip-down model.
#'   \item call: the matched call.
#'   \item na.action: (where relevant) information returned by
#'   \code{model.frame} of the special handling of \code{NA}s.
#'   \item xlevels: the levels of the factor variables.
#'   \item contrasts:  the contrasts used for the factor variables.
#'   \item first_mod: the output from the first-stage regression model.
#'   \item model: full model frame, including all variables.
#'   \item Ytilde: the blipped-down response vector.
#'   \item X: the model matrix for the second stage.
#'   \item M: the model matrix for demediation/blip-down function.
#' }
#' In addition, non-null fits will have components \code{assign},
#' \code{effects}, and \code{qr} from the output of \code{lm.fit} or
#' \code{lm.wfit}, whichever is used.
#' @references Vansteelandt, S. (2009). Estimating Direct Effects in
#' Cohort and Case-Control Studies. Epidemiology, 20(6), 851-860.
#'
#' Acharya, Avidit, Blackwell, Matthew, and Sen, Maya. (2016)
#' "Explaining Causal Effects Without Bias: Detecting and Assessing
#' Direct Effects." American Political Science Review 110:3 pp.
#'   512-529
#' @examples
#' data(ploughs)
#'
#' form_main <- women_politics ~ plow +
#'   agricultural_suitability + tropical_climate + large_animals +
#'   political_hierarchies + economic_complexity +
#'   rugged | years_civil_conflict +
#'   years_interstate_conflict  + oil_pc +
#'   european_descent + communist_dummy + polity2_2000 +
#'   serv_va_gdp2000 | centered_ln_inc + centered_ln_incsq
#'
#' direct <- sequential_g(form_main, ploughs)
#'
#' summary(direct)
#' @export
#' @importFrom stats coef lm.fit lm.wfit model.matrix model.offset
#'   model.response model.weights pt residuals terms update
#'
telescope_match <- function(formula, treatment, data, subset, bc = TRUE, ci = 0.95,
                            L = 5, n_boot = 5000, verbose = TRUE, ...) {

  # store model calls
  cl <- match.call(expand.dots = TRUE)
  mf <- match.call(expand.dots = FALSE)

  m <- match(
    x = c("formula", "data", "subset", "na.action", "weights", "offset"),
    table = names(mf),
    nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE

  ## must be valid formula
  formula <- Formula::as.Formula(formula)
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:3)
  if (inherits(try(terms(formula), silent = TRUE), "try-error")) {
    stop("cannot use dot '.' in formulas")
  }

  mt_x <- terms(formula, data = data, rhs = 1)
  mt_z <- terms(formula, data = data, rhs = 2)
  mt_m <- terms(formula, data = data, rhs = 3)
  mt_xzm <- terms(formula, data = data, rhs = 1:3)
  attr(mt_m, "intercept") <- 0

  xnames <- attr(mt_x, "term.labels")
  mnames <- attr(mt_m, "term.labels")
  znames <- attr(mt_z, "term.labels")

  stopifnot(treatment %in% xnames)
  stopifnot(length(mnames) == 1)

  ## add to mf call
  mf$formula <- formula

  #  finally evaluate model.frame, create data matrix
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # terms object

  XAZ <- model.matrix(mt_xaz, mf, contrasts)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt_x, mf, contrasts)
  M <- model.matrix(mt_m, mf, contrasts)
  A <- mf[,treatment]

  trpos <- which(xnames == treatment)
  s1.vars <- c(xnames, znames)
  ex.list <- rep(FALSE, length(s1.vars))
  ex.list[trpos] <- TRUE


  tm.first <- Match(Y = Y, Tr = M, X = XAZ[,s1.vars], exact = ex.list,
                    M = L, BiasAdjust = FALSE, estimand = "ATT")

    ### Summarize input - First Stage
  if (verbose) {
    cat("First-stage matching: Mediator on pre-treatment, post-treatment\n")
    cat("Number of observations with 'mediator' = 1: ", sum(M), "\n", sep = "" )
    cat("Number of observations with 'mediator' = 0: ", sum(1 - M), "\n", sep = "")
    cat("\n")
  }

  ### Count of # of matches for each M=0 unit in the first stage
  KLm <- table(tm.first$index.control)
  ### Which M=1 units did each M=0 match to?
  first_stage_matchedTo <- tapply(tm.first$index.treated, tm.first$index.control, function(x) x)
  ## Save K_L^m in data
  data$KLm <- 0
  data[as.numeric(names(KLm)), "KLm"] <- KLm


  ## regress de-mediated outcome on Xs by OLS or WLS
  out.first <- lm.fit(XAZ[M == 0,], Y[M == 0])
  reg_y_m0 <- out.first$coefficients %*% XAZ

  y_m0_m_imp <- Y
  y_m0_m_imp[unique(tm.first$index.treated)] <- tapply(tm.first$index.control, tm.first$index.treated, function(x) mean(Y[x]))

  y_m0_r_imp <- reg_y_m0
  y_m0_r_imp[unique(tm.first$index.treated)] <- tapply(tm.first$index.control, tm.first$index.treated, function(x) mean(reg_y_m0[x]))

  Y_blip <- Y
  Y_blip[M == 1] <- y_m0_m_imp[M == 1] + (reg_y_m0[M == 1] - y_m0_r_imp[M == 1])


  trpos <- which(colnames(X) == treatment)
  s2.reg.a0 <- lm.fit(X[A == 0], Y_blip[A == 0])
  s2.reg.a1 <- lm.fit(X[A == 1], Y_blip[A == 1])
  reg_y_a0 <- s2.reg.a0$coefficients %*% X
  reg_y_a1 <- s2.reg.a1$coefficients %*% X

  reg_y_a <- A * reg_y_a1 + (1 - A) * reg_y_a0
  reg_y_1_a <- A * reg_y_a0 + (1 - A) * reg_y_a1


    ### Match controls to treated
  tm.second.a1 <- Match(Y = Y_blip, Tr = A, X = X, estimand = "ATT", M = L)
  KLa0 <- table(tm.second.a1$index.control) ## Count of matched controls - stage 2
  ### Match treateds to control
  tm.second.a0 <- Match(Y = Y_blip, Tr = A, X = X, estimand = "ATC", M = L)
  KLa1 <- table(tm.second.a0$index.treated) ## Count of matched treated - stage 2

  ## Total match counts - stage 2
  data$KLa <- 0
  data[names(KLa0), "KLa"] <- KLa0
  data[names(KLa1), "KLa"] <- KLa1

  SLmatch <- lapply(first_stage_matchedTo, function(x) sum(data$KLa[x])) ## Sum of KLa of matched units
  data$SLm <- 0
  data$SLm[as.integer(names(SLmatch))] <- unlist(SLmatch)


  return(out)
}
