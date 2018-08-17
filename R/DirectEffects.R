#' Perform linear sequential g-estimation to estimate the controlled
#' direct effect of a treatment net the effect of a mediator.
#'
#' @inheritParams stats::lm
#' @param formula formula specification of the direct effect and blip-down models.
#' Should be of the form \code{y ~ tr + x1 + x2 | med} where \code{tr} is the
#' name of the treatment variable, \code{med} is the name of the mediator
#' and \code{x1} and \code{x2} are baseline covariates. Before the \code{|} bar
#' represents the direct effects model and after the bar represents
#' the blip-down model, the latter of which will be used to created
#' the blipped down outcome.
#' @param first_mod an \code{lm} output containing the first-stage
#' regression model. Must contain a coefficient for all variables in
#' the blip-down model in the \code{formula} argument.
#' @param data A dataframe to apply \code{formula} on. Must be identical to the 
#' \code{data} used for \code{first_mod}.
#' @param subset A vector of logicals indicating which rows of \code{data} to keep.
#' The rows of \code{data} that are not in the model matrix of \code{first_mod}
#' will always be dropped, in addition to any user-specified further subsets.
#' @param model logical indicating whether the resulting model frame
#' should be returned.
#' @param y logical indicating whether the blipped-down outcome
#' vector should be returned.
#' @param x logical indicating whether the model matrix of the direct
#' effects model should be returned.
#' @details The \code{sequential_g} function implements the linear
#' sequential g-estimator developed by Vansteelandt (2009) with the
#' consistent variance estimator developed by Acharya, Blackwell, and
#' Sen (2016).
#'
#' The function takes in a first-stage linear model, \code{first_mod},
#' of the \code{lm} class which is assumed to estimate the effect of
#' the mediator on the outcome (conditional on the treatment,
#' intermediate confounders, and baseline confounders). The formula
#' specifies both the second stage, direct effect model of the
#' treatment and baseline confounders, but also contains the
#' specification of the 'blip-down' or 'demediation' function that is
#' used to remove the average effect of the mediator (possibly
#' interacted) from the outcome to create the blipped-down outcome.
#' This blipped-down outcome is the passed to a standard linear model
#' with the covariates as specified for the direct effects model.
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
#'   \item terms: the \code{terms} object used.
#'   \item formula: the \code{formula} object used, possibly modified
#' to drop a constant in the blip-down model.
#'   \item call: the matched call.
#'   \item contrasts:  the contrasts used for the factor variables.
#'   \item levels: the levels of the factor variables.
#'   \item first_mod: the first_mod object used.
#'   \item vcov: covariance matrix of the direct-effects coefficients.
#'   \item model: full model frame (if \code{model = TRUE}).
#'   \item y: the blipped-down response vector (if \code{y = TRUE}).
#'   \item x: the direct effects model matrix (if \code{x = TRUE}).
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
#' ploughs$centered_ln_inc <- ploughs$ln_income - mean(ploughs$ln_income, na.rm = TRUE)
#' ploughs$centered_ln_incsq <- ploughs$centered_ln_inc^2
#'
#' fit_first <- lm(women_politics ~ plow + centered_ln_inc +
#'                 centered_ln_incsq + agricultural_suitability +
#'                 tropical_climate +  large_animals +
#'                 political_hierarchies + economic_complexity +
#'                 rugged + years_civil_conflict +
#'                 years_interstate_conflict  + oil_pc +
#'                 european_descent + communist_dummy + polity2_2000 +
#'                 serv_va_gdp2000, data = ploughs)
#'
#' form_main <- women_politics ~ plow + agricultural_suitability +
#'   tropical_climate +  large_animals + political_hierarchies +
#'   economic_complexity + rugged | centered_ln_inc +
#'   centered_ln_incsq
#'
#' direct <- sequential_g(formula = form_main,
#'                       first_mod = fit_first,
#'                       data = ploughs)
#'
#' summary(direct)
#' @export
#' @importFrom stats coef lm.fit lm.wfit model.matrix model.offset
#'   model.response model.weights pt residuals terms update
sequential_g <- function(formula, first_mod, data, subset, weights, na.action, model = TRUE, y = TRUE, x = FALSE, offset, contrasts = NULL, ...) {

  # store model calls
  cl <- match.call(expand.dots =  TRUE)
  mf <- match.call(expand.dots = FALSE)

  m <- match(x = c("formula", "data", "subset", "na.action", "weights", "offset"),
             table = names(mf), 
             nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE

  ## check if the two models share the same dataset
  if (!identical(first_mod$call$data, cl$data)) {
    stop("data must be the same for both models")
  }
  
  ## must be valid formula
  formula <- Formula::as.Formula(formula)
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  if (inherits(try(terms(formula), silent = TRUE), "try-error")) {
    stop("cannot use dot '.' in formulas")
  }

  ## update formula Y ~ A + M - 1
  f1 <- formula(formula, rhs = 1) # Y ~ A + X
  f2 <- formula(formula, lhs = 0, rhs = 2) # ~ M
  f2 <- update(f2, ~ . - 1) # ~ M - 1, don't model intercept
  formula <- Formula::as.Formula(f1, f2) # Y ~ A + X | M - 1


  ## link mediators across first and second models
  ## ensure that the blip formula doesn't include a constant
  bt <- terms(formula, data = data, rhs = 2) # extract terms from Y ~ M - 1
  bnames <- attr(bt, "term.labels") # their names
  fcoefs <- coef(first_mod) # beta_hat from first model (Y ~ A + M + X)
  ## are all of the Ms found in first stage?
  if (!all(bnames %in% names(fcoefs))) {
    stop("blip.form contains terms not in first_mod")
  }
  bvars <- match(bnames, names(fcoefs), 0L) # which terms in first stage are Ms

  ## add to mf call
  mf$formula <- formula
  
  #  finally evaluate model.frame, create data matrix
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # terms object
  
  ## additional subsetting -- all rows in second stage need to have been 
  ## present in first stage model matrix
  mf <- mf[which(rownames(mf) %in% rownames(model.matrix(first_mod))), , drop = FALSE]

  ### de-mediated Y
  rawY <- model.response(mf, "numeric")
  M <- model.matrix(bt, mf, contrasts)
  gamma <- M %*% fcoefs[bvars] # fitted values from de-mediation function
  Y <- rawY - gamma

  ### weights
  w <- as.vector(model.weights(mf))

  ### offset parameter
  offset <- as.vector(model.offset(mf))

  ### X (including treatment)
  mtX <- terms(formula, data = data, rhs = 1) # Y ~ A + X
  X <- model.matrix(mtX, data = mf, contrasts.arg = contrasts)


  ## regress de-mediated outcome on Xs by OLS or WLS
  if (is.null(w)) {
    out <- lm.fit(X, Y, offset = offset, ...)
  } else {
    out <- lm.wfit(X, Y, w, offset = offset, ...)
  }


  ## Combine all components
  out$terms <- list(direct = terms(formula), blip = bt, seqg = mt)
  out$formula <- formula
  out$call <- cl
  out$na.action <- attr(mf, "na.action")
  out$xlevels <- stats::.getXlevels(mt, mf)
  out$contrasts <- attr(X, "contrasts")
  out$first_mod <- first_mod

  ## Consistent Variance estimator
  out.vcov <- seq.g.vcov(first_mod, out, X1 = model.matrix(first_mod), X2 = X, bnames)
  dimnames(out.vcov) <- list(names(out$coefficients), names(out$coefficients))

  out$vcov <- out.vcov
  if (model)
    out$model <- mf
  if (x)
    out$x <- X
  if (y)
    out$y <- Y

  ## Declare class
  class(out) <- "seqg"

  return(out)
}

#' Consistent Variance estimator
#' 
#' Called internally in sequential_g
#' 
#' @param first_mod first stage model
#' @param direct_mod DirectEffects model with terms
#' @param X1 model matrix of first stage model
#' @param X2 model matrix of DirectEffects model
#' @param med.vars vector of term.label attributes for mediators
#' 
#' @keywords internal
#'  
seq.g.vcov <- function(first_mod, direct_mod, X1, X2, med.vars) {
  n <- NROW(X2)
  Fhat <- crossprod(X2, X1)/n
  Fhat[, !(colnames(X1) %in% med.vars)] <- 0

  ### invert X'X
  p1 <- 1L:first_mod$rank
  R1 <- chol2inv(first_mod$qr$qr[p1, p1, drop = FALSE])

  ## Xepsilon from second stage
  w2 <- direct_mod$weights
  res2 <- residuals(direct_mod)
  efun2 <- if (is.null(w2)) X2 * res2 else w2 * X2 * res2

  ### from first stage
  w1 <- first_mod$weights
  res1 <- residuals(first_mod)
  efun1 <- if (is.null(w1)) X1 * res1 else w1 * X1 * res1

  ### combine previous components to get covariance
  ghat <- t(efun2) + Fhat %*% R1 %*% t(efun1)
  meat <- crossprod(t(ghat))

  ### "bread" of the consisitent variance estimator
  p <- 1L:direct_mod$rank
  bread <- chol2inv(direct_mod$qr$qr[p, p, drop = FALSE])

  ### degrees of freedom correction
  dfc <- (n/(n-max(p)))

  vv <- dfc * (bread %*% meat %*% bread)
  return(vv)
}
