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
sequential_g <- function(formula, data, subset, weights, na.action,
                         offset, contrasts = NULL, verbose = TRUE, ...) {

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
  mt_m <- terms(formula, data = data, rhs = 3)
  mnames <- attr(mt_m, "term.labels")

  ## add to mf call
  mf$formula <- formula

  #  finally evaluate model.frame, create data matrix
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # terms object

  XZM <- model.matrix(mt, mf, contrasts)
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt_x, mf, contrasts)
  M <- XZM[, mnames, drop = FALSE]

  ### weights
  w <- as.vector(model.weights(mf))

  ### offset parameter
  offset <- as.vector(model.offset(mf))

  ## regress de-mediated outcome on Xs by OLS or WLS
  if (is.null(w)) {
    out.first <- lm.fit(XZM, Y, offset = offset, ...)
    mcoefs <- out.first$coefficients[mnames]
    Y.blip <- Y - M %*% mcoefs
    out <- lm.fit(X, Y.blip, offset = offset, ...)
  } else {
    out.first <- lm.wfit(XZM, Y, w, offset = offset, ...)
    mcoefs <- out.first$coefficients[mnames]
    Y.blip <- Y - M %*% mcoefs
    out <- lm.wfit(X, Y.blip, w, offset = offset, ...)
  }
  out.first$aliased <- is.na(out.first$coefficient)
  out.first$terms <- mt
  out.first$contrasts <- attr(XZM, "contrasts")

  ## Combine all components
  out$aliased <- is.na(out$coefficient)

  out$terms <- list(M = mt_m, X = mt_x)
  out$formula <- formula
  out$call <- cl
  out$na.action <- attr(mf, "na.action")
  out$xlevels <- stats::.getXlevels(mt, mf)
  out$contrasts <- list(X = attr(X, "contrasts"), M = attr(M, "contrasts"))

  out$model <- mf
  out$X <- X
  out$M <- M
  out$Ytilde <- Y.blip
  out.first$XZM <- XZM

  out$first_mod <- out.first
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
  Fhat <- crossprod(X2, X1) / n
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

  ### "bread" of the consistent variance estimator
  p <- 1L:direct_mod$rank
  bread <- chol2inv(direct_mod$qr$qr[p, p, drop = FALSE])

  ### degrees of freedom correction
  dfc <- (n / (n - max(p)))

  vv <- dfc * (bread %*% meat %*% bread)
  return(vv)
}
