
#' Bootstrap standard errors
#'
#' Called internally in sequential_g
#'
#' @param seqg Output from sequential_g
#' @param boots The number of bootstrap replicates. Defaults to 1000.
#' @param progress Show a textbar for progress. Defaults to TRUE.
#' @param data Dataframe to be used if none provided in seqg.
#'
#'
#' @keywords internal


boots_g <- function(seqg, boots = 1000, progress = TRUE, data = NULL) {
  acde.boots <- rep(list(NA), times = boots) # holder for boot strap estimates
  if (progress) prog.bar <- utils::txtProgressBar(min = 0, max = boots, style = 3) # start progress bar
  ## must be valid formula
  formula <- Formula::as.Formula(seqg$formula)
  
  for (b in 1:boots) {
    if (progress) utils::setTxtProgressBar(prog.bar, b) # update progress bar
    draw <- sample(1:nrow(seqg$model), replace = TRUE) # vector for sample with replacement
    # bootstrap sampling
    if(!is.null(seqg$model)){ # check if data is contained in seqg$model
    data.draw.first <- seqg$model[draw, ] # data for re-estimation of first model
    data.draw.direct <- seqg$first_mod$model[draw, ] # data for re-estimation of direct effects model
    
    
    # combine data
    shared <- colnames(data.draw.direct)[colnames(data.draw.direct) %in% colnames(data.draw.first)]
    just.direct <- colnames(data.draw.direct)[!(colnames(data.draw.direct) %in% colnames(data.draw.first))]
    just.first <- colnames(data.draw.first)[!(colnames(data.draw.first) %in% colnames(data.draw.direct))]
    data.draw <- cbind(
      data.draw.direct[, shared],
      data.draw.direct[, just.direct],
      data.draw.first[, just.first]
    )
    colnames(data.draw) <- c(shared, just.direct, just.first)
    } else if(is.data.frame(data)){ # if data is not in seqg$model, use data argument
      data.draw <- data[draw,]
    }
  
    ### re-estimate first model on bootstrapped data
    first_mod <- lm(seqg$first_mod$call[[2]], data = data.draw) # estimate first model

    # construct mediator function variables
    mf <- match.call(expand.dots = FALSE)

    m <- match(
      x = c("formula", "data", "subset", "na.action", "weights", "offset"),
      table = names(mf),
      nomatch = 0L
    )
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf$drop.unused.levels <- TRUE
  
    
   
    ## update formula Y ~ A + M - 1
    f1 <- formula(formula, rhs = 1) # Y ~ A + X
    f2 <- formula(formula, lhs = 0, rhs = 2) # ~ M
    f2 <- update(f2, ~. - 1) # ~ M - 1, don't model intercept
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
    mf$data = data.draw
    
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
      boot.direct <- lm.fit(X, Y, offset = offset)
    } else {
      boot.direct <- lm.wfit(X, Y, w, offset = offset)
    }
    
    
    #####
    
   # boot.direct <- sequential_g(seqg$formula, first_mod = boot.first, data = data.draw) # estimate direct effects model
    acde.boots[[b]] <- coef(boot.direct) # store direct effect coefficients
  }
  
  if (progress) close(prog.bar)
  
 # combine list into matrix
  acde.boots <- do.call(rbind, acde.boots)

  # construct output
  out <- list()
  out$acde <- acde.boots # matrix of bootstrap coefficients
  out$acde.sd <- sapply(as.data.frame(acde.boots), sd) # bootstrapped standard errors
  
  
  return(out)
}