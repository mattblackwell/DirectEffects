
#' Bootstrap standard errors
#'
#' Called internally in sequential_g
#'
#' @param seqg Output from sequential_g
#' @param boots The number of bootstrap replicates. Defaults to 1000.
#' @param progress Show a textbar for progress. Defaults to TRUE.
#' @param data Dataframe to be used if none provided in seqg.
#' @param rho Vector to be used if called from cdesens().
#'
#'
#' @keywords internal


boots_g <- function(seqg, boots = 1000, progress = TRUE, data = NULL, rho = NULL) {
  # if called from cdesens, create acde.sens instead of acde.boots object
  if(!is.null(rho)){
    acde.sens <- matrix(NA, nrow = boots, ncol = length(rho)) # bootstrap samples as rows
  } else{
    acde.boots <- rep(list(NA), times = boots) # holder for boot strap estimates
  }
  if (progress) prog.bar <- utils::txtProgressBar(min = 0, max = boots, style = 3) # start progress bar
  ## must be valid formula
  formula <- Formula::as.Formula(seqg$formula)

  for (b in 1:boots) {
    if (progress) utils::setTxtProgressBar(prog.bar, b) # update progress bar
    draw <- sample(1:nrow(data), replace = TRUE) # vector for sample with replacement
    # bootstrap sampling

    data.draw <- data[draw, ]

    ### re-estimate first model on bootstrapped data
    first_mod <- update(seqg$first_mod, data = data.draw) # estimate first model
    fcoefs <- coef(first_mod)

    bnames <- attr(seqg$terms$blip, "term.labels") # M names
    bvars <- match(bnames, names(fcoefs), 0L) # which terms in first stage are Ms



    X <- model.matrix(seqg$terms$ax, data = data.draw, seqg$contrasts)
    M <- model.matrix(seqg$terms$blip, data = data.draw, seqg$contrasts)
    if(!is.null(rho)){

      # residuals
      # epsilon.tilde.i.m: residuals of mediation function
      res.m <- residuals(lm.fit(x = X, y = M))
      
      # epsilon.tilde.i.y: all variables in first model except medvar
      form.first.y.X <- update(seqg$first_mod$terms, paste0(". ~ . -", medvar))
      res.y <- residuals(lm(form.first.y.X, data = data.draw))
      
      rho.tilde <- cor(res.y, res.m)
      
      # for each value of the vector rho, change mediator value
      rho.factor <- rho * sqrt((1 - rho.tilde^2) / (1 - rho^2))
      m.fixed <- coef(first_mod)[medvar] - sd(res.y) * rho.factor / sd(res.m)
      
      # calculate acde at each rho (indexed by r)
      for (r in 1:length(rho)) {
        
        # create ytilde by blipping down with rho
        Ytilde <- model.response(data.draw, "numeric") - m.fixed[r] * (data.draw[[medvar]])
        mf.r <- cbind(Ytilde, data.draw)
        
        # run Ytilde ~ A + X
        sens.direct.r <- lm(form.Ytilde, data = mf.r)
        
        # save final estimate
        acde.sens[b, r] <- coef(sens.direct.r)[trvar]
      } # close rho loop
    } else {
    
    
    
    rawY <- model.response(data.draw, "numeric")
    gamma <- M %*% fcoefs[bvars] # fitted values from de-mediation function
    Y <- rawY - gamma

    ### weights
    w <- as.vector(model.weights(data.draw))

    ### offset parameter
    offset <- as.vector(model.offset(data.draw))

    ## regress de-mediated outcome on Xs by OLS or WLS
    if (is.null(w)) {
      boot.direct <- lm.fit(X, Y, offset = offset)
    } else {
      boot.direct <- lm.wfit(X, Y, w, offset = offset)
    }


    #####

    acde.boots[[b]] <- coef(boot.direct) # store direct effect coefficients
    }
  }

  if (progress) close(prog.bar)

  # if called from cdesens, only output acde.sens
  if(!is.null(rho)){
    out <- acde.sens
    return(out)
  } else {
    # combine list into matrix
    
  acde.boots <- do.call(rbind, acde.boots)

  # construct output
  out <- list()
  out$acde <- acde.boots # matrix of bootstrap coefficients
  out$acde.sd <- sapply(as.data.frame(acde.boots), sd) # bootstrapped standard errors


  return(out)
  } 
  }

