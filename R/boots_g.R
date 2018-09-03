
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
