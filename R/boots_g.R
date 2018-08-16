#' Bootstrap Standard Errors for Sequential G model object
#'
#' Estimate bootstrap standard errors
#'
#' @param seqg Output from sequential_g
#' @param boots The number of bootstrap replicates. Defaults to 1000.
#' @param showbar Show a textbar for progress?
#'
#'
#' @export
#' @return
#' \item{acde.sd}{a vector of bootstrapped standard errors for each named coefficient for the direct effects model.}
#' \item{ate.sd}{a vector of bootstrapped stanadard errors for each named coefficient for the
#' average treatment effects effects model.}
#' \item{pt.sd}{a vector of bootstrapped stanadard errors for each named coefficient for the
#' first stage effects model.}
#'
#'
#' @examples
#' data(civilwar)
#'
#' # First stage
#' fit_first <- lm(onset ~ ethfrac + warl + gdpenl + lpop +
#'                   ncontig + Oil + nwstate + polity2l + relfrac + instab,
#'                 data = civilwar)
#'
#'
#' # main formula: Y ~ A + X | M
#' form_main <- onset ~ ethfrac + lmtnest + ncontig + Oil | instab
#'
#' # estimate CDE
#' direct <- sequential_g(form_main, fit_first, data = civilwar)
#'
#' # bootstrap
#' out_boot <- boots_g(direct)
#'
#' # print bootstrapped standard errors
#' print(out_boot)
#'


boots_g <- function(seqg, boots = 1000, showbar = TRUE) {
  
  acde.boots <- ate.boots <- pt.boots <- rep(list(NA), times = boots) # holder for boot strap estimates
  
  if (showbar) prog.bar <- utils::txtProgressBar(min = 0, max = boots, style = 3) # start progress bar
  
  for (b in 1:boots) {
    if (showbar) utils::setTxtProgressBar(prog.bar, b) # update progress bar

    # bootstrap sampling
    draw <- sample(1:nrow(seqg$model), replace = TRUE) # vector for sample with replacement
    data.draw.first <- seqg$first_mod$model[draw, ] # data for re-estimation of first model
    data.draw.direct <- seqg$model[draw, ] # data for re-estimation of direct effects model

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

    # estimate models
    ate.mod <- lm(formula(seqg$formula, lhs = 1, rhs = 1), data = data.draw) # estimate ATE
    ate.boots[[b]] <- coef(ate.mod) # store ate coefficients
    
    boot.first <- lm(seqg$first_mod$call[[2]], data = data.draw) # estimate first model
    pt.boots[[b]] <- coef(boot.first) # store first model coefficients
    
    boot.direct <- sequential_g(seqg$formula, first_mod = boot.first, data = data.draw) # estimate direct effects model
    acde.boots[[b]] <- coef(boot.direct) # store direct effect coefficients
  }
  
  if (showbar)  close(prog.bar)

  # combine lists into matrices
  ate.boots  <- do.call(rbind, ate.boots)
  acde.boots <- do.call(rbind, acde.boots)
  pt.boots   <- do.call(rbind, pt.boots)

  # construct output
  out <- list()
  out$acde.sd  <- sapply(as.data.frame(acde.boots), sd)
  out$ate.sd   <- sapply(as.data.frame(ate.boots), sd)
  out$pt.sd    <- sapply(as.data.frame(pt.boots), sd)


  return(out)
}
