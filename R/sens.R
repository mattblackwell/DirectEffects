seq.g.var <- function(mod.first, mod.direct, med.vars) {
  require(sandwich)
  mat.x <- model.matrix(mod.direct)
  mat.first <- model.matrix(mod.first)
  n <- nrow(mat.x)
  Fhat <- crossprod(mat.x, mat.first)/n
  Fhat[, !(colnames(mat.first) %in% med.vars)] <- 0
  Mhat.inv <- solve(crossprod(mat.first)/n)
  ghat <- t(estfun(mod.direct)) + Fhat %*% Mhat.inv %*% t(estfun(mod.first))
  meat <- crossprod(t(ghat))/n
  bread <- (t(mat.x)%*%mat.x)/n
  vv <- (n/(n-ncol(mat.x)))*(solve(bread) %*% meat %*% solve(bread))/n
  return(vv)
}



#' Sensitivity of ACDE
#'
#' Estimate how ACDE varies by various levels of unmeasured confounding
#'
#' @param seqg Output from sequential_g
#' @param trvar character vector for treatment (abstract away later)
#' @param medvar character vector for mediator (abstract away later, allowing for more than 1)
#' @param rho A numerical vector of correlations between errors to test for. The 
#'  original model asusmes \env{rho = 0}
#' 
#' @examples
#' data(fear)  
#' # First stage
#' fit_first <- lm(onset ~ ethfrac + warl + gdpenl + lpop + 
#'                   ncontig + Oil + nwstate + polity2l + relfrac + instab,
#'                 data = fear) # lm(Y ~ A + X + M)
#' 
#' rows_use <- rownames(fear) %in% rownames(model.matrix(fit_first)) # listwise deletion
#' 
#' # main formula
#' form_main <- onset ~ ethfrac + lmtnest + ncontig + Oil | instab # Y ~ A + X | M
#' 
#' # estimate CDE
#' direct <- sequential_g(form_main, fit_first, data = fear, subset = rows_use)
#' 
#' # sensitivity 
#' out_sens <- cdesens(direct, trvar = "ethfrac", medvar = "instab")
#' 
#' # show sensitivity
#' plot(out_sens)
#' 
cdesens <- function(seqg, trvar, medvar, rho =  seq(-0.9,0.9, by = 0.05)) {
  data <- seqg$model # model matrix
  
  rho <- sort(rho) # reorder if necessary
  acde.sens <- rep(NA, times = length(rho))
  acde.sens.se <- rep(NA, times = length(rho))
  n <- nrow(seqg$model)
  form.ytilde <- as.formula(paste0("y.tilde ~ ", paste0(names(coefficients(seqg))[-1], collapse = " + "))) # Ytilde ~ A + X
  form.m <- as.formula(paste0(medvar, " ~ ", paste0(names(coefficients(seqg))[-1], collapse = " + "))) # M ~ A + X
  
  res.y <- residuals(seqg$first_mod)
  res.m <- residuals(lm(form.m, data))
  
  rho.tilde <- cor(res.y, res.m) # scalar
  
  # for each value of rho, change mediator value
  m.fixed <- coef(seqg$first_mod)[medvar] - rho*sd(res.y)*sqrt((1-rho.tilde^2)/(1-rho^2))/sd(res.m)
  
  # calculate acde at each rho
  for (i in 1:length(rho)) {
    
    y.tilde <- seqg$model[[1]] - m.fixed[i]*(seqg$model[[medvar]])
    mf.i <- cbind(seqg$model, y.tilde)
    
    sens.direct  <- lm(form.ytilde, data = mf.i) # Ytilde ~ A + X
    acde.sens[i] <- coef(sens.direct)[trvar]
    
    sens.vcov <- seq.g.var(seqg$first_mod, sens.direct, med.vars = medvar)
    acde.sens.se[i] <- sqrt(sens.vcov[trvar, trvar])
  }
  
  out <- list(
    rho = rho,
    acde = acde.sens,
    se = acde.sens.se,
    direct = sens.direct
  )
  
  class(out) <- "cdesens"
  out
}


#' Plot output from cdesens
#' @param x output from \env{cdesens}
plot.cdesens <- function(x, level = 0.975, ...) {
  
  rho <- x$rho
  acde.sens <- x$acde
  ci.hi <- x$acde + qnorm(level) * x$se
  ci.lo <- x$acde - qnorm(level) * x$se
  
  plot(rho,
       acde.sens, 
       type = "n",
       ylim = range(c(ci.lo, ci.hi)),
       xlab = bquote("Correlation between mediator and outcome errors" ~~ (rho)),
       ylab = "Estimated ACDE", bty = "n", las = 1)
  polygon(x = c(rho, rev(rho)), y = c(ci.lo, rev(ci.hi)), col = "grey70", border = NA)
  lines(rho, acde.sens, lwd = 2)
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
}

