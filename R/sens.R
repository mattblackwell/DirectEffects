#' Estimate sensitivity of ACDE estimates under varying levels of unobserved confounding
#'
#' Estimate how the Average Controlled Direct Effect varies by various levels of 
#' unobserved confounding. For each value of unmeasured confounding, summarized as
#' a correlation between residuals, \env{cdesens} computes the ACDE. Standard 
#' errors are computed by a simple bootstrap. 
#'
#' @param seqg Output from sequential_g. The function only supports specifications with one 
#'  mediator variable.
#' @param rho A numerical vector of correlations between errors to test for. The
#'  original model assumes \env{rho = 0}
#' @param boots Number of bootstrap replicates, defaults to 100.
#' @param verbose Whether to show progress and messages, defaults to \env{FALSE}
#'
#' @export
#'
#' @examples
#' data(civilwar)
#'
#' # First stage
#' fit_first <- lm(onset ~ ethfrac + warl + gdpenl + lpop +
#'                   ncontig + Oil + nwstate + polity2l + relfrac + instab,
#'                 data = civilwar)
#'
#' rows_use <- rownames(civilwar) %in% rownames(model.matrix(fit_first)) # listwise deletion
#'
#' # main formula: Y ~ A + X | M
#' form_main <- onset ~ ethfrac + lmtnest + ncontig + Oil | instab
#'
#' # estimate CDE
#' direct <- sequential_g(form_main, fit_first, data = civilwar, subset = rows_use)
#'
#' # sensitivity
#' out_sens <- cdesens(direct)
#'
#' # plot sensitivity
#' plot(out_sens)
#'

cdesens <- function(seqg, rho =  seq(-0.9, 0.9, by = 0.05), boots = 100, 
                    verbose = FALSE) {
  if (!inherits(seqg, what = "seqg")) {
    stop("object should be of class seqg, created from sequential_g()")
  }

  # model matrix
  data <- seqg$model

  rho <- sort(rho) # reorder if necessary

  # containers
  acde.sens <- matrix(NA, nrow = boots, ncol = length(rho)) # bootstrap samples as rows

  # identify treatment and mediator
  trvar <- attr(terms(formula(seqg$formula, lhs = 0, rhs = 1)), "term.labels")[1]
  medvar <- attr(terms(formula(seqg$formula, lhs = 0, rhs = 2)), "term.labels")
  if (length(medvar) > 1) stop("currently only handles one mediator variables")

  # formula
  # ~ A + X
  form.A.X <- formula(seqg$formula, lhs = 0, rhs = 1)
  # Ytilde ~ A + X. will create the Ytilde vector later.
  form.Ytilde <- update(form.A.X, Ytilde ~ .) 

  # start bootstrap (indexed by b)
  for (b in 1:boots) {
    # message at first boot and every 20
    if (verbose & (b == 1))
      cat("Starting Bootstrap estimation:", "\n")
    if (verbose & (b %% 20 == 0))
      cat(glue("sample {b} out of {boots}"), "\n")

    # create bootstrap sample
    b.index <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    data.b <- data[b.index, , drop = FALSE]

    # parts
    AX <- model.matrix(form.A.X, data.b)
    M <- data.b[, medvar, drop = TRUE]

    # re-fit first_mod call with new data
    first_mod.mm <- seqg$first_mod$model[b.index, , drop = FALSE] # first mod bootstrap
    first_mod <- update(seqg$first_mod, . ~ ., data = first_mod.mm)

    # residuals
    # epsilon.tilde.i.m: residuals of mediation function
    res.m <- residuals(lm.fit(x = AX, y = M))

    # epsilon.tilde.i.y: all variables in first model except medvar
    form.first.y.A.X <- update(seqg$first_mod$terms, paste0(". ~ . -", medvar))
    res.y <- residuals(lm(form.first.y.A.X, data = first_mod.mm))

    rho.tilde <- cor(res.y, res.m)

    # for each value of the vector rho, change mediator value
    rho.factor <- rho * sqrt((1 - rho.tilde^2) / (1 - rho^2))
    m.fixed <- coef(first_mod)[medvar] - sd(res.y) * rho.factor / sd(res.m)

    # calculate acde at each rho (indexed by r)
    for (r in 1:length(rho)) {

      # create ytilde by blipping down with rho
      Ytilde <- model.response(seqg$model)[b.index] - m.fixed[r] * (data.b[[medvar]])
      mf.r <- cbind(Ytilde, data.b)

      # run Ytilde ~ A + X
      sens.direct.r <- lm(form.Ytilde, data = mf.r)

      # save final estimate
      acde.sens[b, r] <- coef(sens.direct.r)[trvar]
      
    } # close rho loop
  } # close bootstrap loop

  # mean of bootstraps
  acede.means <- apply(acde.sens, MARGIN = 2, mean)

  # sd of bootstraps
  acde.se <- apply(acde.sens, MARGIN = 2, sd)

  # output
  out <- list(
    rho = rho,
    acde = acede.means,
    se = acde.se
  )

  class(out) <- "cdesens"
  out
}


#' Plot output from cdesens
#' @param x output from \env{cdesens}
#' @param level Level of confidence interval to plot
#' @param ... Other parameters to pass on to \env{plot()}
#' @export
plot.cdesens <- function(x, level = 0.95, ...) {
  rho <- x$rho
  acde.sens <- x$acde
  ci.hi <- x$acde + qnorm(1 - (1 - level) / 2) * x$se
  ci.lo <- x$acde - qnorm(1 - (1 - level) / 2) * x$se

  plot(rho,
    acde.sens,
    type = "n",
    ylim = range(c(ci.lo, ci.hi)),
    xlab = bquote("Correlation between mediator and outcome errors" ~ ~ (rho)),
    ylab = "Estimated ACDE", bty = "n", las = 1, ...
  )
  polygon(x = c(rho, rev(rho)), y = c(ci.lo, rev(ci.hi)), col = "grey70", border = NA)
  lines(rho, acde.sens, lwd = 2)
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
}
