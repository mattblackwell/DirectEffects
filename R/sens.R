#' Sensitivity of ACDE
#'
#' Estimate how ACDE varies by various levels of unmeasured confounding
#'
#' @param seqg Output from sequential_g
#' @param rho A numerical vector of correlations between errors to test for. The
#'  original model assumes \env{rho = 0}
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

cdesens <- function(seqg, rho =  seq(-0.9, 0.9, by = 0.05), boot = 50) {
  if (!inherits(seqg, what = "seqg")) {
    stop("object should be of class seqg, created from sequential_g()")
  }

  # model matrix
  data <- seqg$model

  rho <- sort(rho) # reorder if necessary

  # containers
  acde.sens <- matrix(NA, nrow = boot, ncol = length(rho)) # bootstrap samples as rows
  acde.sens.se <- rep(NA, times = length(rho))

  for (b in boot) {

    # create bootstrap sample
    b.index <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    data.b <- data[b.index, ]

    # identify treatment and mediator
    trvar <- attr(terms(formula(seqg$formula, lhs = 0, rhs = 1)), "term.labels")[1]
    medvar <- attr(terms(formula(seqg$formula, lhs = 0, rhs = 2)), "term.labels")
    if (length(medvar) > 1) stop("currently only handles one mediator variables")

    # formula
    form.A.X <- formula(seqg$formula, lhs = 0, rhs = 1)
    form.Ytilde <- update(form.A.X, Ytilde ~ .) # ytilde ~ A + X

    AX <- model.matrix(form.A.X, data.b)
    M <- data.b[, medvar, drop = TRUE]

    # residuals
    # epsilon.tilde.i.y: all variables in first model except medvar
    form.first.y.A.X <- update(seqg$first_mod$terms, paste0(". ~ . -", medvar))
    res.y <- residuals(lm(form.first.y.A.X, data = seqg$first_mod$model))

    # epsilon.tilde.i.m: residuals of mediation function
    res.m <- residuals(lm.fit(x = AX, y = M))

    rho.tilde <- cor(res.y, res.m)

    # for each value of rho, change mediator value
    rho.factor <- rho * sqrt((1 - rho.tilde^2) / (1 - rho^2))
    m.fixed <- coef(seqg$first_mod)[medvar] - sd(res.y) * rho.factor / sd(res.m)

    # calculate acde at each rho
    for (i in 1:length(rho)) {

      # create ytilde by blipping down with rho
      Ytilde <- seqg$model[[1]] - m.fixed[i] * (seqg$model[[medvar]])
      mf.i <- cbind(Ytilde, seqg$model)

      # Ytilde ~ A + X
      sens.direct.i <- lm(form.Ytilde, data = mf.i)
      acde.sens[b, i] <- coef(sens.direct.i)[trvar]
    }
  } # close bootstrap loop

  acde.sens.se <- apply(acde.sens, MARGIN = 2, sd)

  out <- list(
    rho = rho,
    acde = acde.sens,
    se = acde.sens.se
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
