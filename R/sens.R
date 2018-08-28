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
#' @param bootstrap character of c("none", "standard", "block"), indicating whether to
#' include bootstrap standard errors or block bootstrap. Default is "none".
#' @param boots_n Number of bootstrap replicates, defaults to 100.
#' @param verbose Whether to show progress and messages, defaults to \env{FALSE}
#'
#' @export
#'
#' @importFrom stats formula lm cor sd
#' @import glue
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
#' # sensitivity
#' out_sens <- cdesens(direct)
#'
#' # plot sensitivity
#' plot(out_sens)
#'

cdesens <- function(seqg, rho =  seq(-0.9, 0.9, by = 0.05),
                    bootstrap = c("none", "standard", "block"),
                    boots_n = 1000, verbose = FALSE) {
  if (!inherits(seqg, what = "seqg")) {
    stop("object should be of class seqg, created from sequential_g()")
  }

  if (missing(bootstrap)) bootstrap <- "none"
  # model matrix
  data <- seqg$model

  rho <- sort(rho) # reorder if necessary


  # identify treatment and mediator
  trvar <- attr(terms(formula(seqg$formula, lhs = 0, rhs = 1)), "term.labels")[1]
  medvar <- attr(terms(formula(seqg$formula, lhs = 0, rhs = 2)), "term.labels")
  if (length(medvar) > 1) stop("currently only handles one mediator variables")

  # formula
  # ~ A + X
  form.A.X <- formula(seqg$formula, lhs = 0, rhs = 1)
  # Ytilde ~ A + X. will create the Ytilde vector later.
  form.Ytilde <- update(form.A.X, Ytilde ~ .)

  if (bootstrap != "none") {
    if (bootstrap == "block") warning("block bootstrap not implemented yet, using standard...")

    # containers
    acde.sens <- matrix(NA, nrow = boots_n, ncol = length(rho)) # bootstrap samples as rows

    # start bootstrap (indexed by b)
    for (b in 1:boots_n) {
      # message at first boot and every 20
      if (verbose & (b == 1))
        cat("Starting Bootstrap estimation:", "\n")
      if (verbose & (b %% 20 == 0))
        cat(glue::glue("sample {b} out of {boots_n}"), "\n")

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
    acde.means <- apply(acde.sens, MARGIN = 2, mean)

    # sd of bootstraps
    acde.se <- apply(acde.sens, MARGIN = 2, sd)

  } # close bootstrap if

  if (bootstrap == "none") {
    # parts
    AX <- model.matrix(form.A.X, data)
    M <- data[, medvar, drop = TRUE]

    # residuals
    # epsilon.tilde.i.m: residuals of mediation function
    res.m <- residuals(lm.fit(x = AX, y = M))

    # epsilon.tilde.i.y: all variables in first model except medvar
    form.first.y.A.X <- update(seqg$first_mod$terms, paste0(". ~ . -", medvar))
    res.y <- residuals(lm(form.first.y.A.X, data = stats::model.frame(seqg$first_mod)))

    rho.tilde <- cor(res.y, res.m)

    # for each value of the vector rho, change mediator value
    rho.factor <- rho * sqrt((1 - rho.tilde^2) / (1 - rho^2))
    m.fixed <- coef(seqg$first_mod)[medvar] - sd(res.y) * rho.factor / sd(res.m)

    acde.means <- acde.se <- rep(NA, times = length(rho))
    # calculate acde at each rho (indexed by r)
    for (r in 1:length(rho)) {

      # create ytilde by blipping down with rho
      Ytilde <- model.response(seqg$model) - m.fixed[r] * (data[[medvar]])
      mf.r <- cbind(Ytilde, data)


      # run Ytilde ~ A + X
      sens.direct.r <- lm(form.Ytilde, data = mf.r)
      vcv <- seq.g.vcov(seqg$first_mod, sens.direct.r,
                        X1 = model.matrix(seqg$first_mod), X2 = AX, medvar)

      # save final estimate
      acde.means[r] <- coef(sens.direct.r)[trvar]
      acde.se[r] <- sqrt(diag(vcv)[which(colnames(AX) == trvar)])

    } # close rho loop
  } # close bootstrap=="none" if


  # output
  out <- list(
    rho = rho,
    acde = acde.means,
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
#' @importFrom graphics plot lines polygon abline
#' @importFrom stats qnorm
plot.cdesens <- function(x, level = 0.95, xlim = NULL, ylim = NULL,
                         xlab = NULL, ylab = "Estimated ACDE", bty = "n",
                         ci.col = "grey70", col = "black", lwd = 2,
                         ref.lines = TRUE, ...) {
  rho <- x$rho
  acde.sens <- x$acde
  ci.hi <- x$acde + qnorm(1 - (1 - level) / 2) * x$se
  ci.lo <- x$acde - qnorm(1 - (1 - level) / 2) * x$se

  if (is.null(xlim)) xlim <- range(rho)
  if (is.null(ylim)) ylim <- range(c(ci.lo, ci.hi))
  if (is.null(xlab)) xlab <- bquote("Correlation between mediator and outcome errors" ~ ~ (rho))
  if (is.null(ylab)) ylab <- "Estimated ACDE"

  plot.default(NA, NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       bty = bty, ...)
  polygon(x = c(rho, rev(rho)), y = c(ci.lo, rev(ci.hi)), col = ci.col, border = NA)
  lines(rho, acde.sens, lwd = lwd, col = col)

  if (ref.lines) {
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
  }

  invisible(NULL)
}
