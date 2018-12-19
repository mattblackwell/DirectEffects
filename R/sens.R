#' Estimate sensitivity of ACDE estimates under varying levels of unobserved confounding
#'
#' Estimate how the Average Controlled Direct Effect varies by various levels of
#' unobserved confounding. For each value of unmeasured confounding, summarized as
#' a correlation between residuals, \env{cdesens} computes the ACDE. Standard
#' errors are computed by a simple bootstrap.
#'
#' @param seqg Output from sequential_g. The function only supports specifications with one
#'  mediator variable.
#' @param var A character indicating the name of the variable for
#'   which the estimated ACDE is being evaluated.
#' @param rho A numerical vector of correlations between errors to test for. The
#'  original model assumes \env{rho = 0}
#' @param bootstrap character of c("none", "standard", "block"), indicating whether to
#' include bootstrap standard errors or block bootstrap. Default is "none".
#' @param boots_n Number of bootstrap replicates, defaults to 100.
#' @param verbose Whether to show progress and messages, defaults to
#'   \env{FALSE}
#' @param ... Other parameters to pass on to \env{lm.fit()} when
#'   refitting the model
#'
#' @export
#'
#' @importFrom stats formula lm cor sd
#' @import glue
#' @examples
#' data(civilwar)
#'
#'
#' # main formula: Y ~ A + X | Z | M
#' form_main <- onset ~ ethfrac + lmtnest + ncontig + Oil | warl +
#'   gdpenl + lpop + polity2l + relfrac | instab
#'
#' # estimate CDE
#' direct <- sequential_g(form_main, data = civilwar)
#'
#' # sensitivity
#' out_sens <- cdesens(direct, var = "ethfrac")
#'
#' # plot sensitivity
#' plot(out_sens)
#'

cdesens <- function(seqg, var, rho = seq(-0.9, 0.9, by = 0.05),
                    bootstrap = c("none", "standard", "block"),
                    boots_n = 1000, verbose = FALSE, ...) {
  if (!inherits(seqg, what = "seqg")) {
    stop("object should be of class seqg, created from sequential_g()")
  }
  ## save copy of object to use below

  z <- seqg
  bootstrap <- match.arg(bootstrap)
  # model matrix
  data <- seqg$model

  rho <- sort(rho) # reorder if necessary

  # identify treatment and mediator
  xnames <- attr(seqg$terms$X, "term.labels")
  mnames <- attr(seqg$terms$M, "term.labels")

  if (length(mnames) > 1) stop("currently only handles one mediator variables")
  if (!(var %in% xnames)) stop("'var' not in the set of baseline variables")

  if (bootstrap != "none") {
    if (bootstrap == "block") warning("block bootstrap not implemented yet, using standard...")

    # containers
    acde.sens <- matrix(NA, nrow = boots_n, ncol = length(rho)) # bootstrap samples as rows

    # start bootstrap (indexed by b)
    for (b in 1:boots_n) {
      # message at first boot and every 20
      if (verbose & (b == 1)) {
        cat("Starting Bootstrap estimation:", "\n")
      }
      if (verbose & (b %% 20 == 0)) {
        cat(glue::glue("sample {b} out of {boots_n}"), "\n")
      }

      # create bootstrap sample
      b.index <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
      data.b <- data[b.index,]

      # relevant data matrices
      X <- seqg$X[b.index, ]
      XZM <- seqg$first_mod$XZM[b.index, ]
      M <- seqg$M[b.index, ]
      XZ <- XZM[, !(colnames(XZM) %in% mnames), drop = FALSE]
      Y <- model.response(data.b)
      w <- as.vector(model.weights(data.b))
      offset <- as.vector(model.offset(data.b))

      # first-stage
      if (is.null(w)) {
        out.first <- lm.fit(XZM, Y, offset = offset, ...)
      } else {
        out.first <- lm.wfit(XZM, Y, w, offset = offset, ...)
      }
      bcoefs <- out.first$coefficients[mnames]

      # epsilon.tilde.i.m: residuals of mediation function
      res.m <- residuals(lm.fit(x = X, y = M))
      res.y <- residuals(lm.fit(x = XZ, y = Y))
      rho.tilde <- cor(res.y, res.m)

      # for each value of the vector rho, change mediator value
      rho.factor <- rho * sqrt((1 - rho.tilde^2) / (1 - rho^2))
      m.fixed <- bcoefs - sd(res.y) * rho.factor / sd(res.m)

      # calculate acde at each rho (indexed by r)
      for (r in 1:length(rho)) {
        # create ytilde by blipping down with rho
        Ytilde <- Y - m.fixed[r] * M
        if (is.null(w)) {
          out <- lm.fit(X, Ytilde, offset = offset, ...)
        } else {
          out <- lm.wfit(X, Ytilde, w, offset = offset, ...)
        }
        # save final estimate
        acde.sens[b, r] <- out$coefficients[var]
      } # close rho loop
    } # close bootstrap loop

    # mean of bootstraps
    acde.means <- apply(acde.sens, MARGIN = 2, mean)

    # sd of bootstraps
    acde.se <- apply(acde.sens, MARGIN = 2, sd)
  } # close bootstrap if

  if (bootstrap == "none") {
    # residuals
    # epsilon.tilde.i.m: residuals of mediation function
    XZM <- seqg$first_mod$XZM
    XZ <- XZM[, !(colnames(XZM) %in% mnames), drop = FALSE]
    w <- as.vector(model.weights(seqg$model))
    offset <- as.vector(model.offset(seqg$model))

    res.m <- residuals(lm.fit(x = seqg$X, y = seqg$M))
    res.y <- residuals(lm.fit(x = XZ, y = seqg$Y))
    rho.tilde <- cor(res.y, res.m)

    # for each value of the vector rho, change mediator value
    rho.factor <- rho * sqrt((1 - rho.tilde^2) / (1 - rho^2))
    m.fixed <- coef(seqg$first_mod)[mnames] -
      sd(res.y) * rho.factor / sd(res.m)

    acde.means <- acde.se <- rep(NA, times = length(rho))
    # calculate acde at each rho (indexed by r)
    for (r in 1:length(rho)) {

      # create ytilde by blipping down with rho
      Ytilde <- seqg$Y - m.fixed[r] * seqg$M
      if (is.null(w)) {
        out <- lm.fit(seqg$X, Ytilde, offset = offset, ...)
      } else {
        out <- lm.wfit(seqg$X, Ytilde, w, offset = offset, ...)
      }

      # update the saved copy of the object and recalculate vcov
      z$coefficients <- out$coefficients
      z$residuals <- out$residuals
      z$qr <- out$qr
      vcv <- vcov.seqg(z)

      # save final estimate
      acde.means[r] <- out$coefficients[var]
      acde.se[r] <- sqrt(diag(vcv)[which(colnames(seqg$X) == var)])
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
#' @param level level of confidence interval to plot
#' @param xlim the x limits (x1, x2) of the plot for the sensitivity
#'   analysis parameter, rho. Default is to use the range of
#'   \env{rho}.
#' @param ylim the y limits of the plot for the estimated CDEs.
#'   Default is to show the all of the confidence intervals.
#' @param xlab label for the x axis.
#' @param ylab label for the y axis.
#' @param bty a character string which determined the type of box
#'   which is drawn about plots. Defaults to not drawing a box. See
#'   \link{par} for more information.
#' @param col color for the line indicating the point estimates of the
#'   bias-adjusted ACDE.
#' @param lwd line width for the line indicating the point estimates of the
#'   bias-adjusted ACDE.
#' @param ci.col color for the polygon that shows the confidence
#'   intervals.
#' @param ref.lines a logical indicating whether horizontal and
#'   vertical lines at 0 should be plotted.
#' @param ... Other parameters to pass on to \env{plot()}
#' @export
#' @importFrom graphics plot.default lines polygon abline
#' @importFrom stats qnorm
plot.cdesens <- function(x, level = 0.95, xlim = NULL, ylim = NULL,
                         xlab = NULL, ylab = "Estimated ACDE", bty = "n",
                         col = "black", lwd = 2, ci.col = "grey70",
                         ref.lines = TRUE, ...) {
  rho <- x$rho
  acde.sens <- x$acde
  ci.hi <- x$acde + qnorm(1 - (1 - level) / 2) * x$se
  ci.lo <- x$acde - qnorm(1 - (1 - level) / 2) * x$se

  if (is.null(xlim)) xlim <- range(rho)
  if (is.null(ylim)) ylim <- range(c(ci.lo, ci.hi))
  if (is.null(xlab)) xlab <- bquote("Correlation between mediator and outcome errors" ~ ~(rho))
  if (is.null(ylab)) ylab <- "Estimated ACDE"

  plot.default(NA, NA,
    xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
    bty = bty, ...
  )
  polygon(x = c(rho, rev(rho)), y = c(ci.lo, rev(ci.hi)), col = ci.col, border = NA)
  lines(rho, acde.sens, lwd = lwd, col = col)

  if (ref.lines) {
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
  }

  invisible(NULL)
}
