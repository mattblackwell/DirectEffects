#' Computes standard errors and p-values of DirectEffects estimates
#'
#'
#' @param object An object of class \env{seqg}, computed by
#'   \code{\link{sequential_g}}.
#' @param ... additional arguments affecting the summary produced.
#'
#' @export
#'
summary.seqg <- function(object, ...) {
  z1 <- object$first_mod
  z2 <- object
  p1 <- z1$rank
  p2 <- z2$rank
  rdf <- z2$df.residual
  if (is.null(z2$terms$X)) {
    stop("invalid 'seqg' object: 'terms' incorrect")
  }
  r1 <- z1$residuals
  r2 <- z2$residuals
  n <- length(r2)
  w1 <- z1$weights
  w2 <- z2$weights
  mv <- attr(z2$terms$M, "term.labels")

  p1v <- 1L:p1
  Qr1 <- z1$qr
  R1 <- chol2inv(Qr1$qr[p1v, p1v, drop = FALSE])
  efun1 <- if (is.null(w1)) z1$X * r1 else w1 * z1$X * r1

  p2v <- 1L:p2
  Qr2 <- z2$qr
  R2 <- chol2inv(Qr2$qr[p2v, p2v, drop = FALSE])
  efun2 <- if (is.null(w2)) z2$X * r2 else w2 * z2$X * r2
  X <- z2$X[, !z2$aliased]
  XZM <- z1$XZM[, !z1$aliased]
  Fhat <- crossprod(X, XZM) / n
  Fhat[, !(colnames(XZM) %in% mv)] <- 0

  ghat <- t(efun2) + Fhat %*% R1 %*% t(efun1)
  meat <- crossprod(t(ghat))
  dfc <- (n / (n - p2))
  vcov <- dfc * (R2 %*% meat %*% R2)

  se <- sqrt(diag(vcov))
  est <- z2$coefficients
  tval <- est / se
  pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
  ans <- z2[c("call", "terms", if (!is.null(z2$weights)) "weights")]
  ans$residuals <- r2
  ans$coefficients <- cbind(est, se, tval, pval)
  dimnames(ans$coefficients) <- list(names(z2$coefficients), c("Estimate", "Std. Err.", "t value", "Pr(>|t|)"))
  ans$df <- c(p2, rdf, NCOL(Qr2$qr))
  ans$vcov <- vcov
  dimnames(vcov) <- dimnames(ans$coefficients)[c(1,1)]
  class(ans) <- "summary.seqg"
  ans
}


#' Summary of DirectEffect Bootstrap Estimates
#'
#' @param object An output of class \code{seqg} estimated by \code{\link{boots_g}}.
#' @param level level of intervals to estimate. Defaults to 0.95
#' @param ... additional arguments affecting the summary produced.
#'
#' @export
#' @importFrom glue glue
summary.seqgboots <- function(object, level = 0.95, ...) {

  # lower and upper percentile
  lo <- (1 - level)/2
  hi <- 1 - (1 - level)/2

  # column stats
  est <- colMeans(object)
  se <- apply(object, 2, sd)
  tval <- est / se
  plo <- apply(object, 2, function(x) quantile(x, lo))
  phi <- apply(object, 2, function(x) quantile(x, hi))

  # bind
  coefs <- cbind(est, se, tval, plo, phi)
  dimnames(coefs) <- list(colnames(object),
                          c("Estimate",
                            "Std. Err.",
                            "t value",
                            glue("{round(lo*100, 1)} %"),
                            glue("{round(hi*100, 1)} %")))

  class(coefs) <- "summary.seqgboots"
  coefs
}


#' @export
print.summary.seqg <- function(x, ...) {
  cat("\nt test of coefficients: \n\n")
  stats::printCoefmat(x$coefficients)
  cat("\n")
  invisible(x)
}



#' @export
print.summary.seqgboots <- function(x, ...) {
  cat("\nSummary of bootstrapped coefficients: \n\n")
  stats::printCoefmat(x)
  cat("\n")
  invisible(x)
}
