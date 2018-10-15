#' @export
summary.seqg <- function(object, ...) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (is.null(z$terms$direct)) {
    stop("invalid 'seqg' object: 'terms' incorrect")
  }
  r <- z$residuals
  n <- length(r)
  w <- z$weights
  if (is.null(z$vcov)) {
    se <- z$boots$acde.sd
    se.lab <- "Boot Std. Error"
  } else {
    se <- sqrt(diag(z$vcov))
    se.lab <- "Std. Error"
  }
  est <- z$coefficients
  tval <- est / se
  pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <- cbind(est, se, tval, pval)
  dimnames(ans$coefficients) <- list(names(z$coefficients), c("Estimate", se.lab, "t value", "Pr(>|t|)"))
  ans$df <- c(p, rdf, NCOL(z$qr$qr))
  class(ans) <- "summary.seqg"
  ans
}

#' @export
print.summary.seqg <- function(x, ...) {
  cat("\nt test of coefficients: \n\n")
  stats::printCoefmat(x$coefficients)
  cat("\n")
  invisible(x)
}
