#' @export
confint.seqg <- function(object, parm, level = 0.95, ...) {
  if (!is.null(object$vcov)) {
    return(confint.lm(object, parm, level))
  }
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]

  alpha <- 1 - level
  ltail <- 1 - alpha / 2
  utail <- alpha / 2

  percs <- c(round(alpha/2 * 100, 2), round((1 - alpha/2) * 100, 2))
  out <- matrix(NA_real_, nrow = length(parm), ncol = 2)
  rownames(out) <- parm
  colnames(out) <-  paste(percs, "%")

  preal <- parm[which(parm %in% pnames)]
  lshift <- apply(object$boots$acde[, preal, drop = FALSE], 2, quantile,
                  probs = ltail, na.rm = TRUE)
  ushift <- apply(object$boots$acde[, preal, drop = FALSE], 2, quantile,
                  probs = utail, na.rm = TRUE)

  out[preal, 1] <- 2 * object$coefficients[preal] - lshift
  out[preal ,2] <- 2 * object$coefficients[preal] - ushift
  return(out)
}
