#' @export
vcov.seqg <- function(object, ...) {
  out.vcov <- matrix(NA, ncol(object$X), ncol(object$X))
  dimnames(out.vcov) <- list(names(object$coefficients),
                             names(object$coefficients))
  R <- summary(object)$vcov
  out.vcov[!object$aliased, !object$aliased] <- R
  out.vcov
}
