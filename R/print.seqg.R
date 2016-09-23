print.seqg <- function(object) {
  cat("\nCall:\n")
  print(object$call)

  cat("\n\nCoefficients:\n")
  print(object$coefficients)
  cat("\n")
  invisible(object)
}
