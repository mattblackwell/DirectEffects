#' @export
print.seqg <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\n\nCoefficients:\n")
  print(x$coefficients)
  cat("\n")
  invisible(x)
}
