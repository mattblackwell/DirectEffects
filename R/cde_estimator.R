#' @importFrom generics glance
#' @export
generics::glance

new_cde_estimator <- function(type, args, model_spec) {
  
  check_cde_estimator(type, model_spec)

  out <- list(
    type = type,
    args = args,
    model_spec = model_spec
  )
  out <- add_class(out, c(type, "cde_estimator"))
  out
}


#' @export
print.cde_estimator <- function(x, ...) {
  cat("\n", x$type, "CDE Estimator\n")
  cat("---------------------------\n")
  cat("Causal variables:", paste0(x$treat_names, collapse = ", "), "\n")
  
  cat("\nEstimated Effects:\n")
  print(x$estimates[, c("term", "active", "control", "estimate")])
  cat("\n")
  invisible()
}

#' @export
tidy.cde_estimator <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  est <- x$estimates
  
  est$statistic <- est$estimate / est$std.error
  est$p.value <- 2 * stats::pt(abs(est$statistic), df = est$DF, lower = FALSE)
  est$df <- est$DF
  
  if (conf.int) {
    alpha <- (1 - conf.level) / 2
    c_lo <- stats::qt(alpha, df = est$DF)
    c_hi <- stats::qt(1 - alpha, df = est$DF)
    ci <- cbind(
      est$estimate + est$std.error * c_lo,
      est$estimate + est$std.error * c_hi
    )
    colnames(ci) <- c("conf.low", "conf.high")
    est <- cbind(est, ci)
    col_names <- c(
      "term",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "conf.low",
      "conf.high",
      "df"
    )
  } else {
    col_names <- c(
      "term",
      "estimate",
      "std.error",
      "statistic",
      "p.value",
      "df"
    )
  }

  est$active <- paste0("(", gsub("_", ", ", est$active), ")")
  est$control <- paste0("(", gsub("_", ", ", est$control), ")")
  est$term <- paste0(est$term, " [", est$active, " vs. ", est$control, "]")
  est <- est[, col_names]
  est
}

#' @export 
summary.cde_estimator <- function(object, ...) {
  
  est <- tidy(object)
  rownames(est) <- est$term
  est <- est[, colnames(est) != "term"]
  colnames(est) <- c(
    "Estimate",
    "Std. Error",
    "t value",
    "Pr(>|t|)",
    "CI Lower",
    "CI Upper",
    "DF"
  )
  out <- list()
  out$type <- object$type
  out$crossfit <- object$crossfit
  out$n_folds <- object$n_folds
  out$treat_names <- object$treat_names
  class(out) <- "summary.cde_estimator"
  out$coefficients <- est
  return(out)
}

#' @export
print.summary.cde_estimator <-  function(x,
                                         digits = max(3L, getOption("digits") - 3L),
                                         signif.stars = getOption("show.signif.stars"),
                                         ...) {
  cat("\n", x$type, "CDE Estimator\n")
  cat("---------------------------\n")
  cat("Causal variables:", paste0(x$treat_names, collapse = ", "), "\n")
  cat("Cross-fit:", x$crossfit, "\n")
  if (x$crossfit) cat("Number of folds:", x$n_folds, "\n")
  
  cat("\nEstimated Effects:\n")
  print(
    x$coefficients,
    digits = digits, quote = FALSE, right = TRUE
  )
  cat("\n")

}
