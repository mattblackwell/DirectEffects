#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

new_cde_estimator <- function(type, args, formula, model_spec) {
  
  out <- list(
    type = type,
    args = args,
    formula = formula,
    model_spec = model_spec
  )
  out <- add_class(out, c(type, "cde_estimator"))
  out
}


#' @export
print.cde_estimator <- function(x, ...) {
  cat("\n", x$type, "CDE Estimator\n")
  cat("---------------------------\n")
  tr_names <- unlist(
    lapply(x$model_spec, function(x) as.character(x$treat))
  )
  for (j in seq_along(tr_names)) {
    cat("Causal variable:", tr_names[j], "\n\n")
    if (length(x$model_spec[[j]]$treat_spec)) {
      cat("Treatment model:", deparse(x$model_spec[[j]]$treat_spec$formula), "\n")
      cat("Treatment engine:", x$model_spec[[j]]$treat_spec$engine, "\n\n")
    }
    if (length(x$model_spec[[j]]$outreg_spec)) {
      o_form <- x$model_spec[[j]]$outreg_spec$formula
      cat("Outcome model:", deparse(o_form[-2L]), "\n")
      cat("Outcome engine:", x$model_spec[[j]]$outreg_spec$engine, "\n")
    }
    cat("---------------------------\n")
  }
  invisible()
}


