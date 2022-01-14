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






