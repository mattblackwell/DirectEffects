## for now hard code the egnines and the generalize later since its
## all under the hood.
## here we are going to pass an environment with the data already
## separated into a fit_data df for the fitting and a pred_data df for
## getting predictions. 
fit_model <- function(model, fit_env) {

  args <- as.list(model$args)  
  args$data <- quote(fit_data)
  args$formula <- model$formula
  environment(args$formula) <- fit_env
  
  pred_args <- model$pred_args
  pred_args$newdata <- quote(pred_data)
  pred_args$object <- quote(fit)
  if (model$engine == "lm") {
    fit_call <- rlang::call2("lm", !!!args, .ns = "stats")
    pred_call <- rlang::call2("predict", !!!pred_args, .ns = "stats")
  }

  if (model$engine == "logit") {
    args$family <- quote(binomial(link = "logit"))
    pred_args$type <- "response"
    fit_call <- rlang::call2("glm", !!!args, .ns = "stats")
    pred_call <- rlang::call2("predict", !!!pred_args, .ns = "stats")
  }
  fit_env$fit <- rlang::eval_tidy(fit_call, env = fit_env)
  preds <- rlang::eval_tidy(pred_call, env = fit_env)

  ## preds for weights should return a matrix of predicted
  ## probabilities for each level of the outcome.
  ## TODO: weights for continuous? 
  if (model$engine == "logit") {
    preds <- cbind(1 - preds, preds)
    colnames(preds) <- c("0", "1")
  }
  preds
  
}


fit_fold <- function(object, data, fit_rows, pred_rows) {

  fit_env <- rlang::env()  
  A_fit <- get_treat_df(object, data[fit_rows, ])
  A_pred <- get_treat_df(object, data[pred_rows, ])
  Y_fit <- get_outcome(object, data[fit_rows, ])
  
  N_f <- nrow(A_fit)
  N_p <- nrow(A_pred)
  num_treat <- length(object$model_spec)

  ## move backward through blocks
  block_seq <- rev(seq_len(num_treat))

  
  out <- list()
  if (object$has_ipw) out$ipw_pred <- make_pred_holder(A_pred, type = "ipw")
  if (object$has_outreg) out$outreg_pred <- make_pred_holder(A_pred, type = "outreg")
  fit_env$pred_data <- data[pred_rows, ]
  
  for (j in block_seq) {
    

    
    if (object$has_ipw) {
      if (j > 1) {
        past_fit <- interaction(A_fit[, 1:(j - 1), drop = FALSE], sep = "_")
        past_pred <- interaction(A_pred[, 1:(j - 1), drop = FALSE], sep = "_")
        fit_splits <- split(1:nrow(A_fit), past_fit)
      } else {
        past_fit <- rep(0, times = N_f)
        past_pred <- rep(0, times = N_p)
        fit_splits <- list(1:nrow(A_fit))
      }
      if (!identical(sort(unique(past_fit)), sort(unique(past_pred)))) {
        rlang::abort("lack of overlap not possible with `separate == TRUE`.")
      }
      M <- length(fit_splits)
      #### Fix this
      for (m in seq_len(M)) {
        fit_env$fit_data <- data[fit_rows[fit_splits[[m]]], ]
        treat_fit <- fit_model(object$model_spec[[j]]$treat_spec, fit_env)
        nms <- paste0(
          names(fit_splits)[[m]],
          ifelse(j > 1, "_", ""),
          colnames(treat_fit))
        out$ipw_pred[[j]][, nms] <- treat_fit
      }
    }

    if (object$has_outreg) {
      ## we split on current treatment too for outreg
      past_fit <- interaction(A_fit[, 1L:j, drop = FALSE], sep = "_")
      if (j < num_treat) {
        fut_fit <- interaction(A_fit[, (j + 1):num_treat, drop = FALSE], sep = "_")
      } else {
        fut_fit <- rep(0, times = N_f)
      }
      M_past <- length(unique(past_fit))
      M_fut <- length(unique(fut_fit))
      fit_splits <- split(seq_len(nrow(A_fit)), past_fit)
      for (mp in seq_len(M_past)) {
        for (mf in seq_len(M_fut)) {
          fit_env$fit_data <- data[fit_rows[fit_splits[[mp]]], ]
          if (j == num_treat) {
            fit_env$`.de_y` <- Y_fit[fit_splits[[mp]]]
            strata <- levels(past_fit)[mp]
          } else {
            strata <- paste0(levels(past_fit)[mp], "_", levels(fut_fit)[mf])
            old_fitted <- out$outreg_pred[[j + 1]][, strata]
            fit_env$`.de_y` <- old_fitted[fit_rows[fit_splits[[mp]]]]
          }
          
          outreg_fit <- fit_model(object$model_spec[[j]]$outreg_spec, fit_env)
          out$outreg_pred[[j]][, strata] <- outreg_fit
        }
        
      }
      

    }
  }
  out
}
