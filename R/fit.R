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


fit_fold <- function(object, data, fit_rows, pred_rows, out) {

  fit_env <- rlang::env()  
  A_fit <- get_treat_df(object, data[fit_rows, ])
  A <- get_treat_df(object, data)
  paths <- interaction(A, sep = "_")
  Y <- get_outcome(object, data)
  
  N_f <- nrow(A_fit)
  num_treat <- length(object$model_spec)

  ## move backward through blocks
  block_seq <- rev(seq_len(num_treat))

  
  if (object$has_ipw) out$ipw_pred <- make_pred_holder(A, type = "ipw")
  if (object$has_outreg) {
    out$outreg_pred <- make_pred_holder(A, type = "outreg")
    blipped_y <- make_pred_holder(A, type = "outreg")
  }
  fit_env$pred_data <- data
  
  for (j in block_seq) {
    if (object$has_ipw)  {
      if (j > 1) {
        past_fit <- interaction(A_fit[, 1:(j - 1), drop = FALSE], sep = "_")
        fit_splits <- split(seq_len(nrow(A_fit)), past_fit)
      } else {
        past_fit <- rep(0, times = N_f)
        fit_splits <- list(seq_len(nrow(A_fit)))
      }
      M <- length(fit_splits)      
      
      ## TODO: add check about overlap?
      for (m in seq_len(M)) {
        fit_strata_rows <- fit_rows[fit_splits[[m]]]
        fit_env$fit_data <- data[fit_strata_rows, ]
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
          if (j == num_treat) {
            strata <- levels(past_fit)[mp]
          } else {
            strata <- paste0(levels(past_fit)[mp], "_", levels(fut_fit)[mf])
          }
          
          fit_strata_rows <- fit_rows[fit_splits[[mp]]]
          fit_env$fit_data <- data[fit_strata_rows, ]
          blipped_y[[j]][fit_strata_rows, strata] <- blip_down(object, out, Y, blipped_y, paths, j, strata, fit_strata_rows)
          
          fit_env$`.de_y` <- blipped_y[[j]][fit_strata_rows, strata]
          outreg_fit <- fit_model(object$model_spec[[j]]$outreg_spec, fit_env)
          out$outreg_pred[[j]][, strata] <- outreg_fit
        }
        
      }
      
    }
  }
  out
}


blip_down <- function(object, out, y, b_y, treat, j, strata, rows) {
  num_treat <- length(object$model_spec)
  y <- y[rows]
  if (j == num_treat) {
    return(y)
  }

  if (class(object) %has% c("aipw", "did_aipw")) {
    p_scores <- get_ipw_preds(out, strata)
    p_scores <- p_scores[rows, (j + 1):num_treat, drop = FALSE]

    ## apply + do.call to ensure that weights is always a matrix
    weights <- apply(p_scores, 1, cumprod, simplify = FALSE)
    weights <- do.call(rbind, weights) 
    regs <- get_reg_preds(out, strata)
    regs <- regs[rows, (j + 1):num_treat, drop = FALSE]
    A <- get_path_inds(treat, strata)
    A <- A[rows, j:num_treat, drop = FALSE]
    eps <- cbind(regs, y) - cbind(0, regs)
    weights <- cbind(1, weights)
    if (length(object$args$trim)) {
      weights <- winsorize_matrix(weights, object$args$trim)
    }
    if (object$args$aipw_blip) {
      blipped_y <- rowSums(A * eps / weights)
    } else {
      blipped_y <- regs[, 1L]
    }
    
  }

  if (class(object) %has% "reg_impute") {
    blipped_y <- out$outreg_pred[[j + 1]][rows, strata]
  }

  if (class(object) %has% "telescope_match") {
    blipped_y <- b_y[[j + 1]][rows, strata]
    regs <- get_reg_preds(out, strata)
    regs <- regs[rows, j + 1]
    A <- get_path_inds(treat, strata)
    A <- A[rows, j + 1]
    matches <- out$match_out[[j + 1]]$matches[rows]
    yhat_mr <- unlist(lapply(matches, function(x) mean(regs[x])))
    yhat_m <- unlist(lapply(matches, function(x) mean(blipped_y[x])))
    blipped_y[A == 0] <- yhat_m[A == 0] + (regs[A == 0] - yhat_mr[A == 0])
  }
  
  blipped_y
}
