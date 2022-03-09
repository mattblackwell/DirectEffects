form_engines <- c("lm", "logit", "multinom", "rlasso", "rlasso_logit")
matrix_engines <- c("lasso", "lasso_logit")

## for now hard code the egnines and the generalize later since its
## all under the hood.
## here we are going to pass an environment with the data already
## separated into a fit_data df for the fitting and a pred_data df for
## getting predictions. 
fit_model <- function(model, fit_env) {

  args <- as.list(model$args)  

  if (model$engine %in% form_engines) {
    args$formula <- model$formula
    args$data <- quote(fit_data)
    environment(args$formula) <- fit_env

    ## reorder args list to have formula first. 
    form_pos <- which(names(args) == "formula")
    args_pos <- seq_along(args)
    args_pos <- args_pos[c(form_pos, args_pos[-form_pos])]
    args <- args[args_pos]
  } else if (model$engine %in% matrix_engines) {
    mf <- model.frame(model$formula, data = fit_env$fit_data)
    fit_env$`.x` <- model.matrix(attr(mf, "terms"), data = mf)[, -1, drop = FALSE]
    fit_env$`.y` <- model.response(mf)
    args$x <- quote(`.x`)
    args$y <- quote(`.y`)
    if (var(fit_env$`.y`) == 0) fit_env$`.y` <- fit_env$`.y` + rnorm(length(fit_env$`.y`), 0, 0.0000001)
  } else {
    rlang::abort("engine not implemented")
  }
  
  if (model$engine == "lm") {
    fit_call <- rlang::call2("lm", !!!args, .ns = "stats")
  }
  if (model$engine == "lasso") {
    if (require(glmnet)) {
      fit_call <- rlang::call2("cv.glmnet", !!!args, .ns = "glmnet")
    } else {
      rlang::abort("glmnet package must be installed for `lasso` engine")
    }
  }
  if (model$engine == "rlasso") {
    if (require(hdm)) {
      fit_call <- rlang::call2("rlasso", !!!args, .ns = "hdm")
    } else {
      rlang::abort("hdm package must be installed for `rlasso` engine")
    }
  }
  if (model$engine == "rlasso_logit") {
    if (require(hdm)) {
      fit_call <- rlang::call2("rlassologit", !!!args, .ns = "hdm")
    } else {
      rlang::abort("hdm package must be installed for `rlasso` engine")
    }
  }
  if (model$engine == "lasso_logit") {
    args$family <- quote(binomial())
    if (require(glmnet)) {
      fit_call <- rlang::call2("cv.glmnet", !!!args, .ns = "glmnet")
    } else {
      rlang::abort("glmnet package must be installed for `lasso` engine")
    }
  }

  
  if (model$engine == "logit") {
    args$family <- quote(binomial(link = "logit"))
    fit_call <- rlang::call2("glm", !!!args, .ns = "stats")
  }

  if (model$engine == "multinom") {
    if (require(nnet)) {
      fit_call <- rlang::call2("multinom", !!!args, .ns = "nnet")
    } else {
      rlang::abort("nnet package must be installed for `multinom` engine")
    }
  }
  fit_env$fit <- rlang::eval_tidy(fit_call, env = fit_env)
  fit_env  
}

predict_model <- function(model, fit_env) {
  
  pred_args <- model$pred_args
  pred_args$object <- quote(fit)
  if (model$engine %in% c("lm", "rlasso")) {
    pred_args$newdata <- quote(pred_data)
    pred_call <- rlang::call2("predict", !!!pred_args, .ns = "stats")
  }
  if (model$engine %in% c("logit", "rlasso_logit")) {
    pred_args$newdata <- quote(pred_data)
    pred_args$type <- "response"
    pred_call <- rlang::call2("predict", !!!pred_args, .ns = "stats")
  }
  
  if (model$engine %in% c("lasso", "lasso_logit")) {
    if (require(glmnet)) {
      mf <- model.frame(model$formula[-2L], data = fit_env$pred_data)
      fit_env$`.x` <- model.matrix(attr(mf, "terms"), data = mf)[, -1, drop = FALSE]
      pred_args$newx <- quote(`.x`)
      pred_args$type <- "response"
      pred_args$s <- "lambda.1se"
      pred_call <- rlang::call2("predict", !!!pred_args)
    } else {
      rlang::abort("glmnet package must be installed for `lasso` engine")
    }
  }

  
  if (model$engine == "multinom") {
    if (require(nnet)) {
      pred_args$type <- "probs"
      pred_call <- rlang::call2("predict", !!!pred_args, .ns = "stats")
    } else {
      rlang::abort("nnet package must be installed for `multinom` engine")
    }
  }

  preds <- rlang::eval_tidy(pred_call, env = fit_env)

  ## preds for weights should return a matrix of predicted
  ## probabilities for each level of the outcome.
  ## TODO: weights for continuous? 
  if (model$engine %in% c("logit", "lasso_logit", "rlasso_logit")) {
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
  eval_vals <- get_eval_vals(object, data)
  eval_grid <- expand.grid(eval_vals)
  N_f <- nrow(A_fit)
  num_treat <- length(object$model_spec)

  ## move backward through blocks
  block_seq <- rev(seq_len(num_treat))

  model_fits <- out$model_fits


  
  if (object$has_outreg) {
    blipped_y <- create_blip_list(model_fits, eval_vals, Y)
    outreg_strata <- lapply(model_fits, function(x) colnames(x$outreg_pred))
  }
  
  fit_env$pred_data <- data
  tr_names <- unlist(lapply(
    object$model_spec,
    function(x) rlang::get_expr(x$treat)
  ))
  
  for (j in block_seq) {
    this_spec <- object$model_spec[[j]]
    
    if (object$has_ipw)  {
      if (!this_spec$treat_spec$separate) {

        fit_env$fit_data <- data[fit_rows,]
        fit_env <- fit_model(this_spec$treat_spec, fit_env)
        j_vals <- as.character(eval_vals[[j]])
        
        if (j > 1L) {
          past_vals <- apply(eval_grid[1:(j - 1)], paste0, collapse = "_")
          
          for (mp in seq_along(past_fit)) {
            strata <- paste0(past_fit[mp], "_", j_vals)
            curr_vals <- as.numeric(strsplit(past_vals[mp], "_")[[1]])
            
            fit_env$pred_data <- data
            for (k in seq_along(curr_vals)) {              
              fit_env$pred_data[[tr_names[[k]]]] <- curr_vals[k]
            }
            treat_fit <- predict_model(this_spec$treat_spec, fit_env)
            out$model_fits[[j]]$treat_pred[, nms] <- treat_fit[, j_vals]
            model_fits[[j]]$treat_pred[pred_rows, nms] <- treat_fit[pred_rows, j_vals]

          }
        }
        
        fit_env$pred_data <- data
        
        
        treat_fit <- predict_model(this_spec$treat_spec, fit_env)
        
        nms <- "pi"
        if (this_spec$treat_type == "categorical") {
          nms <- paste0(nms, "_", colnames(treat_fit))
        }
        out$model_fits[[j]]$treat_pred[, nms] <- treat_fit
        model_fits[[j]]$treat_pred[pred_rows, nms] <- treat_fit[pred_rows, ]
      } else {
        if (j > 1L) {
          past_fit <- interaction(A_fit[, 1:(j - 1), drop = FALSE], sep = "_")
          fit_splits <- split(seq_len(nrow(A_fit)), past_fit)
        } else {
          past_fit <- rep(0, times = N_f)
          fit_splits <- list(seq_len(nrow(A_fit)))
        }
        j_vals <- as.character(eval_vals[[j]])
        M <- length(fit_splits)
        
        ## TODO: add check about overlap?
        for (m in seq_len(M)) {
          fit_strata_rows <- fit_rows[fit_splits[[m]]]
          fit_env$fit_data <- data[fit_strata_rows, ]
          fit_env <- fit_model(this_spec$treat_spec, fit_env)
          fit_env$pred_data <- data
          treat_fit <- predict_model(this_spec$treat_spec, fit_env)
          nms <- names(fit_splits)[[m]]
          if (this_spec$treat_type == "categorical") {
            nms <- paste0(nms, ifelse(j > 1, "_", ""), j_vals)
          }
    
          model_fits[[j]]$treat_pred[pred_rows, nms] <- treat_fit[pred_rows, j_vals]
          out$model_fits[[j]]$treat_pred[, nms] <- treat_fit[, j_vals]
        }
      }
    }

    if (object$has_outreg) {
      ## TODO: figure out if we should allow class here
      if (this_spec$outreg_spec$engine_type == "class") {
        rlang::abort("only regression engines currently allowed for outcome models.")
      }
      
      upto_j <- 1L:j
      after_j <- setdiff(1L:num_treat, upto_j)
      before_j <- setdiff(1L:num_treat, j:num_treat)
      histories <- create_history_strings(eval_vals, 1L:num_treat)
      for (h in seq_along(histories)) {
        hist <- histories[h]
        upto_j_hist <- subset_history_string(hist, upto_j)
        aft_j_hist <- subset_history_string(hist, after_j)
        upto_j_hists <- create_history_factor(A_fit, upto_j)
       fit_env$pred_data <- data

        if (this_spec$outreg_spec$separate) {
          fit_splits <- which(upto_j_hists == upto_j_hist)
          these_rows <- fit_rows[fit_splits]
        } else {
          these_rows <- fit_rows
          fit_env$pred_data <- data
          curr_vals <- as.numeric(strsplit(upto_j_hist, "_")[[1]])          
          for (k in seq_len(j)) {
            fit_env$pred_data[[tr_names[[k]]]] <- curr_vals[k]
          }
        }
        fit_env$fit_data <- data[these_rows, ]
        if (j == num_treat) {
          fit_env$fit_data$`.de_y` <- blipped_y[[j]][these_rows, ]
        } else {
          fit_env$fit_data$`.de_y` <- blipped_y[[j]][these_rows, aft_j_hist]
        }
        fit_env <- fit_model(this_spec$outreg_spec, fit_env)
        outreg_fit <- predict_model(this_spec$outreg_spec, fit_env)
        out$model_fits[[j]]$outreg_pred[, hist] <- outreg_fit
        model_fits[[j]]$outreg_pred[pred_rows, hist] <- outreg_fit[pred_rows]
        if (j > 1L) {
          j_forward_hist <- subset_history_string(hist, j:num_treat)
          bef_j_hist <- subset_history_string(hist, before_j)
          bef_j_hists <- create_history_factor(A_fit, before_j)
          blip_splits <- split(seq_len(nrow(A_fit)), bef_j_hists)
          blip_rows <- fit_rows[which(bef_j_hists == bef_j_hist)]
          
          blips <- blip_down(object, out, Y, blipped_y, paths, j, hist, blip_rows)
          blipped_y[[j - 1]][blip_rows, j_forward_hist] <- blips
        }
      }
    }
  }
  model_fits
}

fit_treat <- function(spec, fit_env, data) {
  if (!treat_spec$separate) {
    treat_fit <- fit_model(this_spec$treat_spec, fit_env)
  } else {
    
    
  }
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
    treat_fit <- fit_model(this_spec$treat_spec, fit_env)
    nms <- paste0(
      names(fit_splits)[[m]],
      ifelse(j > 1, "_", ""),
      colnames(treat_fit))
  }  
  treat_fit
}


blip_down <- function(object, out, y, b_y, treat, j, strata, rows) {
  num_treat <- length(object$model_spec)
  y <- y[rows]
  if (class(object) %has% c("aipw", "did_aipw")) {
    p_scores <- get_ipw_preds(out, strata)
    p_scores <- p_scores[rows, j:num_treat, drop = FALSE]

    ## apply + do.call to ensure that weights is always a matrix
    weights <- apply(p_scores, 1, cumprod, simplify = FALSE)
    weights <- do.call(rbind, weights) 
    regs <- get_reg_preds(out, strata)
    regs <- regs[rows, j:num_treat, drop = FALSE]
    A <- get_path_inds(treat, strata)
    A <- cbind(1, A[rows, j:num_treat, drop = FALSE])
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
    blipped_y <- out$model_fits[[j]]$outreg_pred[rows, strata]
  }

  if (class(object) %has% "telescope_match") {
    if (j == num_treat) {
      blipped_y <- y[rows]
    } else {
      after_j_hist <- subset_history_string(strata, (j + 1):num_treat)
      blipped_y <- b_y[[j]][rows, strata]
    }    
    regs <- get_reg_preds(out, strata)
    regs <- regs[rows, j]
    A <- get_path_inds(treat, strata)
    A <- A[rows, j]
    matches <- out$match_out[[j]]$matches[rows]
    yhat_mr <- unlist(lapply(matches, function(x) mean(regs[x])))
    yhat_m <- unlist(lapply(matches, function(x) mean(blipped_y[x])))
    blipped_y[A == 0] <- yhat_m[A == 0] + (regs[A == 0] - yhat_mr[A == 0])
  }
  
  blipped_y
}
