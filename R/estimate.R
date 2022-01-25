#' @export
estimate <- function(object, formula, data, crossfit = TRUE, n_folds, fold_seed = NULL) {

  N <- nrow(data)
  J <- length(object$model_spec)

  A <- get_treat_df(object, data)
  object$has_ipw <- object$type %in% c("ipw", "aipw")
  object$has_outreg <- object$type %in% c("reg_impute", "aipw", "telescope_match")
  object$has_match <- object$type %in% c("telescope_match")
  
  if (missing(n_folds)) {
    if (crossfit) {
      n_folds <- 5L
    } else {
      n_folds <- 1L
    }
  }

  object$outcome <- formula[[2L]]

  

    
  out <- list()
  if (object$has_ipw) out$ipw_pred <- make_pred_holder(A, type = "ipw")
  if (object$has_outreg) out$outreg_pred <- make_pred_holder(A, type = "outreg")


  ## match outside the crossfit loop
  if (object$has_match) {
    if (crossfit) {
      rlang::warn("cannot use crossfitting with matching; setting `crossfit = FALSE`")
      crossfit <- FALSE
    }
    out$match_out <- t_match(object, data)
  }


  if (crossfit) {
    if (length(fold_seed)) set.seed(fold_seed)
    fold_size <- N / n_folds
    f_split <- ceiling(seq_len(N) / fold_size)
    folds <- split(sample(seq_len(N)), f_split)
  } else {
    folds <- list(seq_len(N))
  }  
  out$n_folds <- n_folds
  out$folds <- folds

  for (k in seq_len(n_folds)) {

    if (crossfit) {
      fit_rows <- unlist(folds[-k])
      pred_rows <- unlist(folds[k])
    } else {
      fit_rows <- 1L:N
      pred_rows <- 1L:N
    }
    
    res <- fit_fold(object, data, fit_rows, pred_rows, out)

    ## fill in prediction matrices
    for (j in seq_len(J)) {
      if (object$has_ipw) {
        out$ipw_pred[[j]][pred_rows, ] <- res$ipw_pred[[j]][pred_rows, ]
      }
      if (object$has_outreg) {
        out$outreg_pred[[j]][pred_rows, ] <- res$outreg_pred[[j]][pred_rows, ]
      }
    }

  }
  out <- estimate_cde(object, formula, data, out)
  
  out
}



estimate_cde <- function(object, formula, data, out) {
  eff_vars <- all.vars(formula)[-1L]
  tr_names <- unlist(lapply(object$model_spec, function(x) as.character(x$treat[[2L]])))
  y <- get_outcome(object, data)
  if (!all(eff_vars %in% tr_names)) {
    rlang::abort("Unspecified treatment included in `estimate()` formula.")
  }

  eff_pos <- match(eff_vars, tr_names)
  eff_vars <- eff_vars[eff_pos]
  eff_pos <- sort(eff_pos)
  A <- get_treat_df(object, data)
  N <- nrow(A)
  num_treat <- length(tr_names)
  out$estimates <- empty_est_tab()
  paths <- interaction(A, sep = "_")
  path_levs <- levels(paths)
  path_splits <- split(1L:N, paths)
  
  if (class(object) %has% "ipw") {

    psi <- vector("numeric", length = length(path_levs))
    psi_sq <- vector("numeric", length = length(path_levs))
    names(psi) <- names(psi_sq) <- path_levs
    for (p in path_levs) {
      p_rows <- path_splits[[p]]
      p_scores <- get_ipw_preds(out, p)
      weights <- apply(p_scores, 1, prod)
      if (object$args$hajek) {
        N_p <-  sum(1 / weights[p_rows])
      } else {
        N_p <- length(p_rows)
      }
      psi[p] <- sum(y[p_rows] / weights[p_rows]) / N_p
      psi_sq[p] <- (N / N_p ^ 2) * sum(y[p_rows] ^ 2 / weights[p_rows] ^ 2)      
    }

    for (e in seq_along(eff_vars)) {
      j <- eff_pos[e]
      j_levs <- unique(A[, j])
      out$estimates <- rbind(
        out$estimates,
        compute_ipw(j, j_levs, y, paths, out, object$args, eff_vars[e])
      )
    }
  }

  if (class(object) %has% "reg_impute") {
    for (e in seq_along(eff_vars)) {
      j <- eff_pos[e]
      j_levs <- unique(A[, j])
      out$estimates <- rbind(
        out$estimates,
        compute_reg_impute(j, j_levs, y, paths, out$outreg_pred[[e]], eff_vars[e])
      )
    }
  }

  if (class(object) %has% "aipw") {
    for (e in seq_along(eff_vars)) {
      j <- eff_pos[e]
      j_levs <- unique(A[, j])
      out$estimates <- rbind(
        out$estimates,
        compute_aipw(j, j_levs, y, paths, out, object$args, eff_vars[e])
      )
    }
  }

  if (class(object) %has% "telescope_match") {
    dfs <- unlist(lapply(object$model_spec, function(x) x$outreg_spec$df))
    for (e in seq_along(eff_vars)) {
      j <- eff_pos[e]
      j_levs <- unique(A[, j])
      out$estimates <- rbind(
        out$estimates,
        compute_telescope_match(j, j_levs, y, paths, out, object$args, eff_vars[e], dfs)
      )
    }
  }

  out

}
