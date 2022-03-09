#' @export
estimate <- function(object, formula, data, subset, crossfit = TRUE, n_folds, n_splits = 1L) {

  if (!missing(subset)) {
    rows <- rlang::eval_tidy(rlang::enquo(subset), data)
  } else {
    rows <- 1:nrow(data)
  }
  all_rownames <- rownames(data)
  object$formula <- update.formula(
    object$formula,
    rlang::expr(!!formula[[2L]] ~ . + !!formula[[3L]])
  )
  
  data <- model.frame(object$formula, data = data[rows, , drop = FALSE])
  
  N <- nrow(data)
  J <- length(object$model_spec)


  
  A <- get_treat_df(object, data)
  object$has_ipw <- object$type %in% c("ipw", "aipw", "did_aipw")
  object$has_outreg <- object$type %in% c("reg_impute", "aipw", "did_aipw", "telescope_match")
  object$has_match <- object$type %in% c("telescope_match")
  object$has_blip <- object$type %in% c("sequential_g")

  check_cde_estimator(object)
  
  if (missing(n_folds)) {
    if (crossfit) {
      n_folds <- 5L
    } else {
      n_folds <- 1L
    }
  }

  object$outcome <- formula[[2L]]

  

    
  out <- list()
  out$type <- object$type
  out$formula <- formula
  out$crossfit <- crossfit
  out$n_folds <- n_folds
  out$model_spec <- object$model_spec
  out$N <- N
  out$included <- all_rownames %in% rownames(data)

  eval_vals <- get_eval_vals(object, data)

  ## match outside the crossfit loop
  if (object$has_match) {
    if (crossfit) {
      rlang::warn("cannot use crossfitting with matching; setting `crossfit = FALSE`")
      crossfit <- FALSE
    }
    out$match_out <- t_match(object, data)
  }


  ## do initial fit without cross-fitting
  out$model_fits <- make_model_fits(A, eval_vals, object$model_spec)
  fit_rows <- 1L:N
  pred_rows <- 1L:N
  res <- fit_fold(object, data, fit_rows, pred_rows, out)
  out$model_fits <- res 
  out <- estimate_cde(object, formula, data, out)

  if (crossfit) {
    
    out$n_folds <- n_folds
    out$n_splits <- n_splits
    split_ests <- matrix(NA, nrow = n_splits, ncol = nrow(out$estimates))
    split_ses <- matrix(NA, nrow = n_splits, ncol = nrow(out$estimates))
    for (s in seq_len(n_splits)) {
      fold_size <- N / n_folds
      f_split <- ceiling(seq_len(N) / fold_size)
      folds <- split(sample(seq_len(N)), f_split)
      
      out$model_fits <- make_model_fits(A, eval_vals, object$model_spec)
      for (k in seq_len(n_folds)) {
        fit_rows <- unlist(folds[-k])
        pred_rows <- unlist(folds[k])
        res <- fit_fold(object, data, fit_rows, pred_rows, out)
        out$model_fits <- res
    }
      out <- estimate_cde(object, formula, data, out)
      split_ests[s, ] <- out$estimates[, "estimate"]
      split_ses[s, ] <- out$estimates[, "std.error"]
    }
  }
  med_ests <- apply(split_ests, 2, median)
  split_vars <- apply(split_ses ^ 2, 2, median) + sweep(split_ests, 2, med_ests) ^ 2
  med_ses <- sqrt(apply(split_vars, 2, median))

  out$estimates[, "estimate"] <- med_ests
  out$estimates[, "std.error"] <- med_ses

  class(out) <- class(object)
  class(out)[2L] <- "cde_estimate"
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
  out$treat_names <- tr_names
  A <- get_treat_df(object, data)
  N <- nrow(A)
  num_treat <- length(tr_names)
  out$estimates <- empty_est_tab()


  eval_vals <- get_eval_vals(object, data)
  paths <- interaction(A, sep = "_")
  path_levs <- levels(paths)
  path_splits <- split(1L:N, paths)
  
  if (class(object) %has% "ipw") {
    for (e in seq_along(eff_vars)) {
      j <- eff_pos[e]
      j_levs <- unique(A[, j])
      out$estimates <- rbind(
        out$estimates,
        compute_ipw(j, j_levs, y, paths, out, object$args, tr_names[j], eval_vals)
      )
    }
  }

  if (class(object) %has% "reg_impute") {
    for (e in seq_along(eff_vars)) {
      j <- eff_pos[e]
      j_levs <- eval_vals[[j]]
      if (length(j_levs) > 1L) { 
        out$estimates <- rbind(
          out$estimates,
          compute_reg_impute(j, j_levs, y, paths, out$model_fits[[e]]$outreg_pred, eff_vars[e])
        )
      }
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

  if (class(object) %has% "did_aipw") {
    j <- eff_pos[1L]
    j_levs <- unique(A[, 1L])
    m0_name <- rlang::as_label(object$args$base_mediator)
    m0 <- data[[m0_name]]
    out$estimates <- rbind(
      out$estimates,
      compute_did_aipw(j, j_levs, y, paths, out, object$args, eff_vars[1L], m0)
    )
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

#' @export
print.cde_estimate <- function(x, ...) {
  cat("\n", x$type, "CDE Estimator\n")
  cat("---------------------------\n")
  cat("Causal variables:", paste0(x$treat_names, collapse = ", "), "\n")
  
  cat("\nEstimated Effects:\n")
  print(x$estimates[, c("term", "active", "control", "estimate")])
  cat("\n")
  invisible()
}

#' @export
tidy.cde_estimate <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
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
summary.cde_estimate <- function(object, ...) {
  
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
  out$model_spec <- object$model_spec
  class(out) <- "summary.cde_estimate"
  out$coefficients <- est
  return(out)
}

#' @export
print.summary.cde_estimate <-  function(x,
                                         digits = max(3L, getOption("digits") - 3L),
                                         signif.stars = getOption("show.signif.stars"),
                                         ...) {
  cat("\n", x$type, "CDE Estimator\n")
  cat("---------------------------\n")
  tr_names <- unlist(
    lapply(x$model_spec, function(x) as.character(x$treat[[2L]]))
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
  cat("Cross-fit:", x$crossfit, "\n")
  if (x$crossfit) cat("Number of folds:", x$n_folds, "\n")
  
  cat("\nEstimated Effects:\n")
  print(
    x$coefficients,
    digits = digits, quote = FALSE, right = TRUE
  )
  cat("\n")

}
