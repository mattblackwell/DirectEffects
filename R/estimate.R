#' @export
estimate <- function(object, formula, data, subset, crossfit = TRUE, n_folds, fold_seed = NULL) {

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
        compute_ipw(j, j_levs, y, paths, out, object$args, tr_names[j])
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
