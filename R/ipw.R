#' @export
cde_ipw <- function(hajek = TRUE, trim = c(0.01, 0.99)) {
  args <- list(
    hajek = hajek,
    trim = trim
  )

  new_cde_estimator(
    "ipw",
    args = args,
    formula = NULL,
    model_spec = NULL
  )
}


get_ipw_preds <- function(x, levs) {  
  N <- nrow(x$model_fits[[1L]]$treat_pred)
  J <- length(x$model_fits)
  out <- matrix(NA, nrow = N, ncol = J,
                dimnames = list(NULL, seq_len(J)))
  for (j in seq_len(J)) {    
    label <- subset_history_string(levs, 1L:j)
    colnames(out)[j] <- label
    out[, j] <- x$model_fits[[j]]$treat_pred[, label]
  }
  out
}



compute_ipw <- function(j, j_levs, y, treat, out, args, term_name, eval_vals) {
  num_treat <- length(out$model_fits)
  N <- length(treat)
  j_levs <- sort(j_levs)
  paths <- create_history_strings(eval_vals, 1L:num_treat)
  sp <- strsplit(paths, "_")
  templates <- unique(replace_each(sp, j, NA))
  est_tab <- empty_est_tab()
  for (k in seq_along(templates)) {
    base <- templates[[k]]
    base[j] <- j_levs[1L]
    base <- paste0(base, collapse = "_")
    ctr <- as.numeric(treat == base)
    N_c <- sum(ctr)
    p_ctr <- get_ipw_preds(out, base)
    w_ctr <- apply(p_ctr, 1, prod)
    for (p in seq_along(j_levs[-1L])) {
      plus <- templates[[k]]
      plus[j] <- j_levs[-1L][p]
      plus <- paste0(plus, collapse = "_")
      trt <- as.numeric(treat == plus)
      N_t <- sum(trt)
      p_trt <- get_ipw_preds(out, plus)
      w_trt <- apply(p_trt, 1, prod)
      

      if (length(args$trim)) {
        w_trt <- winsorize(w_trt, args$trim)
        w_ctr <- winsorize(w_ctr, args$trim)
      }
      psi_trt <- trt * y / w_trt
      psi_ctr <- ctr * y / w_ctr
      
      psi <- psi_trt - psi_ctr
      est <- mean(psi)
      est_var <- mean((psi - est)^ 2) / N
      this_est <- data.frame(
        term = term_name,
        active = plus,
        control = base,
        estimate = est,
        std.error = sqrt(est_var),
        DF = N_c + N_t
      )
      est_tab <- rbind(est_tab, this_est)      
    }
  }
  rownames(est_tab) <- NULL
  est_tab
}


winsorize <- function(x, trim) {
  qs <- quantile(x, trim)
  xt <- x
  xt[xt <= qs[1L]] <- qs[1L]
  xt[xt >= qs[2L]] <- qs[2L]
  xt
}

winsorize_matrix <- function(x, trim) {
  for (j in seq_len(ncol(x))) {
    x[, j] <- winsorize(x[, j], trim)
  }
  x
}
