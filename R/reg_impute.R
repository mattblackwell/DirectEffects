#' @export
cde_reg_impute <- function(...) {
  args <- rlang::enquos(...)

  new_cde_estimator(
    "reg_impute",
    args = args,    
    formula = NULL,
    model_spec = NULL
  )
}


get_reg_preds <- function(x, path) {
  num_treat <- length(x$model_fits)
  N <- nrow(x$model_fits[[1L]]$outreg_pred)
  levs <- unlist(strsplit(path, "_"))
  out <- matrix(NA, nrow = N, ncol = num_treat)
  for (j in seq_len(num_treat)) {
    out[, j] <- x$model_fits[[j]]$outreg_pred[, path]
  }
  out
}


compute_reg_impute <- function(j, levs, y, treat, mu_hat, term_name) {
  levs <- sort(levs)
  N <- length(treat)
  paths <- colnames(mu_hat)
  sp <- strsplit(paths, "_")
  templates <- unique(replace_each(sp, j, NA))
  est_tab <- empty_est_tab()
  for (k in seq_along(templates)) {
    base <- templates[[k]]
    base[j] <- levs[1L]
    base <- paste0(base, collapse = "_")
    ctr <- as.numeric(treat == base)
    N_c <- sum(ctr)
    for (p in seq_along(levs[-1L])) {
      plus <- templates[[k]]
      plus[j] <- levs[-1L][p]
      plus <- paste0(plus, collapse = "_")
      trt <- as.numeric(treat == plus)
      N_t <- sum(trt)
      N_b <- N_t + N_c
      ## check here in case eval_vals are not in the observed set. 
      if (N_b > 0 & N_c > 0) {
        psi <- mu_hat[, plus] + (N / N_t) * trt * (y - mu_hat[, plus])
        psi <- psi - mu_hat[, base] - (N / N_c) * ctr * (y - mu_hat[, base])        
      } else {
        ## TODO: check if this is valid for SEs/crossfitting.
        psi <- mu_hat[, plus] - mu_hat[, base]
      }
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
