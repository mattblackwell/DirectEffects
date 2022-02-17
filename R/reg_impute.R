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
  num_treat <- length(x$outreg_pred)
  levs <- unlist(strsplit(path, "_"))
  out <- matrix(NA, nrow = nrow(x$outreg_pred[[1L]]), ncol = num_treat)
  for (j in seq_len(num_treat)) {
    out[, j] <- x$outreg_pred[[j]][, path]
  }
  out
}


compute_reg_impute <- function(j, levs, y, treat, mu_hat, term_name) {
  levs <- sort(levs)
  N <- length(treat)
  paths <- levels(treat)
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
      psi <- mu_hat[, plus] + (N / N_t) * trt * (y - mu_hat[, plus])
      psi <- psi - mu_hat[, base] - (N / N_c) * ctr * (y - mu_hat[, base])
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
