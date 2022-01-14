#' @export
cde_reg_impute <- function(...) {
  args <- rlang::enquos(...)

  new_cde_estimator(
    "reg_impute",
    args = args,
    model_spec = NULL
  )
}



compute_reg_impute <- function(j, levs, y, treat, mu_hat, term_name, N) {
  levs <- sort(levs)
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
      psi <- mu_hat[, plus] + (N_b / N_t) * trt * (y - mu_hat[, plus])
      psi <- psi - mu_hat[, base] - (N_b / N_c) * ctr * (y - mu_hat[, base])
      est <- sum((trt + ctr) * psi) / N_b
      est_var <- sum((trt + ctr) * (psi - est)^ 2) / N_b ^ 2
      this_est <- data.frame(
        term = term_name,
        block_num = j,
        active = format_path(plus),
        control = format_path(base),
        estimate = est,
        std_err = sqrt(est_var)
      )
      est_tab <- rbind(est_tab, this_est)      
    }
  }
  rownames(est_tab) <- NULL
  est_tab

}
