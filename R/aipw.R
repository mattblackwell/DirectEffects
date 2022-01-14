#' @export
cde_aipw <- function(trim = c(0.01, 0.99)) {
  args <- list(
    trim = trim
  )

  new_cde_estimator(
    "aipw",
    args = args,
    model_spec = NULL
  )
}



compute_aipw <- function(j, j_levs, y, treat, out, args, term_name, N) {
  num_treat <- length(out$outreg_pred)
  N <- nrow(out$outreg_pred[[1L]])
  j_levs <- sort(j_levs)
  paths <- colnames(out$outreg_pred[[j]])
  
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
    w_ctr <- t(apply(p_ctr, 1, cumprod))
    w_ctr <- w_ctr[, j:num_treat, drop = FALSE]
    r_ctr <- get_reg_preds(out, base)
    r_ctr <- r_ctr[, j:num_treat, drop = FALSE]
    A_ctr <- get_path_inds(treat, base)
    A_ctr <- cbind(1, A_ctr[, j:num_treat, drop = FALSE])
    eps_ctr <- cbind(r_ctr, y) - cbind(0, r_ctr)

    for (p in seq_along(j_levs[-1L])) {
      plus <- templates[[k]]
      plus[j] <- j_levs[-1L][p]
      plus <- paste0(plus, collapse = "_")
      trt <- as.numeric(treat == plus)
      N_t <- sum(trt)
      p_trt <- get_ipw_preds(out, plus)
      w_trt <- t(apply(p_trt, 1, cumprod))
      w_trt <- w_trt[, j:num_treat, drop = FALSE]
      r_trt <- get_reg_preds(out, plus)
      r_trt <- r_trt[, j:num_treat, drop = FALSE]
      A_trt <- get_path_inds(treat, plus)
      A_trt <- cbind(1, A_trt[, j:num_treat, drop = FALSE])
      eps_trt <- cbind(r_trt, y) - cbind(0, r_trt)
      

      N_b <- N_t + N_c
      if (length(args$trim)) {
        w_trt <- trim_weights(w_trt, args$trim)
        w_ctr <- trim_weights(w_ctr, args$trim)
      }
      w_trt <- cbind(1, w_trt)
      w_ctr <- cbind(1, w_ctr)
      psi_trt <- rowSums(A_trt * eps_trt / w_trt)
      psi_ctr <- rowSums(A_ctr * eps_ctr / w_ctr)
      
      psi <- psi_trt - psi_ctr
      est <- mean(psi)
      est_var <- mean((psi - est)^ 2) / N
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
