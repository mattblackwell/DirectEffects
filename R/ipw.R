#' @export
cde_ipw <- function(hajek = TRUE, trim = c(0.01, 0.99)) {
  args <- list(
    hajek = hajek,
    trim = trim
  )

  new_cde_estimator(
    "ipw",
    args = args,
    model_spec = NULL
  )
}


get_ipw_preds <- function(x, levs) {  
  levs <- unlist(strsplit(levs, "_"))
  if (length(x$ipw_pred) != length(levs)) {
    rlang::abort("levs must match number of causal blocks with ipw predictions")
  }
  out <- matrix(NA, nrow = nrow(x$ipw_pred[[1L]]), ncol = length(levs),
                dimnames = list(NULL, 1L:length(levs)))
  for (j in seq_along(levs)) {    
    label <- paste0(levs[1L:j], collapse = "_")
    colnames(out)[j] <- label
    out[, j] <- x$ipw_pred[[j]][, label]
  }
  out
}



compute_ipw_contrasts <- function(j, levs, psi, psi_sq, term_name, N) {
  levs <- sort(levs)
  paths <- names(psi)
  sp <- strsplit(paths, "_")
  templates <- unique(replace_each(sp, j, NA))
  est_tab <- empty_est_tab()
  for (k in seq_along(templates)) {
    base <- templates[[k]]
    base[j] <- levs[1L]
    base <- paste0(base, collapse = "_")
    
    for (p in seq_along(levs[-1L])) {
      plus <- templates[[k]]
      plus[j] <- levs[-1L][p]
      plus <- paste0(plus, collapse = "_")
      
      est <- psi[plus] - psi[base]
      est_var <- (psi_sq[plus] + psi_sq[base] - est ^ 2) / N
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

trim_weights <- function(x, trim) {
  ## todo: add checks for trim here or upon user input
  for (j in seq_len(ncol(x))) {
    qs <- quantile(x, trim)
    xt <- x[, j]
    xt[xt <= qs[1]] <- qs[1]
    xt[xt >= qs[2]] <- qs[2]
    x[, j] <- xt
  }
  x
}
