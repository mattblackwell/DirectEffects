#' Initialize an AIPW DID-CDE estimator 
#'
#' Initializes the specification of a difference-in-differences
#' estimator for the CDE based on an augmented inverse probability
#' weighting.
#' 
#' @param base_mediator The (unquoted) name of the variable that
#' measures the mediator at baseline.
#' @inheritParams cde_aipw
#' @param on_treated If `FALSE` (the defafult), the effects are
#' average effects conditional on the levels of the baseline mediator.
#' If `TRUE`, the effects are conditional on the treated path. For
#' difference in identficiation, see Details below.
#'
#' @details
#' This function, unlike other CDE estimators in the package, only
#' returns the estimated effects of the first treatment variable.
#' These effects are conditional on the baseline value of the mediator
#' (`base_mediator`) when `on_treated` is `TRUE`. A marginalized CDE
#' estimand is also estimated. When `on_treated` is `FALSE`, these
#' estimates are conditional on the entire "treated" history.
#' Identification  requirements are slightly different between these
#' two cases. When `on_treated` is `FALSE`, the confounders for the
#' mediator cannot be affected by treatment. See Blackwell et al
#' (2022) for more information. 
#' 
#' @export
#' @md
cde_did_aipw <- function(
                        base_mediator,
                        trim = c(0.01, 0.99),
                        aipw_blip = TRUE,
                        on_treated = FALSE
                        ) {
  args <- list(
    base_mediator = rlang::enquo(base_mediator),
    trim = trim,
    aipw_blip = aipw_blip,
    on_treated = on_treated
  )
  new_cde_estimator(
    "did_aipw",
    args = args,    
    formula = NULL,
    model_spec = NULL
  )
}


compute_did_aipw <- function(j, j_levs, y, treat, out, args, term_name, m0) {
  num_treat <- length(out$model_fits)
  N <- length(treat)
  j_levs <- sort(j_levs)
  paths <- colnames(out$model_fits[[j]]$outreg_pred)
  
  sp <- strsplit(paths, "_")
  templates <- unique(replace_each(sp, j, NA))

  est_tab <- empty_est_tab()
  pi_0 <- rep(NA, times = length(templates))
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

    M_0 <- as.numeric(m0 == templates[[k]][2])
    pi_0[k] <- mean(M_0)
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
        w_trt <- winsorize_matrix(w_trt, args$trim)
        w_ctr <- winsorize_matrix(w_ctr, args$trim)
        p_trt <- winsorize_matrix(p_trt, args$trim)
        p_ctr <- winsorize_matrix(p_ctr, args$trim)
      }
      w_trt <- cbind(1, w_trt)
      w_ctr <- cbind(1, w_ctr)
      if (args$on_treated) {
        ## assumes only 2 treatments
        a_trt <- M_0 * A_trt[, 2] * A_trt[, 3]
        a_ctr <- M_0 * A_ctr[, 2] * A_ctr[, 3]
        ipw <-  1 / mean(M_0 * A_trt[, 2] * A_trt[, 3])
        psi_trt <- a_trt * (eps_trt[, 3] + r_trt[, 2] - r_ctr[, 2])
        
        p_ratio <- p_trt[, 2] / p_ctr[, 2]

        psi_ctr <- a_ctr * (p_ratio * eps_ctr[, 3])
        psi <- ipw * (psi_trt - psi_ctr)
        est <- mean(psi)
        est_var <- mean((psi - ipw * a_trt * est) ^ 2) / N
      } else {
        psi_trt <- rowSums(A_trt * eps_trt / w_trt)
        psi_ctr <- rowSums(A_ctr * eps_ctr / w_ctr)
        ipw <- M_0 / mean(M_0)
        psi <- ipw * (psi_trt - psi_ctr)
        est <- mean(psi)
        est_var <- mean((psi - ipw * est)^ 2) / N
      }

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

  marg_est <- sum(est_tab$estimate * pi_0)
  marg_var <- sum(est_tab$std.err ^ 2 * pi_0 ^ 2)
  plus <- templates[[k]]
  plus[j] <- j_levs[-1L][p]
  
  plus <- paste0(plus, collapse = "_")
  this_est <- data.frame(
    term = term_name,
    active = paste0(j_levs[2L], "_*"),
    control = paste0(j_levs[1L], "_*"),
    estimate = marg_est,
    std.error = sqrt(marg_var),
    DF = N
  )
  est_tab <- rbind(est_tab, this_est)
  rownames(est_tab) <- NULL
  est_tab
}
