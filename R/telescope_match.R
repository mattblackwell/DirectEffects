#' Perform telescope matching to estimate the controlled
#' direct effect of a binary treatment net the effect of a binary mediator.
#'
#' @param formula A formula object that specifies the outcome,
#'   baseline covariates, treatment, intermediate covariate, and
#'   mediator to be used in the matching. Each of the last four
#'   variable groups should be separated by \code{|}. See below for
#'   more details.
#' @param data A dataframe containing columns referenced by
#'   \code{formula}.
#' @param caliper A scalar denoting the caliper to be used in matching
#'   in the treatment stage (calipers cannot be used for matching on
#'   the mediator). Observations outside of the caliper are dropped.
#'   Calipers are specified in standard deviations of the covariates.
#'   NULL by default (no caliper).
#' @param L Number of matches to use for each unit. Must be a numeric
#'   vector of eitehr length 1 or 2. If length 1, L sets the number of
#'   matches used in both the first stage (matching on mediator) and
#'   in the second stage (matching on treatment). If length 2, the
#'   first element sets the number of matches used in the first stage
#'   (matching on mediator) and the second element sets the number of
#'   matches used in the second stage (matching on treatment) Default
#'   is 5.
#' @param boot logical indicating whether to conduct inference using
#'   the weighted bootstrap for matching estimators extended from Otsu
#'   and Rai (2017) (\code{TRUE}) or the asymptotic variance estimator
#'   from Blackwell and Strezhnev (2019) (\code{FALSE}). Defaults to
#'   \code{FALSE}.
#' @param nBoot If \code{boot} is \code{TRUE}, number of bootstrap
#'   iterations to use. Default is \code{5000}.
#' @param ci percent level of confidence interval to return. If
#'   \code{boot} is \code{FALSE}, returns symmetric asymptotic
#'   interval centered on the estimated treatment effect. If
#'   \code{boot} is \code{TRUE} returns the \code{(100 - ci)/2} and
#'   \code{100 - (100 - ci)/2} percentiles of the bootstrap
#'   distribution. Must be in the interval \code{(0, 100)}. Defaults
#'   to 95.
#' @param verbose logical indicating whether to display progress
#'   information. Default is \code{TRUE}.
#' @param subset A vector of logicals indicating which rows of
#'   \code{data} to keep.
#' @param contrasts a list to be passed to the \code{contrasts.arg}
#' argument of \code{model.matrix()} when generating the data matrix.
#' @param separate_bc logical indicating whether or not bias
#'   correction regressions should be run separately within levels of
#'   the treatment and mediator. Defaults to \code{TRUE}. If
#'   \code{TRUE}, any interactions between treatment/mediator and
#'   covariates in the specification should be omitted.
#' @inheritParams stats::lm
#'
#' @details The \code{telescope_match} function implements the
#'   two-stage "telescope matching" procedure developed by Blackwell
#'   and Strezhnev (2020).
#'
#'  The procedure first estimates a demediated outcome using a
#'  combination of matching and a regression bias-correction based on
#'  the specification in \code{formula}, which should be specified as
#'  \code{Y ~ X | A | Z | M}, where \code{Y} is the outcome, \code{X}
#'  is a formula of baseline covariates, \code{A} is a single variable
#'  name indicating the binary treatment, \code{Z} is a formula of
#'  intermediate covariates, and \code{M} is a single variable name
#'  indicating the mediator.
#'
#' Under the default \code{separate_bc == TRUE}, the first stage will
#' match for \code{M} on \code{X} and \code{Z} within levels of
#' \code{A}. The bias correction regressions will be linear with a
#' specification of \code{Y ~ X + Z} within levels of \code{A} and
#' \code{M}. In the second stage, we match for \code{A} within levels
#' of \code{X} and run bias correction regressions of \code{Ytilde ~
#' X} within levels of \code{A}, where \code{Ytilde} is the imputed
#' value of \code{Y(A, 0)} from the first stage. This second step
#' estimates the controlled direct effect of treatment when the
#' mediator is 0.
#'
#' When \code{separate_bc} is \code{FALSE}, the bias correction
#' regressions are not broken out by the treatment/mediator and
#' \code{A} and \code{M} are simply included as separate terms in the
#' linear regression. In this setting, interactions between the
#' treatment/mediator and covariates can be added on a selective basis
#' to the \code{X} and \code{Z} specifications.
#'
#'  Matching is performed using the \code{Match()} routine from the
#' \code{Matching} package. By default, matching is L-to-1 nearest
#' neighbor with replacement using Mahalanobis distance.
#'
#'
#' See the references below for more details.
#'
#' @return Returns an object of \code{class} \code{tmatch}. Contains
#' the following components
#' \itemize{
#' \item estimate: Estimated ACDE fixing M=0
#' \item std.err: Estimated asymptotic standard error. \code{NULL} if
#' \code{boot} is \code{TRUE}
#' \item boot.dist: Bootstrap distribution of \code{estimate}.
#' \code{NULL} if  \code{boot} is \code{FALSE}
#' \item conf.low: Lower bound of \code{ci} confidence interval for
#' the estimate
#' \item conf.high: Upper bound of \code{ci} confidence interval for
#' the estimate
#' \item ci.level: Level of the confidence interval
#' \item outcome: Name of outcome variable
#' \item treatment: Name of treatment variable
#' \item mediator: Name of mediator variable
#' \item pre.treatment: Vector of names of pre-treatment confounders
#' (appear in both stage 1  and 2)
#' \item post.treatment: Vector of names of post-treatment confounders
#' (appear only in stage 1)
#' \item s1.formula: Stage 1 bias-correction regression formula
#' (pre-/post-treatment covariates)
#' \item s2.formula: Stage 2 bias-correction regression formula
#' (pre-treatment covariates)
#' \item outcome.vec: Vector of outcomes used in estimation
#' \item treatment.vec: Vector of treatment indicators used in
#' estimation
#' \item mediator.vec: Vector of mediator indicators used in
#' estimation
#' \item L_m: Number of matches found for each unit in the first stage
#' mediator matching  procedure
#' \item L_a: Number of matches found for each unit in the second
#' stage mediator matching  procedure
#' \item KLm: Number of times unit is used as a match in the first
#' stage mediator matching  procedure
#' \item KLa: Number of times unit is used as a match in the second
#' stage treatment matching  procedure
#' \item N: Number of observations
#' \item N_summary: Number of observations in each treatment/mediator
#' combination.
#' \item caliper: Caliper (if any) used in the treatment stage
#' matching to drop distant  observations.
#' }

#' @references Blackwell, Matthew, and Strezhnev, Anton (2020)
#' "Telescope Matching: Reducing Model Dependence
#' in the Estimation of Direct Effects." Working Paper.
#'
#' @examples
#' data(jobcorps)
#'
#' ## Split male/female
#' jobcorps_female <- subset(jobcorps, female == 1)
#'
#' ## Telescope matching formula - First stage (X and Z)
#' tm_form <- exhealth30 ~  schobef + trainyrbef + jobeverbef  |
#' treat | emplq4 + emplq4full | work2year2q
#'
#'
#' ### Estimate ACDE for women holding employment at 0
#' tm_out <-  telescope_match(
#'   tm_form,
#'   data = jobcorps_female,
#'   L = 3,
#'   boot = FALSE,
#'   verbose = TRUE
#' )
#'
#' @export
#' @importFrom Matching Match
#' @importFrom stats as.formula model.frame predict rbinom var
#'
telescope_match <- function(formula, data, caliper = NULL, L = 5,
                            boot = FALSE, nBoot = 5000, ci = 95,
                            verbose = TRUE, subset, contrasts = NULL,
                            separate_bc = TRUE, ...) {

  ########################
  ### Quick pre-processing
  ########################


  mf <- match.call(expand.dots = FALSE)

  m <- match(
    x = c("formula", "data", "subset", "na.action"),
    table = names(mf),
    nomatch = 0L
  )
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE

  ## must be valid formula
  formula <- Formula::as.Formula(formula)
  T <- length(formula)[2] / 2
  if (T < 2 | (T %% 1) != 0) {
    stop(
      "Right-hand side of `formula` must have an even number of groups.",
      "\nℹ Groups should be divded by `|`", call. = FALSE
    )
  }

  ### If length of L is 1, L_a and L_m are the same
  if (!is.numeric(L) | !is.vector(L) | !(length(L)  %in% c(1, T))) {
    stop("Error: `L` must be a numeric vector of length 1 or T", call. = FALSE)
  }

  if (length(L) == 1) {
    L <- rep(L, times = T)
  }

  stopifnot(length(formula)[1] == 1L)
  if (inherits(try(terms(formula), silent = TRUE), "try-error")) {
    stop("cannot use dot '.' in formulas")
  }

  ## add to mf call
  mf$formula <- formula
  #  finally evaluate model.frame, create data matrix
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # terms object

  all_data <- model.matrix(mt, mf, contrasts)
  Y <- model.response(mf, "numeric")

  N <- length(Y)

  X <- list()
  drops <- list()
  A <- matrix(NA, nrow = N, ncol = T)
  anames <- rep("", times = T)

  x_spots <- seq(1, T * 2, by = 2)
  a_spots <- seq(2, T * 2, by = 2)
  for (s in seq_len(T)) {
    ## covariates up to s
    mt_a <- terms(formula, data = data, rhs = a_spots[s])
    if (length(attr(mt_a, "term.labels")) != 1) {
      stop("Error: there must only be single treatment in each period.",
           call. = FALSE)
    }
    anames[s] <- attr(mt_a, "term.labels")
    A[, s] <- all_data[, anames[s]]
    mt_x <- terms(formula, data = data, rhs = x_spots[1:s])
    X[[s]] <- model.matrix(mt_x, mf, contrasts)

    ## treatnent in covariate specifications. need to dropped for
    ## matching
    a_in_x <- which(rownames(attributes(mt_x)$factors) %in% anames[1:s])
    drops[[s]] <- which(attributes(mt_x)$factors[a_in_x, ] == 1)
    if (colnames(X[[s]])[1] == "(Intercept)") drops[[s]] <- c(1, drops[[s]] + 1)

    if (separate_bc & length(a_in_x)) {
      stop(
        "Error: covariate specification can't contain treatment terms",
        "when `separate_bc = TRUE`.",
        call. = FALSE
      )
    }
  }


  ####################
  #### Sanity checks
  ####################

  m_out <- list()
  r_out <- list()
  K <- list()
  K_paths <- list()

  ## matching
  for (s in T:1) {
    A_hist <- A[, 1:s, drop = FALSE]
    m_out[[s]] <- match_at_time(A_hist, X[[s]], L[s], drops[[s]], caliper)

    ## indirect matches are a pain with T > 2 bc an indirect match can
    ## come from t=1 to t=2, t=1 to t=3, or t=1 to t=2 to t=3. using
    ## `combn` here to set up those possibilities and
    ## `find_indirect_matches` "crawls" the path get the indirect
    ## matches for each unit.

    ## NB: Ks are divided by L values.
    path_lengths <- 1:(T - s + 1)
    paths <- lapply(
      path_lengths,
      function(x) combn(s:T, m = x, simplify = FALSE)
    )
    paths <- unlist(paths, recursive = FALSE)
    paths <- paths[unlist(lapply(paths, function(x) x[1] == s))]
    K[[s]] <- lapply(
      paths,
      function(x) find_indirect_matches(m_out, path = x, N) / prod(L[x])
    )
    K_paths[[s]] <- paths
    names(K[[s]]) <- unlist(lapply(paths, paste0, collapse = "_"))
    names(K[[s]]) <- paste0("K_", names(K[[s]]))
  }

  K <- as.data.frame(unlist(K, recursive = FALSE))

  ## bias correction
  A_levs <- expand.grid(rep(list(c(0, 1)), times = T - 1))
  tau_i <- matrix(NA, nrow = N, ncol = nrow(A_levs))
  tau_raw <- rep(0, times = nrow(A_levs))
  tau_se <- rep(0, times = nrow(A_levs))
  mu_hat <- matrix(NA, nrow = N, ncol = T)
  Yt <- matrix(NA, nrow = N, ncol = T * nrow(A_levs))
  for (j in seq_len(nrow(A_levs))) {
    A_j <- A_levs[j, ]
    Ytilde <- Y
    for (s in T:1) {
      A_hist <- A[, 1:s, drop = FALSE]
      r_out[[s]] <- regress_at_time(Ytilde, A_hist, X[[s]], separate_bc)
      yhat_r <- A_j[s - 1] * r_out[[s]]$yhat_r_1 +
        (1 - A_j[s - 1]) * r_out[[s]]$yhat_r_0

      yhat_mr_0 <- lapply(
        m_out[[s]]$matches,
        function(x) mean(r_out[[s]]$yhat_r_0[x])
      )
      yhat_mr_1 <- lapply(
        m_out[[s]]$matches,
        function(x) mean(r_out[[s]]$yhat_r_1[x])
      )
      yhat_mr <- A_j[s - 1] * unlist(yhat_mr_1) +
        (1 - A_j[s - 1]) * unlist(yhat_mr_0)

      yhat_m <- unlist(lapply(m_out[[s]]$matches, function(x) mean(Ytilde[x])))

      not_a_j <- A[, s] != A_j[s - 1]

      Ytilde[not_a_j] <- yhat_m[not_a_j] + (yhat_r[not_a_j] - yhat_mr[not_a_j])
      Yt[, j + (T - s)] <- Ytilde

      if (s > 1) {
        mu_hat[, s] <- yhat_r
      }
    }
    mu_hat[, 1] <- A[, 1] * r_out[[1]]$yhat_r_1 +
      (1 - A[, 1]) * r_out[[1]]$yhat_r_0

    cdes <- calculcate_cdes(Y, A, K, mu_hat, r_out, A_j)
    tau_i[, j] <- cdes$tau_i
    tau_raw[j] <- cdes$tau_raw
    tau_se[j] <- cdes$tau_se
  }
  tau <- colMeans((2 * A[, 1] - 1) * tau_i)


  ### Return output
  out <- list(formula = formula, m_out = m_out, K = K, r_out = r_out, tau = tau,
                 tau_raw = tau_raw, tau_se = tau_se, Ytilde = Yt)

  class(out) <- "tmatch"
  return(out)

}


## given the matching output and a set of periods `path = (k, ..., j)`
## find the number of times units are used indirectly as a match in
## period `k` through the given path. E.g., if `path = c(1, 2 ,3)`,
## then this returns the number of times each unit is matched to a
## unit j in period 3, j is matched to h in period 2, and h is matched
## to g in period 1, where j, h, and g are all different.
find_indirect_matches <- function(m_out, path, N) {
  K <- rep(0, times = N)
  path <- sort(path, decreasing = TRUE)
  d <- m_out[[path[1]]]$donors
  path <- path[-1]
  for (j in path) {
    d <- lapply(d, function(x) unlist(m_out[[j]]$donors[as.character(x)]))
  }
  n_ind_matches <- lapply(d, length)

  pos <- as.numeric(names(n_ind_matches))
  K[pos] <- unlist(n_ind_matches)
  return(K)
}

## working assumption: the A matrix is n \times k where the first
## column is t=1 and the last column is t=T-k+1.
match_at_time <- function(A, X, L, drops, caliper) {
  k <- ncol(A)
  Ak <- A[, k]
  if (!isTRUE(all.equal(unique(Ak)[order(unique(Ak))], c(0, 1)))) {
    stop(
      "Error: treatments must only take on values 0 or 1",
      "\n ✖ treatment ", k, " non-binary values."
    )
  }
  X_m <- cbind(X[, -drops], A[, -k])

  exacts <- c(
    rep(FALSE, times = ncol(X[, -drops])),
    rep(TRUE, times = k - 1)
  )

  tm_att <- Matching::Match(
    Tr = Ak, X = X_m,
    estimand = "ATT",
    exact = exacts,
    caliper = caliper,
    M = L, ties = FALSE,
    Weight = 2
  )
  tm_atc <- Matching::Match(
    Tr = Ak, X = X_m,
    estimand = "ATC",
    exact = exacts,
    caliper = caliper,
    M = L, ties = FALSE,
    Weight = 2
  )

  ## donation lists: where each unit was used as a match
  donors_1 <- split(tm_att$index.treated, tm_att$index.control)
  donors_0 <- split(tm_atc$index.control, tm_atc$index.treated)
  donors <- c(donors_1, donors_0)
  donors <- donors[as.character(sort(as.numeric(names(donors))))]

  ## matched sets: matched for each unit
  matches_1 <- split(tm_att$index.control, tm_att$index.treated)
  matches_0 <- split(tm_atc$index.treated, tm_atc$index.control)
  matches <- c(matches_1, matches_0)
  matches <- matches[as.character(sort(as.numeric(names(matches))))]


  out <- list(matches = matches, donors = donors)

  return(out)

}


## calculates fitted values for all units under treatment at time `k`
regress_at_time <- function(y, A, X, separate_bc) {
  N <- nrow(A)
  k <- ncol(A)
  Ak <- A[, k]
  A_past <- A[, -k, drop = FALSE]
  if (ncol(A_past) == 0) A_past <- rep(0, times = N)

  if (separate_bc) {
    A_fact <- do.call(interaction, list(as.data.frame(A_past), sep = "_"))
    A_split <- split(1:N, A_fact)
    preds_0 <- lapply(
      X = A_split,
      FUN = strata_reg_predict,
      y = y, x = X, Ak = Ak, Ak_lev = 0
    )
    preds_1 <- lapply(
      X = A_split,
      FUN = strata_reg_predict,
      y = y, x = X, Ak = Ak, Ak_lev = 1
    )

    n_coefs <- length(A_split) * ncol(X) * 2
    yhat_r_0 <- unsplit(preds_0, A_fact)
    yhat_r_1 <- unsplit(preds_1, A_fact)
  } else {
    X_des <- cbind(X, A)
    cf <- lm.fit(y = y, x = X_des)$coefficients
    cf <- cf[!is.na(cf)]
    X_des[, ncol(X_des)] <- 0
    yhat_r_0 <- X_des[, names(cf)] %*% cf
    X_des[, ncol(X_des)] <- 1
    yhat_r_1 <- X_des[, names(cf)] %*% cf
    n_coefs <- ncol(X_des)
  }
  out <- list(yhat_r_0 = yhat_r_0, yhat_r_1 = yhat_r_1, n_coefs = n_coefs)
  return(out)
}


## this function allows us to send strata based on the past treatment
## values (where the rows indicate observations within those strata)
## and calculate predicted values from the regression among those with
## `A_t_lev` current treatment status.
strata_reg_predict <- function(rows, y, x, Ak, Ak_lev) {
  rows_lev <- rows[which(Ak[rows] == Ak_lev)]

  cf <- lm.fit(y = y[rows_lev], x = x[rows_lev, ])$coefficients
  cf <- cf[!is.na(cf)]
  pred <- x[rows, names(cf)] %*% cf

  return(as.vector(pred))
}

calculcate_cdes <- function(Y, A, K, mu_hat, r_out, A_j) {
  N <- length(Y)
  T <- ncol(A)

  ## calculate simple matching parts
  sm_part <- (1 + rowSums(K)) * Y
  A_fut <- A[, -1, drop = FALSE]
  A_fut_str <-  apply(A_fut, 1, paste0, collapse = "")
  i_j <- as.numeric(A_fut_str == paste0(A_j, collapse = ""))

  tau_i <- i_j * sm_part
  tau_raw <- mean((2 * A[, 1] - 1) * tau_i)

  mu_hat_A <- A[, 1] * r_out[[1]]$yhat_r_1 +
    (1 - A[, 1]) * r_out[[1]]$yhat_r_0

  mu_hat_1_A <- A[, 1] * r_out[[1]]$yhat_r_0 +
      (1 - A[, 1]) * r_out[[1]]$yhat_r_1

  K_term <- i_j * rowSums(K) - (1 - i_j)
  tau_i <- tau_i - K_term * mu_hat_A - mu_hat_1_A

  p <- sapply(r_out, function(x) x$n_coefs)
  dfc_1 <- sum(i_j) / (sum(i_j) - p[T])
  epsilon <- Y - mu_hat[, T]
  tau_var <- dfc_1 * mean(i_j * (1 + rowSums(K))^2 * epsilon ^ 2)

  for (s in 2:T) {
    ## get indices for A_{s+1}:T == A_j
    if (s < T) {
      A_fut <- A[, (s + 1):T, drop = FALSE]
      A_fut_str <-  apply(A_fut, 1, paste0, collapse = "")
      i_j <- A_fut_str == paste0(A_j[(s - 1):(T - 1)], collapse = "")
      i_j <- as.numeric(i_j)
    } else {
      i_j <- 1
    }
    A_jt <- as.numeric(A[, s] == A_j[s - 1])

    ## calculate K weight for this period
    fut_patt <- paste0("[", paste0(s:T, collapse = ""), "]")
    past_K <- grep(fut_patt, names(K), invert = TRUE)
    fut_K <- grep(fut_patt, names(K))
    K_term <- (A_jt * rowSums(K[fut_K]) - (1 - A_jt) * (1 + rowSums(K[past_K])))

    eta <- mu_hat[, s] - mu_hat[, s - 1]
    tau_i <- tau_i - as.numeric(i_j) * K_term * eta

    mu_lag_A <- A[, s - 1] * r_out[[s - 1]]$yhat_r_1 +
      (1 - A[, s - 1]) * r_out[[s - 1]]$yhat_r_0
    eta_v <- mu_hat[, s] - mu_lag_A
    dfc_s <- sum(i_j) / (sum(i_j) - sum(p[1:s]))
    tau_var <- tau_var +
      dfc_s * mean(i_j * (1 + rowSums(K[past_K])) ^ 2 * eta_v ^ 2)
  }

  tau <- mean((2 * A[, 1] - 1) * tau_i)

  tau_x <- r_out[[1]]$yhat_r_1 - r_out[[1]]$yhat_r_0
  dfc_t <- N / (N - sum(p) - 1)
  tau_var <- tau_var + dfc_t * mean((tau_x - tau) ^ 2)
  tau_se <- sqrt(tau_var / N)

  return(list(tau = tau, tau_i = tau_i, tau_raw = tau_raw, tau_se = tau_se))
}


boot_tm <- function(obj, R = 100, ci_alpha = 0.05) {

    ## De-mean
    tau.norm <- obj$tau.i - obj$tau

    ## Bootstrap iterations
    W.bern <- sapply(
      1:R,
      function(x) rbinom(N, 1, prob = (sqrt(5) - 1) / (2 * sqrt(5)))
    )
    Wstar <- (((sqrt(5) + 1) / 2) * W.bern +
                ((-sqrt(5) + 1) / 2) * (1 - W.bern)) / N

    ## Apply bootstrap weights to each "observation"
    Tstar <- (t(Wstar) %*% tau.norm)

    ## Add in tau
    Tstar <- Tstar + obj$tau

    ### Get quantiles for the CI
    ci.low <- quantile(Tstar, ci_alpha / 2)
    ci.high <- quantile(Tstar, 1 - ci_alpha / 2)

  return(list(ci.low = ci.low, ci.high = ci.high))

}

#' Summarize telescope match objects
#'
#' @details \code{summary} method for \code{tmatch} objects returned by
#' \code{telescope_match}
#'
#' @param object an object of class \code{tmatch} -- results from a
#' call to \code{telescope_match}
#' @param ... additional arguments affecting the summary produced.
#'
#' @details Returns a summary data frame containing the estimate,
#' standard error and confidence  interval from the `telescope_match`
#' object.
#'
#' @return Returns an object of \code{class} \code{summary.tmatch}.
#' Contains the following components
#' \itemize{
#' \item outcome: Name of outcome variable
#' \item treatment: Name of treatment variable
#' \item mediator: Name of mediator variable
#' \item pre.treatment: Vector of names of pre-treatment confounders
#' (appear in both stage 1 and 2)
#' \item post.treatment: Vector of names of post-treatment confounders
#' (appear only in stage 1)
#' \item sizes: Number of observations in each treatment-mediator combination
#' \item L_m: Number of units matched in the mediator (first) stage
#' \item L_a: Number of units matched in the treatment (second) stage
#' \item estimate: Point estimate of the ACDE
#' \item se.type: Character indicating the type of standard error
#' estimator used in the
#' telescope_matching routine, either 'Asymptotic" or 'Bootstrap'
#' \item std.err: Estimated standard error of the ACD
#' \item ci.level: Confidence interval level
#' \item conf.low: Lower bound of the `ci.level` asymptotically normal
#' confidence interval
#' \item conf.high: Upper bound of the `ci.level` asymptotically
#' normal confidence interval
#' }
#' @export
summary.tmatch <- function(object, ...) {

  summary_obj <- NULL

  ### Names of treatments/outcome
  summary_obj$outcome <- object$outcome
  summary_obj$treatment <- object$treatment
  summary_obj$mediator <- object$mediator

  ### Covariates
  summary_obj$pre.treatment <- object$pre.treatment
  summary_obj$post.treatment <- object$post.treatment

  ### Sample size
  summary_obj$sizes <- object$N_summary

  ### Number of units matched in each stage
  summary_obj$L_m <- object$L_m
  summary_obj$L_a <- object$L_a

  ### Estimates
  if (!is.null(object$std.err)) {
    summary_obj$se.type <- "Asymptotic"
    summary_obj$std.err <- object$std.err
  } else {
    summary_obj$se.type <- "Bootstrap"
    summary_obj$std.err <- sd(object$boot.dist)
  }

  summary_obj$estimate <- object$estimate
  summary_obj$ci.level <- object$ci.level
  summary_obj$conf.low <- object$conf.low
  summary_obj$conf.high <- object$conf.high

  class(summary_obj) <- "summary.tmatch"

  return(summary_obj)

}

print.summary.tmatch <- function(object, digits = max(3, getOption("digits") - 3)) {
  cat("Telescope matching results\n")
  cat("----------------------------\n")
  cat("Variables:\n")
  cat(paste("Outcome: ", object$outcome, "\n", sep = ""))
  cat(paste("Treatment: ", object$treatment, "\n", sep = ""))
  cat(paste("Mediator: ", object$mediator, "\n", sep = ""))
  cat("----------------------------\n")
  cat(paste(
    "Pre-treatment Covariates: ",
    paste(object$pre.treatment, collapse = ", "),
    "\n", sep = ""))
  cat(paste(
    "Post-treatment Covariates: ",
    paste(object$post.treatment, collapse = ", "),
    "\n", sep = ""))
  cat("----------------------------\n")
  cat(
    paste(
      "Number of units matched to each observation in first stage (mediator): ",
      object$L_m, "\n", sep = ""
    )
  )
  cat(paste(
    "Number of units matched to each observation in second stage (treatment): ",
    object$L_a, "\n", sep = ""
  ))
  cat("----------------------------\n")
  cat("Results:\n")
  if (object$se.type == "Asymptotic") {
    cat(paste(
      "Estimated ACDE: ",
      format(object$estimate, digits = digits), "\n", sep = ""
    ))
    cat(paste(
      "Asymptotic Standard Error: ",
      format(object$std.err, digits = digits), "\n", sep = ""
    ))
    cat(paste(
      "Asymptotic ", format(object$ci.level, digits = digits),
      "% confidence interval: [",
      format(object$conf.low, digits = digits), ", ",
      format(object$conf.high, digits = digits), "]\n", sep = ""
    ))
  } else if (object$se.type == "Bootstrap") {
    cat(paste(
      "Estimated ACDE: ",
      format(object$estimate, digits = digits), "\n", sep = ""
    ))
    cat(paste(
      "Bootstrap Standard Error: ",
      format(object$std.err, digits = digits), "\n", sep = ""
    ))
    cat(paste(
      "Bootstrap ", format(object$ci.level, digits = digits),
      "% confidence interval (percentile): [",
      format(object$conf.low, digits = digits), ", ",
      format(object$conf.low, digits = digits), "]\n", sep = ""
    ))
  } else {
    cat("Something went wrong with the summary - SEs not Asymptotic or Bootstrap\n")
  }
  cat("----------------------------\n")
  cat("Sample sizes (N):\n")
  cat(paste(
    "Treatment = 0, Mediator = 0: ",
    object$sizes$N[
      object$sizes[[object$treatment]] == 0 &
        object$sizes[[object$mediator]] == 0],
    "\n", sep = ""
  ))
  cat(paste(
    "Treatment = 0, Mediator = 1: ",
    object$sizes$N[
      object$sizes[[object$treatment]] == 0 &
        object$sizes[[object$mediator]] == 1],
    "\n", sep = ""
  ))
  cat(paste(
    "Treatment = 1, Mediator = 0: ",
    object$sizes$N[
      object$sizes[[object$treatment]] == 1 &
        object$sizes[[object$mediator]] == 0],
    "\n", sep = ""
  ))
  cat(paste(
    "Treatment = 1, Mediator = 1: ",
    object$sizes$N[
      object$sizes[[object$treatment]] == 1 &
        object$sizes[[object$mediator]] == 1],
    "\n", sep = ""
  ))
  cat("----------------------------\n")

}

#' Balance diagnostics for Telescope Match objects
#'
#' @details Provides matching balance diagnostics for \code{tmatch}
#' objects returned by \code{telescope_match}
#'
#' @param object an object of class \code{tmatch} -- results from a
#'   call to \code{telescope_match}
#' @param vars a formula object containing either the treatment or the
#'   mediator as the dependent variable (which denotes whether
#'   first-stage or second-stage balance diagnostics are returned) and
#'   the covariates for which balance diagnostics are requested as the
#'   independent variables. Each covariate or function of covariates
#'   (e.g. higher-order polynomials or interactions) should be
#'   separated by a +.
#' @param data the data frame used in the call to
#'   \code{telescope_match}
#'
#' @return Returns a data frame with the following columns.
#' \itemize{
#'
#' \item variable: Name of covariate
#'
#' \item before_0: Pre-matching average of the covariate in the
#' mediator == 0 (if first stage balance) or treatment == 0 (if second
#' stage balance) condition
#'
#' \item before_1: Pre-matching average of the covariate in the
#' mediator == 1 (if first stage balance) or treatment == 1 (if second
#' stage balance) condition
#'
#' \item after_0: Post-matching average of the covariate in the
#' mediator == 0 (if first stage balance) or treatment == 0 (if second
#' stage balance) condition
#'
#' \item after_1: Post-matching average of the covariate in the
#' mediator == 1 (if first stage balance) or treatment == 1 (if second
#' stage balance) condition
#'
#' \item before_sd: standard deviation of the outcome (pre-Matching)
#'
#' \item before_diff: Pre-matching covariate difference between
#' mediator arms (if first stage balance) or treatment arms (if second
#' stage balance).
#'
#' \item before_std_diff: Pre-matching standardized covariate
#' difference between mediator arms (if first stage balance) or
#' treatment arms (if second stage balance), Equal to Before_Diff/SD.
#'
#' \item after_diff: Post--matching covariate difference between
#' mediator arms (if first stage balance) or treatment arms (if second
#' stage balance).
#'
#' \item after_std_diff: Post-matching standardized covariate
#' difference between mediator arms (if first stage balance) or
#' treatment arms (if second stage balance), Equal to Before_Diff/SD.
#' }
#'
#' @export
#' @importFrom stats as.formula model.frame weighted.mean

balance.tmatch <- function(object, vars, data) {

  ##################
  ### Sanity checks
  ##################

  ### Is the class a 'tmatch'
  if (class(object) != "tmatch") {
    stop("Error: 'object' not of class 'tmatch'")
  }

  ### Does N of data match number of obs
  if (nrow(data) != length(object$included)) {
    stop("Error: number of rows in data not consistent with object")
  }

  data_tr <- data[[object$treatment]][object$included]
  if (!all.equal(as.numeric(data_tr),
                 as.numeric(object$treatment.vec))) {
    stop("Error: treatment vector in data and object don't match")
  }

  ##########################
  ### Validating the formula
  if (object$treatment == as.character(vars[2])) {
    ### Balance diagnostics for treatment
    get.balance <- "treatment"
  } else if (object$mediator == as.character(vars[2])) {
    ### Balance diagnostics for mediator
    get.balance <- "mediator"
  } else {
    stop("Error: Left-hand side of 'vars' must be either the 'treatment' or the 'mediator' from 'object'")
  }

  ###########################
  ### Validating the data frame

  ## subset to the selected rows used in matching
  data <- data[object$included, ]

  ## Do model.frame first to parse any functions
  mf <- tryCatch({
    model.frame(vars, data)
  },
  error = function(e) {
    stop("Could not extract all variables in 'formula' from data")
  })

  ## Do model.matrix to get functions/interactions/expansions
  X <- tryCatch({
    model.matrix(vars, mf)
  },
  error = function(e) {
    stop("Could not extract all variables in 'formula' from data")
  })

  ## Strip out the "(Intercept) column
  X <- X[, colnames(X) != "(Intercept)"]
  A <- model.response(mf)

  before_0 <- colMeans(X[A == 0, ])
  before_1 <- colMeans(X[A == 1, ])
  before_sd <- apply(X, 2, sd)
  before_diff <- before_1 - before_0
  before_std_diff <- before_diff / before_sd

  ### Generate weights on each observation
  if (get.balance == "mediator") {
    ### If it's the mediator, all M=1 units get KLm+1, all M=0 get KLm
    obs.weights <- (object$KLm / object$L_m + 1) * A +
      (object$KLm / object$L_m) * (1 - A)
  } else if (get.balance == "treatment") {
    ### If it's the treatment, all A=1 units get KLa+1, all A=0 get KLa+1
    obs.weights <- (object$KLa / object$L_a + 1)
  }

  ## get weighted means of each variable
  after_0 <- apply(X[A == 0, ], 2, weighted.mean, obs.weights[A == 0])
  after_1 <- apply(X[A == 1, ], 2, weighted.mean, obs.weights[A == 1])
  after_diff <- after_1 - after_0
  after_std_diff <- after_diff / before_sd

  bal_tab <- data.frame(
    variable = colnames(X),
    before_0 = before_0, before_1 = before_1,
    after_0 = after_0, after_1 = after_1,
    before_sd = before_sd,
    before_diff = before_diff, before_std_diff,
    after_diff = after_diff, after_std_diff
  )


  ## Output
  return(bal_tab)
}

#' Histograms of matching weights
#'
#' @details Provides histograms of the number of times each unit is
#' used as a match given a  \code{tmatch} object returned by
#' \code{telescope_match}
#'
#' @param object an object of class \code{tmatch} -- results from a
#' call  to \code{telescope_match}
#' @param stage a character vector equal to either 'mediator' or
#' 'treatment'. If equal to 'mediator', returns a histogram
#' of matching weights for units with mediator = 0. If equal to
#' 'treatment', returns a histogram of matching weights for
#' all units.
#'
#' @returns Outputs a `plot()` object containing the histogram of match counts
#'
#' @export
#' @importFrom graphics hist

plotDiag.tmatch <- function(object, stage = "mediator") {

  ##################
  ### Sanity checks
  ##################

  ### Is the class a 'tmatch'
  if (class(object) != "tmatch") {
    stop("Error: 'object' not of class 'tmatch'")
  }

  if (stage != "mediator" & stage != "treatment") {
    stop("Error: 'stage' must be either 'mediator' or 'treatment'.")
  }

  #######################
  ### Main output
  #######################

  ### If it's mediator, plot the histogram of K_Lm for mediator == 0
  if (stage == "mediator") {
    plot_title <- paste("Matching weights for first stage (mediator)\namong M = 0")
    hist(
      object$KLm[object$mediator.vec == 0] / object$L_m,
      main = plot_title,
      xlab = "Number of times unit is matched"
    )
    abline(v = 1, col = "red", lty = 2, lwd = 2)
  } else if (stage == "treatment") {
    plot_title <- paste("Matching weights for second stage (treatment)")
    hist(
      object$KLa / object$L_a,
      main = plot_title,
      xlab = "Number of times unit is matched"
    )
    abline(v = 1, col = "red", lty = 2, lwd = 2)
  }

  return()

}

##' Calculate telescope matching weights
##'
##' Calculates the weight given to each unit in the telescope matching
##' estimation procedure.
##'
##' @param object an object of class `tmatch` from [telescope_match]
##' output.
##' @param data a data.frame that was passed to [telescope_match] to
##' produce `object`.
##' @return A list containing two vectors:
##' \item{treat_w}{A n-length vector containing the matching weights
##' for the treatment or first period matching.}
##' \item{med_w}{A n-length vector containing the matching weights
##' for the mediator or second period matching.}
##' @author Matthew Blackwell
##' @export
tmatch_weights <- function(object, data) {
  ### Is the class a 'tmatch'
  if (class(object) != "tmatch") {
    stop("Error: 'object' not of class 'tmatch'")
  }

  ### Does N of data match number of obs
  if (nrow(data) != length(object$included)) {
    stop("Error: number of rows in data not consistent with object")
  }

  data_tr <- data[[object$treatment]][object$included]
  if (!all.equal(as.numeric(data_tr),
                 as.numeric(object$treatment.vec))) {
    stop("Error: treatment vector in data and object don't match")
  }
  a_w <- rep(NA, nrow(data))
  m_w <- rep(NA, nrow(data))

  M <- object$mediator.vec
  a_w[object$included] <- (object$KLa / object$L_a + 1)
  m_w[object$included] <- (object$KLm / object$L_m + 1) * M +
    (object$KLm / object$L_m) * (1 - M)

  return(list(treat_w = a_w, med_w = m_w))
}
