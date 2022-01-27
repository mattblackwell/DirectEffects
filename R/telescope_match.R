#' Perform telescope matching to estimate the controlled
#' direct effect of a binary treatment net the effect of binary
#' mediators
#'
#' @param formula A formula object that specifies the covariates and
#' treatment variables (or mediators) in causal ordering from oldest
#' to newest with each group separated by \code{|}. See below for
#'   more details.
#' @param data A dataframe containing variables referenced by
#'   \code{formula}.
#' @param caliper A scalar denoting the caliper to be used in matching
#'   in the treatment stage (calipers cannot be used for matching on
#'   the mediator). Observations outside of the caliper are dropped.
#'   Calipers are specified in standard deviations of the covariates.
#'   NULL by default (no caliper).
#' @param L Number of matches to use for each unit. Must be a numeric
#'   vector of either length 1 or 2. If length 1, L sets the number of
#'   matches used in both the first stage (matching on mediator) and
#'   in the second stage (matching on treatment). If length 2, the
#'   first element sets the number of matches used in the first stage
#'   (matching on mediator) and the second element sets the number of
#'   matches used in the second stage (matching on treatment) Default
#'   is 5.
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
#'   and Strezhnev (2021).
#'
#'  The procedure first estimates a demediated outcome using a
#'  combination of matching and a regression bias-correction. The
#'  \code{data.frame} passed to \code{data} should be in the wide
#'  format so that each row corresponds to a single unit and
#'  treatments and covariates from different time periods appear as
#'  different columns. The \code{formula} argument specifies both the
#'  causal ordering of the variables and the regression specifications
#'  for the bias correction. It should be of the form \code{Y ~ X1 |
#'  A1 | X2 | A2}, where \code{Y} is the outcome, \code{X1} is a
#'  formula of baseline covariates, \code{A1} is a single variable
#'  name indicating the binary treatment in the first period,
#'  \code{X2} is a formula of covariates in period 2, and \code{A2} is
#'  a single variable name indicating treatment in period 2 (which is
#'  also sometimes called the mediator). Note that it is possible to
#'  add more covariate/treatment pairs for additional time periods.
#'
#' Under the default \code{separate_bc == TRUE}, the function will
#' match for each treatment/mediator based on the the covariates up to
#' that point within levels of past treatments (so for \code{A2} this
#' matching finds units with similar values of \code{X1} and \code{X2}
#' and the same value of \code{A1}). Once this matching is complete,
#' the function moves backward through treatments and imputes
#' potential outcomes using matches and bias-correction regressions,
#' which regress the current imputed potential outcome on the past
#' covariates, within levels of the treatment history up to the
#' current period. The functional form comes from the specification in
#' \code{formula}. Controlled direct effects of \code{A1} are
#' estimated for every possible combination of future treatments.
#'
#' When \code{separate_bc} is \code{FALSE}, the bias correction
#' regressions are not broken out by the treatments/mediators and
#' those variables are simply included as separate regressors as
#' specified in \code{formula}. In this setting, interactions between
#' the treatment/mediator and covariates can be added on a selective
#' basis to the covariate block (\code{X1} or \code{X2} and so on)
#' specifications.
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
#'   \item \code{call}: the matched call.
#'   \item \code{formula}: formula used to fit the model.
#'   \item \code{m_out}: list of matching solutions at each time
#'     point. Each member of the list has a `matches` list giving the
#'     units matched to that unit, a `donors` list with the units to
#'     which the unit is matched, and a `tr` vector which is just the
#'     treatment vector being matched.
#'   \item \code{K}: data.frame of indicating how many times a unit
#'     has been used as a match, directly in each period and indirectly
#'     across periods.
#'   \item \code{L}: vector of matching ratios used in each period.
#'   \item \code{r_out}: nested list of regression imputations used in the
#'     bias correction. The first level of the list varies across
#'     different controlled direct effects (different sequences of future
#'     treatments/mediators). Each of these is a list of time periods and
#'     each of these time periods is a list of `yhat_r_0` and `yhat_r_1`
#'     that give the regression predictions for the potential outcomes at
#'     that time point when the treatment at that time point is 0 or 1,
#'     respectively, along with `n_coefs` giving the number of
#'     coefficients estimated in those models.
#'   \item \code{tau}: vector of bias-corrected estimates of average
#'     controlled direct effects for different vectors of future
#'     treatments/mediators.
#'   \item \code{tau_raw}: vector of standard matching estimates of average
#'     controlled direct effects for different vectors of future
#'     treatments/mediators without using bias correction.
#'   \item \code{tau_se}: vector of estimated standard errors for the average
#'     controlled direct effects estimates for different vectors of future
#'     treatments/mediators. 
#'   \item \code{tau_i}: matrix of individuals contributions to the ACDE
#'     estimates (units on rows, different ACDEs on columns). Used for
#'     weighted bootstrap.
#'   \item \code{included}: logical vector indicating if each row of
#'     \code{data} was included in estimating \code{tau}.
#'   \item \code{effects}: data frame where each row describes the
#'     different ACDEs in \code{tau}. The \code{active} column describes
#'     the which variable's direct effect is being assessed and the rest
#'     of the columns describe the fixed values of the future
#'     treatments/mediators for that ACDE.
#'   \item \code{a_names}: character vector with the names of the
#'     treatment/mediator variables used in estimation.
#'   \item \code{caliper}: caliper (if any) used in matching to drop
#'     distant observations.
#' }

#' @references Blackwell, Matthew, and Strezhnev, Anton (2020)
#' "Telescope Matching: Reducing Model Dependence
#' in the Estimation of Direct Effects." Journal of the Royal
#' Statistical Society (Series A). \doi{10.1111/rssa.12759}
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
#' @importFrom utils combn
#'
telescope_match <- function(formula, data, caliper = NULL, L = 5,
                            verbose = TRUE, subset, contrasts = NULL,
                            separate_bc = TRUE, ...) {

  ########################
  ### Quick pre-processing
  ########################


  cl <- match.call(expand.dots = TRUE)
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
      "\n Groups should be divded by `|`", call. = FALSE
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
  included <- rownames(data) %in% rownames(mf)
  
  all_data <- model.matrix(mt, mf, contrasts)
  Y <- model.response(mf, "numeric")

  N <- length(Y)

  X <- list()
  drops <- list()
  A <- matrix(NA, nrow = N, ncol = T)
  a_names <- rep("", times = T)

  x_spots <- seq(1, T * 2, by = 2)
  a_spots <- seq(2, T * 2, by = 2)
  for (s in seq_len(T)) {
    ## covariates up to s
    mt_a <- terms(formula, data = data, rhs = a_spots[s])
    if (length(attr(mt_a, "term.labels")) != 1) {
      stop("Error: there must only be single treatment in each period.",
           call. = FALSE)
    }
    a_names[s] <- attr(mt_a, "term.labels")
    A[, s] <- all_data[, a_names[s]]
    mt_x <- terms(formula, data = data, rhs = x_spots[1:s])
    X[[s]] <- model.matrix(mt_x, mf, contrasts)

    ## treatnent in covariate specifications. need to dropped for
    ## matching
    a_in_x <- which(rownames(attributes(mt_x)$factors) %in% a_names[1:s])
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
  r_out_big <- list()
  K <- list()
  K_paths <- list()

  ## matching
  if (verbose) message("Beginning matching...")
  for (s in T:1) {
    if (verbose) message("Matching ", a_names[s], "...")
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
    K_paths[[s]] <- unlist(
      lapply(
        paths,
        function(x) paste0(a_names[x], collapse = ":")
      )
    )
    
    names(K[[s]]) <- unlist(lapply(paths, paste0, collapse = ":"))
    names(K[[s]]) <- paste0("K_", names(K[[s]]))
  }

  K <- as.data.frame(unlist(K, recursive = FALSE))
  
  ## bias correction
  if (verbose) message("Beginning bias correction...")
  A_levs <- expand.grid(rep(list(c(0, 1)), times = T - 1))
  tau_i <- matrix(NA, nrow = N, ncol = nrow(A_levs))
  tau_raw <- rep(0, times = nrow(A_levs))
  tau_se <- rep(0, times = nrow(A_levs))
  mu_hat <- matrix(NA, nrow = N, ncol = T)
  Yt <- matrix(NA, nrow = N, ncol = T * nrow(A_levs))
  for (j in seq_len(nrow(A_levs))) {
    A_j <- A_levs[j, ]
    Ytilde <- Y
    r_out <- list()
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

    cdes <- calculate_cdes(Y, A, K, mu_hat, r_out, A_j)
    tau_i[, j] <- cdes$tau_i
    tau_raw[j] <- cdes$tau_raw
    tau_se[j] <- cdes$tau_se
    r_out_big[[j]] <- r_out
  }
  tau <- colMeans((2 * A[, 1] - 1) * tau_i)
  
  names(K) <- unlist(K_paths)

  effects <- data.frame(active = a_names[1], A_levs)
  names(effects)[-1] <- a_names[-1]
  
  ### Return output
  out <- list(call = cl, formula = formula, m_out = m_out, K = K, L = L, 
              r_out = r_out_big, tau = tau, tau_raw = tau_raw, tau_se = tau_se,
              tau_i = tau_i, included = included, effects = effects,
              a_names = a_names)

  class(out) <- "tmatch"
  return(out)

}


## given the matching output and a set of periods `path = (k, ..., j)`
## find the number of times units are used indirectly as a match in
## period `k` through the given path. E.g., if `path = c(1, 2 ,3)`,
## then this returns the number of times each unit is matched to a
## unit j in period 3, j is matched to h in period 2, and h is matched
## to g in period 1, where j, h, and g are all different.

## old tm code used row names, new assembly line code uses positions. 
find_indirect_matches <- function(m_out, path, N) {
  K <- rep(0, times = N)
  path <- sort(path, decreasing = TRUE)
  d <- m_out[[path[1]]]$donors
  path <- path[-1]
  for (j in path) {   
    if (length(names(d))) {
      d <- lapply(d, function(x) unlist(m_out[[j]]$donors[as.character(x)]))
    } else {
      d <- lapply(d, function(x) unlist(m_out[[j]]$donors[x]))
    }
  }
  n_ind_matches <- lapply(d, length)
  if (length(names(n_ind_matches))) {
    pos <- as.numeric(names(n_ind_matches))
  } else {
    pos <- 1L:N
  }

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
      "\n treatment ", k, " non-binary values."
    )
  }
  X_m <- cbind(X[, -drops], A[, -k])

  exacts <- c(
    rep(FALSE, times = ncol(X[, -drops, drop = FALSE])),
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

  out <- list(matches = matches, donors = donors, tr = Ak)
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

calculate_cdes <- function(Y, A, K, mu_hat, r_out, A_j) {
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

  ## we have to adjust for the fact that we run regressions for all
  ## treatment combinations, but the CDE is only for one sequence
  ## (for t>1)
  p <- sapply(r_out, function(x) x$n_coefs)
  p <- p / (2 ^ (0:(T - 1)))
  
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
      i_j <- rep(1, N)
    }
    A_jt <- as.numeric(A[, s] == A_j[s - 1])

    ## calculate K weight for this period
    fut_patt <- paste0("[", paste0(s:T, collapse = ":"), "]")
    past_K <- grep(fut_patt, names(K), invert = TRUE)
    fut_K <- grep(fut_patt, names(K))
    K_term <- (A_jt * rowSums(K[fut_K]) - (1 - A_jt) * (1 + rowSums(K[past_K])))

    eta <- mu_hat[, s] - mu_hat[, s - 1]
    tau_i <- tau_i - as.numeric(i_j) * K_term * eta

    mu_lag_A <- A[, s - 1] * r_out[[s - 1]]$yhat_r_1 +
      (1 - A[, s - 1]) * r_out[[s - 1]]$yhat_r_0
    eta_v <- mu_hat[, s] - mu_lag_A
    dfc_s <- sum(i_j) / (sum(i_j) - sum(p[(s - 1) : T]))
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


#' Bootstrap Uncertainty Estimates for Telescope Matching
#' 
#' Performs a weighted bootstrap procedure for the output of
#' \code{\link{telescope_match}}. 
#'
#' @param obj A \code{tmatch} object, computed by \code{\link{telescope_match}}.
#' @param boots The number of bootstrap replicates. Defaults to 1000.
#' @param ci_alpha alpha value for the bootstrapped confidence
#' intervals. Corresponds to a 100 * (1-alpha) confidence interval. 
#'
#' @return An data.frame with columns `ci_low` and `ci_high` which
#' contain the bootstrapped confidence intervals for the estimated
#' ACDEs in \code{obj$tau}.
#'
#' @examples
#' \donttest{
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
#' out.boots <- boots_tm(tm_out)
#'
#' out.boots
#' }
#' @export
#' @importFrom stats rbinom quantile
#'
boots_tm <- function(obj, boots = 1000, ci_alpha = 0.05) {

  N <- nrow(obj$tau_i)
  J <- length(obj$tau)
  ## Bootstrap iterations
    W.bern <- sapply(
      1:boots,
      function(x) rbinom(N, 1, prob = (sqrt(5) - 1) / (2 * sqrt(5)))
    )
    Wstar <- (((sqrt(5) + 1) / 2) * W.bern +
                ((-sqrt(5) + 1) / 2) * (1 - W.bern)) / N

  ci_high <- ci_low <- rep(NA, times = J)
  for (j in seq_len(J)) {
    tau_norm <- obj$tau_i[, j]  - obj$tau[j]
    Tstar <- (t(Wstar) %*% tau_norm)
    Tstar <- Tstar + obj$tau[j]
    ci_low[j] <- quantile(Tstar, ci_alpha / 2)
    ci_high[j] <- quantile(Tstar, 1 - ci_alpha / 2)
  }

  return(data.frame(ci_low = ci_low, ci_high = ci_high))
}


#' @export 
print.tmatch <- function(x, ...) {
  cat("\nTelescope matching output\n\n")
  cat("Call:\n")
  print(x$call)

  cat("\nActive treatment: ", x$effects$active[1], "\n")
  cat("Controlled treatment(s):", names(x$effects)[-1])

  cat(
    "\n\nEstimated controlled direct effects of ",
    x$effects$active[1],
    ":\n", sep = ""
  )
  df <- data.frame(
    x$effects[, -1, drop = FALSE],
    estimate = x$tau)
  print(df)
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
#' @details Returns a summary data frame containing the estimate and
#' standard errors from the `telescope_match` object.
#'
#' @return Returns an object of \code{class} \code{summary.tmatch}.
#' Contains the following components
#' \itemize{
#'   \item \code{call}: matched call.
#'   \item \code{m_summary}: data.frame summarizes the matching
#'     ratios (\code{{ratio}}), number of units \code{n_1, n_0},
#'     and number of matched units (\code{matched_1, matched_0}) for
#'     each treatment/mediator (\code{term}).
#'   \item \code{K}: \code{K} data frame from the \code{object}
#'     telescope matching output.
#'   \item \code{L}: \code{L} vector from the \code{object}
#'     telescope matching output.
#'   \item \code{a_names}: character vector of the names of the
#'     treatment/mediator variables used in matching. 
#'   \item \code{estimates}: matrix of estimated ACDEs with and
#'     without bias correction and the estimated standard errors.
#' }
#' @export
summary.tmatch <- function(object, ...) {

  
  an <- object$a_names
  J <- length(an)
  sum_mat <- matrix(NA, nrow = J, ncol = 4)
  for (j in 1:J) {
    tr <- object$m_out[[j]]$tr
    sum_mat[j, ] <- c(
      sum(tr), sum(1 - tr),
      sum(object$K[[an[j]]] * tr > 0),
      sum(object$K[[an[j]]] * (1 - tr) > 0)
    )
  }
  colnames(sum_mat) <- c(
    "n_1", "n_0", "matched_1", "matched_0"
  )

  
  out <- list()
  out$call <- object$call
  out$m_summary <- data.frame(term = an, ratio = object$L, sum_mat)

  out$K <- object$K
  out$L <- object$L
  out$a_names <- an

  cf <- cbind(
    object$effects[, -1, drop  = FALSE],
    "Estimate" = object$tau,
    "Estimate (no BC)" = object$tau_raw,
    "Std. Err." = object$tau_se
  )

  out$estimates <- as.matrix(cf)

  A_levs <- expand.grid(rep(list(c(0, 1)), times = J - 1))
  lab_1 <- paste0("(1, ", apply(A_levs, 1, paste0, collapse = ", "), ")")
  lab_0 <- paste0("(0, ", apply(A_levs, 1, paste0, collapse = ", "), ")")
  labs <- paste0(lab_1, " vs. ", lab_0)
  rownames(out$estimates) <- labs
  
  class(out) <- "summary.tmatch"

  return(out)
}

#' @export
print.summary.tmatch <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("Telescope matching results\n")
  cat("----------------------------\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nActive treatment: ", x$a_names[1], "\n")
  cat("Controlled treatment(s):", x$a_names[-1], "\n\n")

  cat("Matching summary:\n")
  m_sum <- x$m_summary
  names(m_sum) <- c(
    "Term", "Matching Ratio L:1", "N == 1", "N == 0",
    "Matched == 1", "Matched == 0"
  )
  print(m_sum)

  cat("\n\nSummary of units matching contributions:\n")
  print(t(sapply(x$K, summary)))
  
  cat(
    "\n\nEstimated controlled direct effects of ",
    x$a_names[1],
    ":\n", sep = ""
  )

  print(
    x$estimates,
    digits = digits, quote = FALSE, right = TRUE
  )
  cat("\n")
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
#' @param comparison a binary indicator for if the function should
#' return the balance for the treated group (`1`), for the control
#' group (`0`), or for overall combined balanced (`NULL`, the
#' default). 
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

balance.tmatch <- function(object, vars, data, comparison = NULL) {

  ##################
  ### Sanity checks
  ##################

  ### Is the class a 'tmatch'
  if (class(object) != "tmatch") {
    stop("`object` not of class 'tmatch'", call. = FALSE)
  }

  ### Does N of data match number of obs
  if (nrow(data) != length(object$included)) {
    stop(
      "number of rows in data not consistent with object",
      .call = FALSE
    )
  }

  ##########################
  ### Validating the formula
  a_names <- names(object$K)

  if (!(as.character(vars[2]) %in% a_names)) {
    stop(
      "Left-hand side of `vars` must be a treatment variable from `object`.",
      call. = FALSE
    )
  }

  this_K <- object$K[, a_names == as.character(vars[2])]

  
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

  if (!is.null(comparison)) {
    obs.weights <- comparison * (this_K * (1 - A) + A) +
      (1 - comparison) * (this_K * A + (1 - A))
  } else {
    obs.weights <- this_K + 1
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
#' used as a match given a \code{tmatch} object returned by
#' \code{telescope_match}
#'
#' @param object an object of class \code{tmatch} -- results from a
#' call  to \code{telescope_match}
#' @param stage a character vector equal to the name of one treatment
#' from the `object`. 
#'
#' @returns Outputs a `plot()` object containing the histogram of match counts
#'
#' @export
#' @importFrom graphics hist

plotDiag.tmatch <- function(object, stage) {

  ##################
  ### Sanity checks
  ##################

  ### Is the class a 'tmatch'
  if (class(object) != "tmatch") {
    stop("`object` not of class 'tmatch'", call. = FALSE)
  }


  a_names <- names(object$K)
  if (missing(stage)) stage <- a_names[1]
  
  if (!(stage %in% a_names)) {
    stop("`stage` must be a treatment in `object`.", call. = FALSE)
  }

  #######################
  ### Main output
  #######################

  this_K <- object$K[, a_names == stage]
  
  plot_title <- paste0("Matching weights for ", stage)
  hist(
    this_K,
    main = plot_title,
    xlab = "Number of times unit is matched"
  )
  abline(v = 1, col = "red", lty = 2, lwd = 2)

  invisible()

}


#' @export
cde_telescope_match <- function(...) {
  args <- rlang::enquos(...)

  new_cde_estimator(
    "telescope_match",
    args = args,
    model_spec = NULL
  )
}


t_match <- function(object, data) {

  fit_env <- rlang::env()  
  A <- get_treat_df(object, data)
  paths <- interaction(A, sep = "_")
  Y <- get_outcome(object, data)
  num_treat <- length(object$model_spec)
  N <- length(Y)
  
  ## move backward through blocks
  block_seq <- rev(seq_len(num_treat))

  
  out <- list()
  out$match_out <- make_match_holder(A)
  L <- unlist(lapply(object$model_spec, function(x) rlang::eval_tidy(x$treat_spec$args$L)))
  for (j in block_seq) {
    if (j > 1) {
      past_fit <- interaction(A[, 1:(j - 1), drop = FALSE], sep = "_")
      past_splits <- split(seq_len(nrow(A)), past_fit)
    } else {
      past_fit <- rep(0, times = nrow(A))
      past_splits <- list(seq_len(nrow(A)))
    }
    M <- length(past_splits)

    ## TODO: add check about overlap?
    for (m in seq_len(M)) {
      past_strata_rows <- past_splits[[m]]
      fit_env$fit_data <- data[past_strata_rows, ]
      match_out <- match_model(object$model_spec[[j]]$treat_spec, fit_env)
      match_out$matches <- lapply(
        match_out$matches,
        function(x) past_strata_rows[x]
      )
      match_out$donors <- lapply(
        match_out$donors,
        function(x) past_strata_rows[x]
      )
      
      match_rows <- past_strata_rows[as.numeric(names(match_out$matches))]
      out$match_out[[j]]$matches[match_rows] <- match_out$matches
      donor_rows <- past_strata_rows[as.numeric(names(match_out$donors))]
      out$match_out[[j]]$donors[donor_rows] <- match_out$donors
      out$match_out[[j]]$L <- match_out$L

    }
    out$match_out[[j]]$K <- compute_match_weights(object, out$match_out, j)
  }
  out$match_out
}


match_model <- function(model, fit_env, ...) {

  fit_env$X <- model.matrix(model$formula, data = fit_env$fit_data)
  fit_env$X <- fit_env$X[, which(colnames(fit_env$X) != "(Intercept)")]
  fit_env$tr <- model.response(model.frame(model$formula, fit_env$fit_data))
  args <- as.list(model$args)
  args$X <- quote(X)
  args$Tr <- quote(tr)
  names(args)[names(args) == "L"] <- "M"
  args$Weight <- 2
  args$ties <- FALSE
  att_call <- rlang::call2("Match", !!!args, estimand = "ATT", .ns = "Matching")
  atc_call <- rlang::call2("Match", !!!args, estimand = "ATC", .ns = "Matching")
  tm_att <- rlang::eval_tidy(att_call, env = fit_env)
  tm_atc <- rlang::eval_tidy(atc_call, env = fit_env)

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
  return(list(donors = donors, matches = matches, L = rlang::eval_tidy(args$M)))

}

compute_match_weights <- function(object, out, j) {
  num_treat <- length(object$model_spec)
  N <- length(out[[num_treat]]$matches)
  a_names <- unlist(lapply(
    object$model_spec,
    function(x) rlang::as_label(x$treat)
  ))
  path_lengths <- 1:(num_treat - j + 1)
  paths <- lapply(
    path_lengths,
    function(x) combn(j:num_treat, m = x, simplify = FALSE)
  )
  paths <- unlist(paths, recursive = FALSE)
  paths <- paths[unlist(lapply(paths, function(x) x[1] == j))]
  L <- unlist(lapply(
    object$model_spec,
    function(x) rlang::eval_tidy(x$treat_spec$args$L)
  ))
  K <- lapply(
    paths,
    function(x) find_indirect_matches(out, path = x, N) / prod(L[x])
  )
  K_paths <- unlist(
    lapply(
      paths,
      function(x) paste0(a_names[x], collapse = ":")
    )
  )
    
  names(K) <- unlist(lapply(paths, paste0, collapse = ":"))
  names(K) <- paste0("K_", names(K))
  K
}

make_match_holder <- function(A) {
  out <- vector("list", ncol(A))
  for (j in seq_along(out)) {
    out[[j]]$matches <- vector("list", nrow(A))
    out[[j]]$donors <- vector("list", nrow(A))
  }
  out
}

compute_telescope_match <- function(j, j_levs, y, treat, out, args, term_name, dfs) {
  num_treat <- length(out$outreg_pred)
  N <- length(treat)
  j_levs <- sort(j_levs)
  paths <- colnames(out$outreg_pred[[j]])
  sp <- strsplit(paths, "_")
  templates <- unique(replace_each(sp, j, NA))
  A <- strsplit(as.character(treat), "_", fixed = TRUE)
  A <- do.call(rbind, lapply(A, as.numeric))
  K <- lapply(out$match_out[j:num_treat], function(x) as.data.frame(x$K))
  K <- do.call(cbind, K)
  W <- 1 + rowSums(K)
  sm_part <- W * y

  est_tab <- empty_est_tab()
  num_fut_strata <- 2 ^ rev(seq(j, num_treat)) - j + 1
  dfs[j:num_treat] <- dfs[j:num_treat] * num_fut_strata
  for (k in seq_along(templates)) {
    base <- templates[[k]]
    base[j] <- j_levs[1L]
    
    base <- paste0(base, collapse = "_")
    ctr <- as.numeric(treat == base)
    N_c <- sum(ctr)
    r_ctr <- get_reg_preds(out, base)
    A_ctr <- get_path_inds(treat, base)
    eps_ctr <- cbind(r_ctr[, -1L, drop = FALSE], y) - r_ctr
    
    plus <- templates[[k]]
    plus[j] <- j_levs[2L]
    plus <- paste0(plus, collapse = "_")
    trt <- as.numeric(treat == plus)
    N_t <- sum(trt)
    r_trt <- get_reg_preds(out, plus)
    A_trt <- get_path_inds(treat, plus)
    eps_trt <- cbind(r_trt[, -1L, drop = FALSE], y) - r_trt
    
    A_j <- A_trt[, j]
    if (j > 1) {
      strata_past <- paste0(templates[[k]][1:(j - 1)], collapse = "_")      
      past <- apply(A[, 1:(j - 1), drop = FALSE], 1, paste0, collapse = "_")
      this_past <- as.numeric(past == strata_past)
    } else {
      this_past <- rep(1, N)
    }
    
    
    strata_rows <- as.numeric(treat %in% c(base, plus))
    N_s <- sum(strata_rows)
    tau_i1 <- trt * sm_part
    tau_i0 <- ctr * sm_part
    tau_i <- tau_i1 - tau_i0
    tau_raw <- mean(tau_i[this_past == 1L])

    K_term <- strata_rows * rowSums(K) - (1 - strata_rows)
    B_0 <- A_j * (r_ctr[, j]  + K_term * r_trt[, j])
    B_0 <- B_0 - (1 - A_j) * (r_trt[, j] + K_term * r_ctr[, j])
    tau_i <- tau_i - B_0

    dfc_1 <- N_s / (N_s - dfs[num_treat])
    tau_v <- W ^ 2 * (trt * eps_trt[, num_treat] ^ 2 + ctr * eps_ctr[, num_treat] ^ 2)
    tau_var <- dfc_1 * mean(tau_v[this_past == 1L])
    if (j == num_treat) {
      fut <- numeric(0)
    } else {
      fut <- (j + 1):num_treat
    }

    
    
    for (s in fut) {
      A_fut <- A[, s:num_treat, drop = FALSE]
      A_fut_str <-  apply(A_fut, 1, paste0, collapse = "_")
      strata_fut <- paste0(templates[[k]][s:num_treat], collapse = "_")
      D <- as.numeric(A_fut_str == strata_fut)

      fut_patt <- paste0("[", paste0(s:num_treat, collapse = ":"), "]")
      past_K <- grep(fut_patt, names(K), invert = TRUE)
      fut_K <- grep(fut_patt, names(K))
      W_s <- D * (rowSums(K[fut_K])) - (1 - D) * (1 + rowSums(K[past_K]))

      bias_s <- A_j * eps_trt[, s - 1] * W_s - (1 - A_j) * eps_ctr[, s - 1] * W_s
      tau_i <- tau_i - bias_s

      dfc_s <- sum(D * this_past) / (sum(D * this_past) - sum(dfs[(s - 1):num_treat]))
      tau_v <- dfc_s * (W - rowSums(K[fut_K])) ^ 2 * (A_j * eps_trt[, s - 1] ^ 2 + (1 - A_j) * eps_ctr[, s - 1] ^ 2)
      tau_var <- tau_var + mean(tau_v[this_past == 1L])
    }
    
    tau_est <- mean(tau_i[this_past == 1])
    dfc_x <- sum(this_past) / (sum(this_past) - sum(dfs) - 1)
    tau_x <- dfc_x * (r_trt[, j] - r_ctr[, j] - tau_est) ^ 2
    tau_var <- tau_var + mean(tau_x[this_past == 1L])
    this_est <- data.frame(
        term = term_name,
        block_num = j,
        active = format_path(plus),
        control = format_path(base),
        estimate = tau_est,
        std_err = sqrt(tau_var / sum(this_past))
    )
      est_tab <- rbind(est_tab, this_est)

  }
  rownames(est_tab)
  est_tab

}
