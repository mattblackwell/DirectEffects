add_class <- function(x, cls) {
  class(x) <- unique(c(cls, class(x)))
  x
}

direct_effects <- rlang::new_environment()
direct_effects$models <- NULL

get_de_env <- function () {
  e <- utils::getFromNamespace("direct_effects", ns = "DirectEffects")
  e
}

get_from_env <- function(keys) {
  de_env <- get_de_env()
  rlang::env_get(de_env, keys, default = NULL)
}

bind_in_env <- function(...) {
  de_env <- get_de_env()
  rlang::env_bind(de_env, ...)
}

set_env <- function(key, val) {
  if (length(key) != 1L || !is.character(key)) {
    rlang::abort("`key` can only be a single character value.")
  }
  de_env <- get_de_env()
  x <- val
  names(x) <- key
  rlang::env_bind(de_env, !!!x)
}



get_treat_df <- function(object, data) {
  
  model_spec <- object$model_spec
  vars <- as.list(rlang::set_names(seq_along(data), names(data)))
  locs <- lapply(model_spec, function(x) rlang::eval_tidy(x$treat, vars))
  locs <- unlist(locs)

  data[, locs, drop = FALSE]
}

get_outcome <- function(object, data) {
  vars <- as.list(rlang::set_names(seq_along(data), names(data)))
  locs <- rlang::eval_tidy(object$outcome, vars)
  data[[locs]]
}

make_pred_holder <- function(A, type) {
  J <- length(A)
  out <- vector(mode = "list", length = J)
  N <- nrow(A)  
  for (j in seq_len(J)) {
    if (type == "ipw") {
      levs <- levels(interaction(A[, 1L:j, drop = FALSE], sep = "_"))
    } else if (type == "outreg") {
      levs <- levels(interaction(A, sep = "_"))
    } else {
      rlang::abort("type of prediction holder incorrect")
    }

    num_past <- length(levs)
    out[[j]] <- matrix(NA, nrow = N, ncol = num_past)
    colnames(out[[j]]) <- levs
  }
  out
}


make_pred_holder <- function(A, eval_vals, j, type, engine_type, separate) {
  N <- nrow(A)
  num_treat <- length(eval_vals)
  
  if (type == "treat") {
    if (engine_type == "categorical") { ## classification
      levs <- apply(expand.grid(eval_vals[1L:j]), 1, paste0, collapse = "_")
    } else {
      levs <- apply(expand.grid(eval_vals[1L:(j - 1)]), 1, paste0, collapse = "_")
    }
  } else if (type == "outreg") {
    levs <- apply(expand.grid(eval_vals), 1, paste0, collapse = "_")
  } else {
    rlang::abort("type of prediction holder incorrect")
  }
  num_cols <- length(levs)
  out <- matrix(NA, nrow = N, ncol = num_cols)
  colnames(out) <- levs
  

  out
}

make_model_fits <- function(A, eval_vals, model_spec) {
  model_fits <- vector("list", length = length(model_spec))
  for (j in seq_along(model_spec)) {
    block <- model_spec[[j]]
    
    if (!is.null(block$treat_spec)) {
      model_fits[[j]]$treat_pred <- make_pred_holder(
        A, eval_vals, j, "treat",
        block$treat_type,
        block$treat_spec$separate
      )
    }
    if (!is.null(block$outreg_spec)) {
      model_fits[[j]]$outreg_pred <- make_pred_holder(
        A, eval_vals, j, "outreg",
        block$treat_type,
        block$outreg_spec$separate
      )
      model_fits[[j]]$blip_pred <- make_pred_holder(
        A, eval_vals, j, "outreg",
        block$treat_type,
        block$outreg_spec$separate
      )
    }
  }
  model_fits
}


split_length <- function(x, split) {
  unlist(lapply(strsplit(x, split), length))
}

make_effect_labels <- function(labs, j) {
  splits <- strsplit(labs, "_")
  out <- unlist(lapply(
    splits,
    function(x) {
      x[j] <- "*"
      paste0(x, collapse = "_")
    }
  ))
  out
}

replace_each <- function(x, ind, val) {
  lapply(
    x,
    function(x) {
      x[ind] <- val
      x
    }
  )
}

format_path <- function(paths) {
  paths <- strsplit(paths, "_")
  out <- lapply(
    paths,
    function(x) {
      int <- paste0(x, collapse = ", ")
      paste0("(", int, ")")
    }
  )
  out <- unlist(out)
  out
}

format_term <- function(term, template, plus, base, var_names) {

  vars <- var_names[!is.na(template)]
  levs <- template[!is.na(template)]
  j_str <- paste0(plus, " vs ", base)
  lev_str <- paste0(vars, " = ", levs, collapse = ", ")
  out <- paste0(term, " [", j_str, ", ", lev_str, "]")
  out
}

empty_est_tab <- function() {  
  out <- data.frame(
    term = character(0),
    active = character(0),
    control = character(0),
    estimate = numeric(0),
    std.error = numeric(0),
    DF = numeric(0)
  )
}

check_cde_estimator <- function(object) {

  
  ## check separate coherency.
  if (object$has_ipw) {
      treat_seps <- unlist(lapply(
        object$model_spec,
        function(x) x$treat_spec$separate
      ))
      treat_engine_types <- unlist(lapply(
        object$model_spec,
        function(x) x$treat_spec$engine_type
      ))
      if (any(treat_engine_types == "reg")) {
        first_reg <- which(treat_engine_types == "reg")[1L]
        if (any(treat_seps[first_reg:length(treat_seps)])) {
          rlang::abort("cannot have `separate == TRUE` with continuous treatments in past`")
        }
      }
  }
  TRUE
}


`%has%` <- function(lhs, rhs) any(rhs %in% lhs)


get_path_inds <- function(A, path) {
  path <- unlist(strsplit(path, "_"))
  data_paths <- strsplit(as.character(A), "_")
  out <- matrix(NA, nrow = length(A), ncol = length(path))

  for (j in seq_along(path)) {
    c_paths <- lapply(
      data_paths,
      function(x) paste0(x[1L:j], collapse = "_")
    )
    out[, j] <- as.numeric(c_paths == paste0(path[1L:j], collapse = "_"))
  }
  out
}


create_blip_list <- function(x, eval_vals, y) {
  J <- length(x)
  N <- nrow(x[[1L]]$outreg_pred)
  n_cols <- lapply(x, function(h) ncol(h$outreg_pred))
  blipped_y <- vector("list", length = J)
  for (j in seq_along(blipped_y)) {
    if (j == J) {
      blipped_y[[j]] <- matrix(y, nrow = N, ncol = 1,
                               dimnames = list(NULL, ""))
    } else {
      fut_inds <- (j + 1):J
      fut <- apply(expand.grid(eval_vals[fut_inds]), 1, paste0, collapse = "_")
      blipped_y[[j]] <- matrix(NA, nrow = N, ncol = length(fut),
                               dimnames = list(NULL, fut))
    }
  }
  blipped_y
}


subset_history_string <- function(x, inds) {
  if (length(x) > 1) stop("only length 1")
  out <- strsplit(x, "_")[[1L]]
  out <- paste0(out[inds], collapse = "_")
  out
}

## from a list of values, create a 
create_history_strings <- function(x, inds) {
  grid <- expand.grid(x[inds])
  out <- apply(grid, 1, paste0, collapse = "_")
  out
}

create_history_factor <- function(A, inds) {
  A <- A[, inds, drop = FALSE]
  int <- interaction(A, sep = "_")
  int
}



#' @importFrom rlang `%||%`
get_eval_vals <- function(object, data) {
  J <- length(object$model_spec)
  out <- vector("list", J)
  for (j in seq_len(J)) {
    this_spec <- object$model_spec[[j]]
    tr_name <- rlang::get_expr(this_spec$treat)
    obs_vals <- sort(unique(data[[tr_name]]))
    out[[j]] <- this_spec$eval_vals %||% obs_vals

    ## don't let eval_vals outside the observed values if there is a
    ## treat_spec specified. IPW doesn't work in this case
    if (!is.null(this_spec$treat_spec)) {
      if (!all(out[[j]] %in% obs_vals)) {
        rlang::abort("`eval_vals` must be among observed values when `treat_spec` specified")
      }
    }
  }
  out
}

check_fold_size <- function(object, data, folds) {

  n_folds <- length(folds)
  fold_balanced <- rep(NA, n_folds)
  for (k in seq_len(n_folds)) {
    fit_rows <- unlist(folds[-k])
    pred_rows <- unlist(folds[k])
    A_fit <- get_treat_df(object, data[fit_rows, ])
    A <- get_treat_df(object, data)
    paths <- interaction(A, sep = "_")
    fit_paths <- interaction(A_fit, sep = "_")  
    num_paths <- length(unique(paths))
    num_fit_paths <- length(unique(fit_paths))
    count_fit_paths <- table(fit_paths)
    fold_balanced[k] <- (num_paths == num_fit_paths) & (all(count_fit_paths > 3))
  }
  return(all(fold_balanced))
}
