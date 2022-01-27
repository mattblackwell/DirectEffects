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

empty_est_tab <- function() {
  data.frame(
    term = character(0),
    block_num = integer(0),
    active = character(0),
    control = character(0),
    estimate = numeric(0),
    std_err = numeric(0)
  )
}

check_cde_estimator <- function(...) TRUE


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
