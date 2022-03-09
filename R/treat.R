reg_engines <- c("lm", "lasso")
class_engines <- c("logit", "multinom")
match_engines <- c("Matching")

#' @export
set_treatment <- function(object,
                          treat,
                          formula = NULL,
                          treat_type = "categorical",
                          eval_vals = NULL) {
  
  model_spec <- object$model_spec
  last_block <- length(model_spec)
  model_spec[[last_block + 1]] <- list(
    treat = rlang::enquo(treat),
    formula = formula,
    treat_type = treat_type,
    eval_vals = eval_vals
  )
  model_spec <- add_class(model_spec, "model_spec")
  
  new_cde_estimator(
    type = object$type,
    args = object$args,
    formula = object$formula,    
    model_spec = model_spec
  )
  
}
#' @export
treat_model <- function(object,
                      formula,
                      engine,
                      separate = TRUE,
                      include_past = TRUE,
                      ...) {
  
  model_spec <- object$model_spec
  last_block <- length(model_spec)

  

  if (missing(formula)) {
    if (is.null(model_spec[[last_block]]$formula)) {
      rlang::abort("`formula` missing in `treat_model` without being set in `set_treatment`.")
    }
    form <- model_spec[[last_block]]$formula
    tr_name <- rlang::get_expr(model_spec[[last_block]]$treat)
    formula <- update.formula(form, rlang::expr(!!tr_name ~ .))
  }

  if (include_past & last_block > 1L) {
    if (is.null(model_spec[[last_block - 1]]$treat_spec)) {
      rlang::abort("`include_past == TRUE` requires past model specifications.")
    }
    last_form <- model_spec[[last_block - 1]]$treat_spec$formula

    if (!separate) {
      formula <- update.formula(
        formula,
        rlang::expr(. ~ . + !!last_form[[2L]] + !!last_form[[3L]])
      )
    } else {
      formula <- update.formula(
        formula,
        rlang::expr(. ~ .  + !!last_form[[3L]])
      )
    }
    if (is.null(object$formula)) {
      object$formula <- as.formula(rlang::expr( ~ !!formula[[2L]] + !!formula[[3L]]))
    } else {
      object$formula <- update.formula(
        object$formula,
        rlang::expr(. ~ . + !!formula[[2L]] + !!formula[[3L]])
      )      
    }
    
  }
  
  args <- rlang::enquos(...)
  if (engine == "Matching" & is.null(args$L)) {
    args$L <- quo(3)
  }

  ## override treat_type
  if (engine %in% class_engines) {
    model_spec[[last_block]]$treat_type <- "categorical"
    engine_type <- "class"
  } else if (engine %in% reg_engines) {
    model_spec[[last_block]]$treat_type <- "continuous"
    engine_type <- "reg"
  } else if (engine %in% match_engines) {
    model_spec[[last_block]]$treat_type <- "categorical"
    engine_type <- "match"
  }
  
  this_treat <- list(
    formula = formula,
    engine = engine,
    engine_type = engine_type,
    separate = separate,
    include_past = include_past,
    args = args 
  )

  if (length(model_spec) == 0) {
    rlang::abort("`treat_model` called without `set_treat`.")
  }

  model_spec[[last_block]]$treat_spec <- this_treat

  new_cde_estimator(
    type = object$type,
    args = object$args,
    formula = object$formula,
    model_spec = model_spec
  )

}

## set treatment should take on different arguments based on
## estimation strategy. 

#' @export
outreg_model <- function(object,
                         formula,
                         engine,
                         separate = TRUE,
                         include_past = TRUE,
                         ...) {

  model_spec <- object$model_spec
  last_block <- length(model_spec)

  if (length(model_spec) == 0) {
    rlang::abort("`outreg_model` called without `set_treat`.")
  }
  
  tr_name <- rlang::as_label(model_spec[[last_block]]$treat)
  
  if (missing(formula)) {
    if (is.null(model_spec[[last_block]]$formula)) {
      rlang::abort("`formula` missing in `outreg_model` without being set in `set_treatment`.")
    }
    formula <- model_spec[[last_block]]$formula
  }

  if (include_past & last_block > 1L) {
    if (is.null(model_spec[[last_block - 1]]$outreg_spec)) {
      rlang::abort("`include_past == TRUE` requires past model specifications.")
    }
    last_form <- model_spec[[last_block - 1]]$outreg_spec$formula
    if (!separate) {
      last_tr <- rlang::get_expr(model_spec[[last_block - 1]]$treat)
      formula <- update.formula(
        formula,
        rlang::expr(`.de_y` ~ . + !!last_tr + !!last_form[[3L]])
      )
    } else {
      formula <- update.formula(
        formula,
        rlang::expr(`.de_y` ~ .  + !!last_form[[3L]])
      )
    }    
  } else {
    formula <- update.formula(formula, `.de_y` ~ .)
  }
  if (is.null(object$formula)) {
      object$formula <- as.formula(rlang::expr( ~  !!formula[[3L]]))
  } else {
    object$formula <- update.formula(
      object$formula,
      rlang::expr(. ~ . + !!formula[[3L]])
    )      
  }

  if (separate & tr_name %in% all.vars(formula)) {
    rlang::abort("`outreg_model` cannot contain treatment variable if `separate == TRUE`.")
  }
  if (engine %in% reg_engines) {
    engine_type <- "reg"
  } else {
    engine_type <- "class"
  }  
  this_outreg <- list(
    formula = formula,
    engine = engine,
    engine_type = engine_type,
    separate = separate,
    include_past = include_past,
    df = length(attr(terms(formula), "term.labels")),
    args = rlang::enquos(...)
  )

  if (is.null(model_spec)) {
    model_spec <- list()
    model_spec[[1]]$outreg_spec <- this_outreg
  } else {
    ## if current block has an outreg_spec, create new block
    ## but how do we know treatment variable with just outreg
    ## specifications? 
    
    if (names(model_spec[[last_block]])[2] == "outreg_spec") {
      model_spec[[last_block + 1]]$outreg_spec <- this_outreg
    } else {
      model_spec[[last_block]]$outreg_spec <- this_outreg
    }
  }

  new_cde_estimator(
    type = object$type,
    args = object$args,
    formula = object$formula,
    model_spec = model_spec
  )
}
