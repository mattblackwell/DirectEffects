#' @export
set_treatment <- function(object, treat, formula = NULL) {
  
  model_spec <- object$model_spec
  last_block <- length(model_spec)
  model_spec[[last_block + 1]] <- list(
    treat = rlang::enquo(treat),
    formula = formula
  )
  model_spec <- add_class(model_spec, "model_spec")
  
  new_cde_estimator(
    type = object$type,
    args = object$args,
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
        expr(. ~ . + !!last_form[[2L]] + !!last_form[[3L]])
      )
    } else {
      formula <- update.formula(
        formula,
        expr(. ~ .  + !!last_form[[3L]])
      )
    }
    
  }
  
  this_treat <- list(
    formula = formula,
    engine = engine,
    separate = separate,
    include_past = include_past,
    args = rlang::enquos(...)
  )

  if (length(model_spec) == 0) {
    rlang::abort("`treat_model` called without `set_treat`.")
  }

  model_spec[[last_block]]$treat_spec <- this_treat

  new_cde_estimator(
    type = object$type,
    args = object$args,
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
        expr(`.de_y` ~ . + !!last_tr + !!last_form[[3L]])
      )
    } else {
      formula <- update.formula(
        formula,
        expr(`.de_y` ~ .  + !!last_form[[3L]])
      )
    }    
  } else {
    formula <- update.formula(formula, `.de_y` ~ .)
  }

  if (separate & tr_name %in% all.vars(formula)) {
    rlang::abort("`outreg_model` cannot contain treatment variable if `separate == TRUE`.")
  }
  
  this_outreg <- list(
    formula = formula,
    engine = engine,
    separate = separate,
    include_past = include_past,
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
    model_spec = model_spec
  )
}
