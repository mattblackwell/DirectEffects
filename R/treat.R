reg_engines <- c("lm", "lasso", "rlasso", "ranger_reg")
class_engines <- c("logit", "multinom", "lasso_logit", "rlasso_logit", "lasso_multinom", "ranger_class")
match_engines <- c("Matching")

#' Specifiy a treatment variable for a controlled direct effect
#'
#' This function specifies a treatment variable in the sequence of
#' treatment variables that define the controlled direct effect of
#' interest. 
#' 
#' 
#' @param object A `cde_estimator` object that may or may have
#' previous treatment variables specified/
#' @param treat Name of the treatment variable (not quoted).
#' @param formula One-sided formula giving the covariates that are
#' pre-treatment to this treatment, but post-treatment to any previous
#' treatment. Unless overridden by the arguments to
#' [DirectEffects::treat_model()] or [DirectEffects::outreg_model()],
#' this formula will be the specification used in the modeling of the
#' propensity scores or outcome regressions. 
#' @param treat_type A string indicating the type of variable this is.
#' Takes either the values `"categorical"` or `"regression"` (the
#' latter is not yet implemented). 
#' of 
#' @param eval_vals A numeric vector of values of this variable to
#' evaluate the controlled direct effecct. If `NULL` (the default),
#' this will be set to all observed values of the variable. 
#' @return An updated `cde_estimator` with this information about the
#' treatment specified. 
#' @author Matthew Blackwell
#' @export
#' @md
set_treatment <- function(object,
                          treat,
                          formula = NULL,
                          treat_type = "categorical",
                          eval_vals = NULL) {
  
  model_spec <- object$model_spec
  last_block <- length(model_spec)
  model_spec[[last_block + 1]] <- list(
    treat = rlang::ensym(treat),
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


#' Specify the propensity score model for a CDE treatment
#'
#' Specifies the functional form and estimation engine for a treatment
#' previously specified by [DirectEffects::set_treatment()]. 
#'
#' 
#' @param object A `cde_estimator` object that contains output from a
#' previous call to [DirectEffects::set_treatment()].
#' @param formula A formula specifying the design matrix of the
#' covariates. Passed to fitting engine or used with
#' [stats::model.frame()] and [stats::model.matrix()] to create the
#' design matrix for fitting engines that do not take formulas. 
#' @param engine String indicating the name of the fitting engine. 
#' @param separate Logical indicating whether the fitting algorithm
#' should be applied separately to each history of the treatment
#' variables up to this point (default) or not. 
#' @param include_past A logical value where `TRUE` indicates that
#' formulas passed to previous `treat_model` calls should be appended
#' to the formula given. 
#' @param ... Other arguments to be passed to the engine algorithms.
#' @author Matthew Blackwell
#' @export
#' @md
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
    args$L <- rlang::quo(3)
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
  } else {
    msg <- sprintf("engine type `%s` not supported", engine)
    rlang::abort(message)
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

#' Specify the outcome regression model for a CDE treatment
#'
#' Specifies the functional form and estimation engine for an outcome
#' regression of a treatment previously specified by
#' [DirectEffects::set_treatment()] and the past history of
#' covariates.
#'
#' @inheritParams treat_model
#' @export
#' @md
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
