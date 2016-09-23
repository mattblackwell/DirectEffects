sequential.g <- function(formula, first.mod, data, subset, na.action, weights, offset, contrasts = NULL, model = TRUE, y = TRUE, x = FALSE, ...) {
  cl <- match.call()
  if (!identical(first.mod$call$data, cl$data)) {
    stop("data must be the same for both models")
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  formula <- as.Formula(formula)
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  if (inherits(try(terms(formula), silent = TRUE), "try-error")) {
    stop("cannot use dot '.' in formulas")
  }
  ## ensure that the blip formula doesn't include a constant
  f1 <- formula(formula, rhs = 1)
  f2 <- formula(formula, lhs = 0, rhs = 2)
  f2 <- update(f2, ~ . - 1)
  formula <- as.Formula(f1, f2)
  fcoefs <- coef(first.mod)
  bt <- terms(formula, data = data, rhs = 2)
  bnames <- attr(bt, "term.labels")
  if (!all(bnames %in% names(fcoefs))) {
    stop("blip.form contains terms not in first.mod")
  }
  bvars <- match(bnames, names(fcoefs), 0L)
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mtX <- terms(formula, data = data, rhs = 1)
  Y <- model.response(mf, "numeric")
  M <- model.matrix(bt, mf, contrasts)
  Y <- Y - M %*% fcoefs[bvars]
  w <- as.vector(model.weights(mf))
  offset <- as.vector(model.offset(mf))
  X <- model.matrix(mtX, mf, contrasts)
  if (is.null(w)) {
    out <- lm.fit(X, Y, offset = offset, ...)
  } else {
    out <- lm.wfit(X, Y, w, offset = offset, ...)
  }
  out$terms <- list(direct = terms(formula), blip = bt, seqg = mt)
  out$formula <- formula
  out$call <- cl
  out$na.action <- attr(mf, "na.action")
  out$levels <- stats::.getXlevels(mt, mf)
  out$contrasts <- attr(X, "contrasts")
  out$first.mod <- first.mod
  X1 <- model.matrix(first.mod)
  n <- NROW(X)
  Fhat <- crossprod(X, X1)/n
  Fhat[, !(colnames(X1) %in% bnames)] <- 0
  p1 <- 1L:first.mod$rank
  R1 <- chol2inv(first.mod$qr$qr[p1, p1, drop = FALSE])
  efun <- if (is.null(w)) X * out$residual else w * X * out$residual
  w1 <- first.mod$weights
  res1 <- residuals(first.mod)
  efun1 <- if (is.null(w1)) X1 * res1 else w1 * X1 * res1
  ghat <- t(efun) + Fhat %*% R1 %*% t(efun1)
  meat <- crossprod(t(ghat))
  p <- 1L:out$rank
  bread <- chol2inv(out$qr$qr[p, p, drop = FALSE])
  dfc <- (n/(n-max(p)))
  out$vcov <- dfc * (bread %*% meat %*% bread)

  if (model)
    out$model <- mf
  if (x)
    out$x <- X
  if (y)
    out$y <- Y
  class(out) <- "seqg"
  return(out)
}
