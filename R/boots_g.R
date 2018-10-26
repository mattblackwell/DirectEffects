
#' Bootstrap standard errors
#'
#' Called internally in sequential_g
#'
#' @param seqg Output from sequential_g
#' @param boots The number of bootstrap replicates. Defaults to 1000.
#' 
#' @export
#'
boots_g <- function(seqg, boots = 1000) {
  
  acde.boots <- foreach(b = 1:boots, .combine = "rbind") %do% {
    draw <- sample(1:nrow(seqg$model), replace = TRUE) 
    seqg.draw <- sequential_g(seqg$formula, seqg$model[draw, ])
    coef(seqg.draw)
  }
  
  out <- list(
    acde.mean = apply(acde.boots, 2, function(x) mean(x, na.rm = TRUE)),
    acde.sd = apply(acde.boots, 2, function(x) sd(x, na.rm = TRUE))
  )
  return(out)
}
