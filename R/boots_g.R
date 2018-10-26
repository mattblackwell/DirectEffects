
#' Bootstrap standard errors
#'
#' Called internally in sequential_g
#'
#' @param seqg Output from sequential_g
#' @param boots The number of bootstrap replicates. Defaults to 1000.
#' 
#' @examples 
#' 
#' data(ploughs)
#' form <- women_politics ~ plow +
#'  agricultural_suitability + tropical_climate + large_animals + rugged | 
#'  years_civil_conflict + years_interstate_conflict  + oil_pc +
#'  european_descent + communist_dummy + polity2_2000 | 
#'  centered_ln_inc + centered_ln_incsq
#' s1 <- sequential_g(form, ploughs)
#'  
#' out.boots <- boots_g(s1)
#' 
#' out.boots
#' 
#' @export
#'
boots_g <- function(seqg, boots = 1000) {
  
  acde.boots <- foreach(b = 1:boots, .combine = "rbind") %do% {
    draw <- sample(1:nrow(seqg$model), replace = TRUE) 
    seqg.draw <- sequential_g(seqg$formula, seqg$model[draw, ])
    coef(seqg.draw)[1:2] # intercept and treatment
  }
  
  out <- list(
    acde.mean = apply(acde.boots, 2, function(x) mean(x, na.rm = TRUE)),
    acde.sd = apply(acde.boots, 2, function(x) sd(x, na.rm = TRUE))
  )
  return(out)
}
